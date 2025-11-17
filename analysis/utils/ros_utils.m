% ROS_UTILS: Utility functions for parsing Omniroute ROS bag data
%
% Message Structure Reference for EASE Topics:
%
%   All EASE messages contain 8 × int16 registers: INT0–INT7
%   These are treated as a 16-byte array in little-endian order.
%
%   Byte layout:
%     Bytes 0–1   → uint16  MsgID
%     Bytes 2–3   → uint8   MsgType, uint8 ErrType
%     Bytes 4–5   → uint8   ArgLen,  uint8 Arg0
%     Bytes 6–13  → uint8   Args[1–8] (up to 9 total args)
%     Bytes 14–15 → uint8   Footer
%
%   - The 'RawINT' field stores the original 8×int16 array
%   - The 'Args' field contains the extracted argument bytes (based on ArgLen)
%
% Specific MsgTypes:
%   - 6  → WALL_MOVE: Args = 8 uint8 values, one per chamber.
%         Each byte encodes 8 wall bits (1 = up, 0 = down), wall1 = LSB
%
%   - 14 → GANTRY_MOVE_REL and
%   - 17 → GANTRY_MOVE_ABS: Args = 6 bytes packed as 3 × int16 (X,Y,Z), little-endian
%
% Functions in this library assume:
%   - Wall move commands always include 8 full chamber bytes
%   - Gantry move commands always include exactly 6 bytes

classdef ros_utils
    % Static utility class for ROS bag parsing and inspection
    methods(Static)

        %------------------------------------------------------------------
        function save_table(save_path, tbl)
            % SAVE_TABLE Save a table to both .mat and .csv formats.

            % Validate input
            if ~istable(tbl)
                warning('Input must be a table.');
            end
            if ~ischar(save_path) && ~isstring(save_path)
                error('Save path must be a character vector or string.');
            end

            % Save .mat
            mat_file = strcat(save_path, '.mat');
            save(mat_file, 'tbl');

            % Save .csv
            csv_file = strcat(save_path, '.csv');
            try
                writetable(tbl, csv_file);
            catch ME
                warning('Failed to write CSV: %s', ME.message);
            end
        end

        %------------------------------------------------------------------
        function topics = list_topics(bag)
            % Returns a cell array of topic names in the bag
            info = bag.AvailableTopics;
            topics = info.Properties.RowNames;
        end

        %------------------------------------------------------------------
        function info = get_topic_info(bag, topic_name)
            % Returns the metadata for a single topic
            try
                info = bag.AvailableTopics(topic_name, :);
            catch
                error('Topic "%s" not found in bag.', topic_name);
            end
        end

        %------------------------------------------------------------------
        function tf = check_topic(bag, topic)
            % CHECK_TOPIC Returns true if the topic exists in the bag
            %   Prints a warning if topic is not found

            topics = ros_utils.list_topics(bag);
            tf = any(strcmp(topic, topics));

            if ~tf
                warning('Topic "%s" not found in bag.', topic);
            end
        end

        %------------------------------------------------------------------
        function filtered_bag = filter_topics(bag, topic_list)
            % Returns a rosbag object filtered to only the specified topics
            if ischar(topic_list)
                topic_list = {topic_list};
            end
            filtered_bag = select(bag, 'Topic', topic_list);
        end

        %------------------------------------------------------------------
        function display_summary(bag)
            % Prints summary of all topics and types
            fprintf('Topics in bag file:\n');
            info = bag.AvailableTopics;
            for i = 1:height(info)
                topic_name = info.Properties.RowNames{i};
                msg_type = char(info.MessageType(i));  % convert categorical to char
                num_msgs = info.NumMessages(i);
                fprintf('  %-50s  [%s]  (%d messages)\n', ...
                    topic_name, msg_type, num_msgs);
            end
        end

        %------------------------------------------------------------------
        function event_tbl = parse_event_topic(bag, topic)
            % PARSE_EVENT_TOPIC Parses an omniroute_operation/Event topic
            %   Extracts timestamp and event string from each message.

            if ~ros_utils.check_topic(bag, topic)
                event_tbl = table();
                return;
            end

            sel  = select(bag, 'Topic', topic);
            msgs = readMessages(sel, 'DataFormat', 'struct');
            n    = numel(msgs);

            Time          = zeros(n,1);
            PublisherTime = zeros(n,1);
            Event         = strings(n,1);

            for i = 1:n
                m = msgs{i};
                % Canonical analysis time: bag record/receive time
                Time(i) = double(sel.MessageList.Time(i));

                % Publisher-provided timestamp carried inside the message
                PublisherTime(i) = double(m.Timestamp.Sec) + double(m.Timestamp.Nsec) * 1e-9;

                Event(i) = string(m.Event_);
            end

            event_tbl = table(Time, PublisherTime, Event);
        end

        %------------------------------------------------------------------
        function string_tbl = parse_std_msgs_topic(bag, topic)
            % PARSE_STD_MSGS_TOPIC Parses any std_msgs/String topic

            if ~ros_utils.check_topic(bag, topic)
                string_tbl = table();
                return;
            end

            sel = select(bag, 'Topic', topic);
            msgs = readMessages(sel, 'DataFormat', 'struct');
            n = numel(msgs);

            Time  = zeros(n,1);
            Event = strings(n,1);

            for i = 1:n
                Time(i)  = double(sel.MessageList.Time(i));
                Event(i) = string(msgs{i}.Data);
            end

            string_tbl = table(Time, Event);
        end

        %------------------------------------------------------------------
        function ease_tbl = parse_ease_topic(bag, topic)
            % PARSE_EASE_TOPIC Parses EASE messages from a ROS topic into a structured table.
            % Supports both write and read topics, automatically filtering repeated MsgIDs for reads.
            %
            %   Each message is 8×int16 → 16 bytes, little-endian:
            %     Bytes 0–1   → MsgID   (uint16)
            %     Byte    2   → MsgType (uint8)
            %     Byte    3   → ErrType (uint8)
            %     Byte    4   → ArgLen  (uint8)
            %     Bytes 5–13  → Args[0–8] (uint8)
            %     Bytes 14–15 → Footer  (2×uint8)

            if ~ros_utils.check_topic(bag, topic)
                ease_tbl = table();
                return;
            end

            % Read messages
            sel = select(bag, 'Topic', topic);
            msgs = readMessages(sel, 'DataFormat', 'struct');
            n = numel(msgs);

            % Preallocate
            Time    = zeros(n,1);
            MsgID   = zeros(n,1,'uint16');
            MsgType = zeros(n,1,'uint8');
            ErrType = zeros(n,1,'uint8');
            ArgLen  = zeros(n,1,'uint8');
            Args    = cell(n,1);
            RawINT  = zeros(n,8,'int16');

            for i = 1:n
                m = msgs{i};
                reg_ints = int16([m.INT0, m.INT1, m.INT2, m.INT3, m.INT4, m.INT5, m.INT6, m.INT7]);
                RawINT(i,:) = reg_ints;

                % Convert 8×int16 → 16×uint8 byte stream
                byteStream = zeros(16,1,'uint8');
                for j = 1:8
                    u16 = typecast(reg_ints(j), 'uint16');
                    byteStream((j-1)*2 + (1:2)) = typecast(u16, 'uint8');  % little-endian
                end

                % Extract fields from known byte offsets
                MsgID(i)   = typecast(byteStream(1:2), 'uint16');  % bytes 0–1
                MsgType(i) = byteStream(3);                        % byte 2
                ErrType(i) = byteStream(4);                        % byte 3
                ArgLen(i)  = byteStream(5);                        % byte 4

                if ArgLen(i) > 0 && (5 + ArgLen(i)) <= 16
                    Args{i} = byteStream(6 : 6 + ArgLen(i) - 1);
                else
                    Args{i} = [];
                end

                Time(i) = double(sel.MessageList.Time(i));
            end

            ease_tbl = table(Time, MsgID, MsgType, ErrType, ArgLen, Args, RawINT);

            % Filter repeated MsgIDs in read topics
            if contains(topic, '/Esmacat_read_')
                [~, first_idx] = unique(MsgID, 'stable');  % keep first occurrence of each MsgID
                ease_tbl = ease_tbl(first_idx, :);
            end
        end

        %------------------------------------------------------------------
        function wall_state_tbl = decode_wall_state_table(maze_write_tbl)
            %DECODE_WALL_STATE_TABLE Parses wall move commands from maze_write_tbl
            %   Returns a table with one row per MsgType==6 (wall move) command,
            %   including Time, MsgID, 8 chamber wall states, and wall movement deltas

            % Filter for wall move messages
            isWallMove = maze_write_tbl.MsgType == 6;
            filteredTbl = maze_write_tbl(isWallMove, :);

            nRows = height(filteredTbl);
            Time  = filteredTbl.Time;
            MsgID = filteredTbl.MsgID;
            MsgType = filteredTbl.MsgType;

            % Preallocate cell arrays for each chamber
            Ch = cell(8, 1);
            for ch = 1:8
                Ch{ch} = repmat({false(8,1)}, nRows, 1);
            end

            for i = 1:nRows
                args = filteredTbl.Args{i};
                if length(args) >= 8
                    for ch = 1:8
                        byteVal = args(ch);
                        wallBits = logical(bitget(byteVal, 1:8))';  % LSB to MSB
                        Ch{ch}{i} = wallBits;
                    end
                end
            end

            % Compute wall movement deltas
            WallsUp = zeros(nRows,1);
            WallsDown = zeros(nRows,1);
            prevState = false(8,8);  % 8 chambers × 8 walls

            for i = 1:nRows
                currState = false(8,8);
                for ch = 1:8
                    currState(ch,:) = Ch{ch}{i}';
                end
                upMoves = currState & ~prevState;
                downMoves = ~currState & prevState;
                WallsUp(i) = sum(upMoves(:));
                WallsDown(i) = sum(downMoves(:));
                prevState = currState;
            end

            % Assemble output table
            wall_state_tbl = table(Time, MsgID, MsgType, ...
                Ch{1}, Ch{2}, Ch{3}, Ch{4}, Ch{5}, Ch{6}, Ch{7}, Ch{8}, ...
                WallsUp, WallsDown, ...
                'VariableNames', {'Time', 'MsgID', 'MsgType', ...
                'Ch1', 'Ch2', 'Ch3', 'Ch4', 'Ch5', 'Ch6', 'Ch7', 'Ch8', ...
                'WallsUp', 'WallsDown'});
        end

        %------------------------------------------------------------------
        function pose_tbl = parse_pose_topic(bag, topic)
            % PARSE_POSE_TOPIC Extracts flat pose data from a PoseStamped topic.
            %   Canonical timebase:
            %       Time       = bag receive time (MessageList.Time) [seconds]
            %   Additional stamp:
            %       HeaderTime = publisher-stamped header time (Header.Stamp) [seconds]
            %   Other fields:
            %       X,Y,Z, qX,qY,qZ,qW, FrameId

            if ~ros_utils.check_topic(bag, topic)
                pose_tbl = table();
                return;
            end

            sel  = select(bag, 'Topic', topic);
            msgs = readMessages(sel, 'DataFormat', 'struct');
            n    = numel(msgs);

            % Preallocate
            Time      = zeros(n,1);
            HeaderTime= zeros(n,1);
            X         = zeros(n,1);
            Y         = zeros(n,1);
            Z         = zeros(n,1);
            qX        = zeros(n,1);
            qY        = zeros(n,1);
            qZ        = zeros(n,1);
            qW        = zeros(n,1);
            FrameId   = strings(n,1);

            for i = 1:n
                m = msgs{i};

                % Canonical analysis time: bag record/receive time
                Time(i) = double(sel.MessageList.Time(i));

                % Publisher-provided header stamp
                HeaderTime(i) = double(m.Header.Stamp.Sec) + double(m.Header.Stamp.Nsec) * 1e-9;

                FrameId(i) = string(m.Header.FrameId);

                X(i)  = m.Pose.Position.X;
                Y(i)  = m.Pose.Position.Y;
                Z(i)  = m.Pose.Position.Z;

                qX(i) = m.Pose.Orientation.X;
                qY(i) = m.Pose.Orientation.Y;
                qZ(i) = m.Pose.Orientation.Z;
                qW(i) = m.Pose.Orientation.W;
            end

            pose_tbl = table(Time, HeaderTime, X, Y, Z, qX, qY, qZ, qW, FrameId);
        end

        %------------------------------------------------------------------
        function gantry_move_tbl = decode_gantry_move_table(gantry_write_tbl, gantry_pose_tbl)
            % DECODE_GANTRY_MOVE_TABLE Extracts gantry movement coordinates aligned to pose frame
            %   Applies static offset derived from first row of gantry_pose_tbl.
            %
            %   Inputs:
            %     gantry_write_tbl - table from parse_ease_topic for /Esmacat_write_gantry_ease
            %     gantry_pose_tbl  - table with tracked gantry poses (used to infer static offset)
            %
            %   Output:
            %     gantry_move_tbl - table with columns:
            %       Time      - timestamp of move command
            %       MsgID     - unique message ID
            %       MsgType   - type of movement (14 = REL, 17 = ABS)
            %       X         - target X position in meters (pose-aligned)
            %       Y         - target Y position in meters (pose-aligned)

            % Compute fixed offset from homed pose
            offset_x = gantry_pose_tbl.X(1);
            offset_y = gantry_pose_tbl.Y(1);

            % Select only movement messages (REL or ABS)
            isMove = gantry_write_tbl.MsgType == 14 | gantry_write_tbl.MsgType == 17;
            filteredTbl = gantry_write_tbl(isMove, :);

            % Preallocate
            nRows = height(filteredTbl);
            Time     = filteredTbl.Time;
            MsgID    = filteredTbl.MsgID;
            MsgType  = filteredTbl.MsgType;
            X = NaN(nRows, 1);
            Y = NaN(nRows, 1);

            % Decode and swap X/Y (convert mm → m) and apply offset
            for i = 1:nRows
                args = filteredTbl.Args{i};
                if length(args) == 6
                    % Original EASE encoding: [Y, X, Z] as 3 × int16
                    raw_y = typecast(uint8(args(1:2)), 'int16');  % true Y
                    raw_x = typecast(uint8(args(3:4)), 'int16');  % true X

                    X(i) = double(raw_x) / 1000 + offset_x;  % meters
                    Y(i) = double(raw_y) / 1000 + offset_y;  % meters
                end
            end

            % Construct output table without Z
            gantry_move_tbl = table(Time, MsgID, MsgType, X, Y);
        end

        %------------------------------------------------------------------
        function latency_tbl = match_ease_write_read(write_tbl, read_tbl)
            % MATCH_EASE_WRITE_READ Matches EASE write messages to read confirmations
            %   Computes round-trip latency between write and first matching read.
            %
            %   Inputs:
            %       write_tbl - table from parse_ease_topic (e.g., Esmacat_write_*)
            %       read_tbl  - table from parse_ease_topic (e.g., Esmacat_read_*)
            %
            %   Output:
            %       latency_tbl - table with columns:
            %           TimeWrite, TimeRead, MsgDT, MsgID, MsgType

            if isempty(write_tbl) || isempty(read_tbl)
                latency_tbl = [];
                return
            end

            nWrite = height(write_tbl);

            TimeWrite = write_tbl.Time;
            MsgID     = write_tbl.MsgID;
            MsgType   = write_tbl.MsgType;

            TimeRead  = NaN(nWrite,1);
            MsgDT     = NaN(nWrite,1);

            for i = 1:nWrite
                id  = MsgID(i);
                typ = MsgType(i);

                % Find first matching read (same MsgID and MsgType)
                idx = find(read_tbl.MsgID == id & read_tbl.MsgType == typ, 1, 'first');

                if ~isempty(idx)
                    t_read = read_tbl.Time(idx);
                    TimeRead(i) = t_read;
                    MsgDT(i) = t_read - TimeWrite(i);
                end
            end

            latency_tbl = table(TimeWrite, TimeRead, MsgDT, MsgID, MsgType);
        end

        %------------------------------------------------------------------
        function result_tbl = analyze_gantry_move_v1(move_tbl, pose_tbl)
            % ANALYZE_GANTRY_MOVE Match gantry pose data to move commands to compute latency and accuracy
            %
            % Inputs:
            %   move_tbl - Table with gantry move commands (fields: Time, MsgID, X, Y in meters)
            %   pose_tbl - Table with tracked gantry poses (fields: Time, X, Y in meters)
            %
            % Output:
            %   result_tbl - Table with:
            %       MsgID       : EASE command identifier
            %       TimeWrite   : Timestamp when command was issued
            %       TargetX     : Target X position (meters)
            %       TargetY     : Target Y position (meters)
            %       TimeArrived : Time of closest approach to target
            %       Latency     : TimeArrived - TimeWrite (seconds)
            %       Accuracy    : Final distance to target (meters)

            nMoves = height(move_tbl);
            MsgID = move_tbl.MsgID;
            TimeWrite = move_tbl.Time;

            % Target coordinates already in meters
            TargetX = move_tbl.X;
            TargetY = move_tbl.Y;

            % Preallocate output vectors
            TimeArrived = NaN(nMoves, 1);
            Latency     = NaN(nMoves, 1);
            Accuracy    = NaN(nMoves, 1);

            % Pose data
            pose_times = pose_tbl.Time;
            pose_xy = [pose_tbl.X, pose_tbl.Y];

            for i = 1:nMoves
                % Define time window: current to next move
                t_start = TimeWrite(i);
                if i < nMoves
                    t_end = TimeWrite(i+1);
                else
                    t_end = inf;
                end

                % Get pose samples in this window
                in_window = pose_times >= t_start & pose_times < t_end;
                if ~any(in_window)
                    continue;
                end

                % Compute Euclidean distance to target
                target = [TargetX(i), TargetY(i)];
                pos = pose_xy(in_window, :);
                dists = vecnorm(pos - target, 2, 2);

                % Find closest approach
                [min_dist, idx] = min(dists);
                match_idx = find(in_window, 1, 'first') - 1 + idx;
                TimeArrived(i) = pose_times(match_idx);
                Latency(i)     = TimeArrived(i) - TimeWrite(i);
                Accuracy(i)    = min_dist;
            end

            % Construct output table
            result_tbl = table(MsgID, TimeWrite, TargetX, TargetY, ...
                TimeArrived, Latency, Accuracy);
        end

        function result_tbl = analyze_gantry_move(move_tbl, pose_tbl)
            % ANALYZE_GANTRY_MOVE Match gantry pose data to move commands to compute latency and accuracy
            %
            % Inputs:
            %   move_tbl - Table with gantry move commands (fields: Time, MsgID, X, Y in meters)
            %   pose_tbl - Table with tracked gantry poses (fields: Time, X, Y in meters)
            %
            % Output:
            %   result_tbl - Table with:
            %       MsgID       : EASE command identifier
            %       TimeWrite   : Timestamp when command was issued
            %       TargetX     : Target X position (meters)
            %       TargetY     : Target Y position (meters)
            %       TimeArrived : First time pose is within tolerance of target
            %       Latency     : TimeArrived - TimeWrite (seconds)
            %       Accuracy    : Final distance to target (meters)

            % --- Parameters ---
            arrival_tol = 0.02;  % meters

            nMoves = height(move_tbl);
            MsgID = move_tbl.MsgID;
            TimeWrite = move_tbl.Time;
            TargetX = move_tbl.X;
            TargetY = move_tbl.Y;

            TimeArrived = NaN(nMoves, 1);
            Latency     = NaN(nMoves, 1);
            Accuracy    = NaN(nMoves, 1);

            pose_times = pose_tbl.Time;
            pose_xy = [pose_tbl.X, pose_tbl.Y];

            for i = 1:nMoves
                t_start = TimeWrite(i);
                if i < nMoves
                    t_end = TimeWrite(i+1);
                else
                    t_end = t_start + 6;
                end

                in_window = pose_times >= t_start & pose_times < t_end;
                if ~any(in_window)
                    continue;
                end

                target = [TargetX(i), TargetY(i)];
                pos = pose_xy(in_window, :);
                dists = vecnorm(pos - target, 2, 2);

                % Accuracy = closest approach (same as before)
                [min_dist, ~] = min(dists);
                Accuracy(i) = min_dist;

                % TimeArrived = first within tolerance
                idx_within_tol = find(dists <= arrival_tol, 1, 'first');
                if ~isempty(idx_within_tol)
                    match_idx = find(in_window, 1, 'first') - 1 + idx_within_tol;
                    TimeArrived(i) = pose_times(match_idx);
                    Latency(i) = TimeArrived(i) - TimeWrite(i);
                end
            end

            result_tbl = table(MsgID, TimeWrite, TargetX, TargetY, ...
                TimeArrived, Latency, Accuracy);
        end

        % ----------------------------- PROJECTION UTILS ----------------------------
        function rosout_tbl = parse_rosout(bag, topic)
            % PARSE_ROSOUT Parse rosgraph_msgs/Log from /rosout.
            %   Canonical analysis time is the bag receive time (MessageList.Time).
            %
            % Output table fields:
            %   Time   : double (s) bag record time
            %   Level  : uint8   (1=DEBUG, 2=INFO, 4=WARN, 8=ERROR, 16=FATAL)
            %   Name   : string  node name (e.g., /projection_display)
            %   Msg    : string  log text

            if nargin < 2 || isempty(topic), topic = '/rosout'; end
            if ~ros_utils.check_topic(bag, topic)
                rosout_tbl = table();
                return;
            end

            sel  = select(bag, 'Topic', topic);
            msgs = readMessages(sel, 'DataFormat', 'struct');
            n    = numel(msgs);

            Time  = zeros(n,1);
            Level = zeros(n,1,'uint8');
            Name  = strings(n,1);
            Msg   = strings(n,1);

            for i = 1:n
                m = msgs{i};
                Time(i)  = double(sel.MessageList.Time(i));
                Level(i) = uint8(m.Level);
                Name(i)  = string(m.Name);
                Msg(i)   = string(m.Msg);
            end

            rosout_tbl = table(Time, Level, Name, Msg);
        end

        % --------------------------------------------------------------------------
        function image_cmd_tbl = parse_int32_multiarray(bag, topic)
            % PARSE_INT32_MULTIARRAY Parse std_msgs/Int32MultiArray (e.g., /projection_image).
            %   Stores Data as a row vector per message.
            %
            % Output table fields:
            %   Time   : double (s) bag record time
            %   Data   : cell array, each cell is 1×K int32 row vector
            %   CmdIdx : uint32, monotonically increasing command index (1..N)

            if nargin < 2 || isempty(topic), topic = '/projection_image'; end
            if ~ros_utils.check_topic(bag, topic)
                image_cmd_tbl = table();
                return;
            end

            sel  = select(bag, 'Topic', topic);
            msgs = readMessages(sel, 'DataFormat', 'struct');
            n    = numel(msgs);

            Time   = zeros(n,1);
            Data   = cell(n,1);
            CmdIdx = zeros(n,1,'uint32');

            for i = 1:n
                m = msgs{i};
                Time(i)   = double(sel.MessageList.Time(i));
                % Flatten Data to row; ensure int32
                vec = int32(m.Data(:)).';
                Data{i}   = vec;
                CmdIdx(i) = uint32(i);
            end

            image_cmd_tbl = table(Time, Data, CmdIdx);
        end

        % --------------------------------------------------------------------------
        function proj_latency_tbl = compute_projection_last_render_latency(image_cmd_tbl, rosout_tbl, opts)
            % COMPUTE_PROJECTION_LAST_RENDER_LATENCY
            %   For each /projection_image command, find subsequent /rosout logs
            %   of the form "[projector_image] Timer stopped, Delta: <sec>" and
            %   take the last timestamp in that command's window as completion.
            %
            % Inputs:
            %   image_cmd_tbl : from parse_int32_multiarray (fields: Time, CmdIdx)
            %   rosout_tbl    : from parse_rosout (fields: Time, Msg)
            %   opts          : struct with optional fields:
            %       .timer_name     : default 'projector_image'
            %       .regex_timer    : default '^\[projector_image\]\s*Timer stopped, Delta:\s*([0-9]*\.?[0-9]+)'
            %       .min_per_cmd    : expected number of stops per command (default 1)
            %
            % Output table fields (one row per image command):
            %   TimeCmd            : command time (s)
            %   CmdIdx             : uint32 command index
            %   TimeLast           : completion time = latest matching rosout time in window
            %   LatencyLast_s      : TimeLast - TimeCmd (s)
            %   NStopsInWindow     : number of matching logs found in the window
            %   StopTimes          : cell{1×1} double vector of individual stop times
            %   StopDeltas_s       : cell{1×1} double vector (parsed Delta per stop, if present)
            %
            % Notes:
            %   - Window k is [TimeCmd(k), TimeCmd(k+1)), last command uses [TimeCmd(end), +inf).
            %   - If no stop logs are found in a window, TimeLast and LatencyLast_s are NaN.

            arguments
                image_cmd_tbl table
                rosout_tbl table
                opts.timer_name (1,1) string = "projector_image"
                opts.regex_timer (1,1) string = "^\[projector_image\]\s*Timer stopped, Delta:\s*([0-9]*\.?[0-9]+)"
                opts.min_per_cmd (1,1) double = 1
            end

            if isempty(image_cmd_tbl)
                proj_latency_tbl = table();
                return;
            end

            % Pre-filter rosout to candidate lines to reduce work
            msg = rosout_tbl.Msg;
            has_phrase = contains(msg, "Timer stopped, Delta:");
            % Also require the bracketed timer name (defensive)
            has_name   = contains(msg, "[" + opts.timer_name + "]");
            cand_idx   = find(has_phrase & has_name);

            cand_times = rosout_tbl.Time(cand_idx);
            cand_msgs  = msg(cand_idx);

            % Prepare outputs
            n = height(image_cmd_tbl);
            TimeCmd        = image_cmd_tbl.Time;
            CmdIdx         = image_cmd_tbl.CmdIdx;
            TimeLast       = nan(n,1);
            LatencyLast_s  = nan(n,1);
            NStopsInWindow = zeros(n,1);
            StopTimes      = cell(n,1);
            StopDeltas_s   = cell(n,1);

            % Build window edges
            T = TimeCmd(:);
            T_next = [T(2:end); inf];

            % Regex for Delta capture
            expr = char(opts.regex_timer); % keep given regex verbatim

            for k = 1:n
                t0 = T(k); t1 = T_next(k);
                in_win = cand_times >= t0 & cand_times < t1;

                if ~any(in_win)
                    NStopsInWindow(k) = 0;
                    StopTimes{k}      = [];
                    StopDeltas_s{k}   = [];
                    continue;
                end

                times_k = cand_times(in_win);
                msgs_k  = cand_msgs(in_win);

                % Parse optional Delta
                deltas = nan(size(times_k));
                for j = 1:numel(msgs_k)
                    tok = regexp(msgs_k(j), expr, 'tokens', 'once');
                    if ~isempty(tok)
                        val = str2double(tok{1});
                        if isfinite(val), deltas(j) = val; end
                    end
                end

                % Record
                StopTimes{k}      = double(times_k(:)).';
                StopDeltas_s{k}   = deltas(:).';
                NStopsInWindow(k) = numel(times_k);

                % Completion = latest stop timestamp
                t_last = max(times_k);
                TimeLast(k)      = t_last;
                LatencyLast_s(k) = t_last - t0;
            end

            proj_latency_tbl = table(TimeCmd, CmdIdx, TimeLast, LatencyLast_s, ...
                NStopsInWindow, StopTimes, StopDeltas_s);
        end

        % --------------------------------------------------------------------------
        function audio_latency_tbl = compute_audio_start_latency(sound_cmd_tbl, rosout_tbl)
            % COMPUTE_AUDIO_START_LATENCY Compute latency from audio command to playback start.
            %
            %   audio_latency_tbl = ros_utils.compute_audio_start_latency(sound_cmd_tbl, rosout_tbl)
            %
            %   Inputs:
            %       sound_cmd_tbl - Table from parse_std_msgs_topic for '/sound_cmd'
            %                       Fields: Time (double), Event (string)
            %       rosout_tbl    - Table from parse_rosout for '/rosout'
            %                       Fields: Time (double), Msg (string)
            %
            %   Output:
            %       audio_latency_tbl - Table with:
            %           TimeCmd        : Command issue time (s)
            %           SoundName      : String name of sound
            %           TimeStart      : Playback start time (s)
            %           LatencyStart_s : TimeStart - TimeCmd (s)
            %
            %   Notes:
            %       - Matches the first "[sound_generator] Played sound:" log after each command.
            %       - Requires both topics to be present; returns empty table otherwise.
            %

            if isempty(sound_cmd_tbl) || isempty(rosout_tbl)
                audio_latency_tbl = table([], [], [], [], ...
                    'VariableNames', {'TimeCmd','SoundName','TimeStart','LatencyStart_s'});
                return;
            end

            % Extract command times and names directly (no guessing)
            TimeCmd   = sound_cmd_tbl.Time(:);
            SoundName = string(sound_cmd_tbl.Event(:));

            % Filter rosout for audio start logs
            logTimes = rosout_tbl.Time;
            logMsgs  = rosout_tbl.Msg;
            is_audio = startsWith(logMsgs, "[sound_generator] Played sound:");
            logTimesA = logTimes(is_audio);
            logMsgsA  = logMsgs(is_audio);

            % Preallocate results
            TimeStart = nan(size(TimeCmd));

            % Match each command to first playback log (prefer name match)
            for k = 1:numel(TimeCmd)
                if strlength(SoundName(k)) > 0
                    j = find(logTimesA >= TimeCmd(k) & endsWith(logMsgsA, SoundName(k)), 1, 'first');
                    if isempty(j)
                        j = find(logTimesA >= TimeCmd(k), 1, 'first');
                    end
                else
                    j = find(logTimesA >= TimeCmd(k), 1, 'first');
                end

                if ~isempty(j)
                    TimeStart(k) = logTimesA(j);
                end
            end

            LatencyStart_s = TimeStart - TimeCmd;
            audio_latency_tbl = table(TimeCmd, SoundName, TimeStart, LatencyStart_s);
        end

        % --------------------------------------------------------------------------
        function gate_latency_tbl = compute_gate_move_latency(maze_write_tbl, maze_read_tbl, gate_testing_tbl)
            % COMPUTE_GATE_MOVE_LATENCY Compute round-trip latency for gate (wall) moves.
            %
            %   gate_latency_tbl = ros_utils.compute_gate_move_latency( ...
            %       maze_write_tbl, maze_read_tbl, gate_testing_tbl)
            %
            %   Inputs:
            %       maze_write_tbl   - Table from parse_ease_topic for '/Esmacat_write_maze_ease'
            %                          Must contain MsgType, MsgID, Time columns.
            %       maze_read_tbl    - Table from parse_ease_topic for '/Esmacat_read_maze_ease'
            %                          Must contain MsgType, MsgID, Time columns.
            %       gate_testing_tbl - Table from parse_std_msgs_topic for '/gate_testing'
            %                          Must contain Time (double) and Event (string).
            %
            %   Output:
            %       gate_latency_tbl - Table with:
            %           TimeWrite  : Command issue time (s)
            %           TimeRead   : First matching read confirmation time (s)
            %           Latency_s  : TimeRead - TimeWrite (s)
            %           MsgID      : EASE command identifier
            %           MsgType    : EASE message type (6 = wall move)
            %           MoveDir    : 'up' or 'down'
            %
            %   Notes:
            %       - Filters write commands to MsgType == 6 (wall move).
            %       - Matches by MsgID and MsgType using match_ease_write_read.
            %       - Orders output exactly as 'gate_test_wall_move_up'/'down' events occurred.
            %

            % Validate inputs
            if isempty(maze_write_tbl) || isempty(maze_read_tbl) || isempty(gate_testing_tbl)
                gate_latency_tbl = table([], [], [], [], [], [], ...
                    'VariableNames', {'TimeWrite','TimeRead','Latency_s','MsgID','MsgType','MoveDir'});
                return;
            end

            % Filter to wall move commands
            write_wall_tbl = maze_write_tbl(maze_write_tbl.MsgType == 6, :);
            read_wall_tbl  = maze_read_tbl(maze_read_tbl.MsgType == 6, :);

            % Compute base latencies
            base_tbl = ros_utils.match_ease_write_read(write_wall_tbl, read_wall_tbl);
            if isempty(base_tbl)
                gate_latency_tbl = table([], [], [], [], [], [], ...
                    'VariableNames', {'TimeWrite','TimeRead','Latency_s','MsgID','MsgType','MoveDir'});
                return;
            end

            % Rename MsgDT to Latency_s
            base_tbl.Properties.VariableNames{'MsgDT'} = 'Latency_s';

            % Identify only up/down events from gate_testing_tbl
            is_move_event = startsWith(gate_testing_tbl.Event, 'gate_test_wall_move_');
            move_events   = gate_testing_tbl(is_move_event, :);

            % Extract move direction ('up' or 'down')
            move_dir = erase(move_events.Event, 'gate_test_wall_move_');  % leaves 'up'/'down'

            % For each move event, find the matching write command
            nMoves = height(move_events);
            rows_out = [];
            dirs_out = strings(nMoves,1);

            for k = 1:nMoves
                t_evt  = move_events.Time(k);
                dir_evt = move_dir(k);

                % Match to nearest write command at or after event time
                idx = find(abs(base_tbl.TimeWrite - t_evt) < 0.5, 1, 'first'); % 0.5s tolerance
                if isempty(idx)
                    idx = find(base_tbl.TimeWrite >= t_evt, 1, 'first');
                end

                if ~isempty(idx)
                    rows_out(end+1) = idx; %#ok<AGROW>
                    dirs_out(length(rows_out)) = dir_evt;
                end
            end

            % Build output table in event order
            gate_latency_tbl = base_tbl(rows_out, :);
            gate_latency_tbl.MoveDir = dirs_out(:);
        end


    end
end
