classdef log_utils
    % LOG_UTILS Utility class for managing per-script log files.

    methods (Static)

        function fid = start(logFile)
            % Initializes the log file at the given full file path and stores the active file id.
            if nargin < 1 || isempty(logFile)
                error('log_utils:MissingLogFile', 'A full log file path must be provided.');
            end

            % Ensure parent directory exists
            logDir = fileparts(logFile);
            if ~exist(logDir, 'dir')
                mkdir(logDir);
            end

            % Close any existing active log
            old_fid = log_utils.active_fid_manager('get');
            if ~isempty(old_fid) && old_fid > 0
                fclose(old_fid);
            end

            % Open the new log file
            fid = fopen(logFile, 'w');
            if fid == -1
                error('log_utils:OpenFailed', 'Failed to open log file: %s', logFile);
            end

            log_utils.active_fid_manager('set', fid);
        end

        function printf(varargin)
            % Writes formatted output to the active log file and to the MATLAB command window.
            fid = log_utils.active_fid_manager('get');
            if isempty(fid) || fid < 0
                error('log_utils:NoActiveLog', 'Log not initialized; call log_utils.start first.');
            end
            if nargin < 1
                return;
            end

            fmt = varargin{1};
            args = varargin(2:end);

            % Write to log file.
            fprintf(fid, fmt, args{:});

            % Write to command window.
            fprintf(fmt, args{:});
        end

        function close()
            % Closes the active log file and clears the stored file id.
            fid = log_utils.active_fid_manager('get');
            if ~isempty(fid) && fid > 0
                fclose(fid);
            end
            log_utils.active_fid_manager('clear');
        end

    end

    methods (Static, Access = private)

        function fid = active_fid_manager(action, fid_in)
            % Manages the persistent active file id for the logger.
            persistent active_fid

            switch action
                case 'set'
                    active_fid = fid_in;
                    fid = active_fid;
                case 'get'
                    if isempty(active_fid)
                        fid = [];
                    else
                        fid = active_fid;
                    end
                case 'clear'
                    active_fid = [];
                    fid = [];
                otherwise
                    error('log_utils:InvalidAction', 'Invalid action for active_fid_manager: %s', action);
            end
        end

    end

end
