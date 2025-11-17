classdef path_utils
    % PATH_UTILS Utility for resolving project paths, environment values, and shared utils.

    methods (Static)

        function root = repo_root()
            % Returns absolute path to the repository root directory.
            persistent cached_root
            if isempty(cached_root)
                this_file   = mfilename('fullpath');
                this_dir    = fileparts(this_file);      % .../analysis/utils
                analysis_dir = fileparts(this_dir);      % .../analysis
                cached_root  = fileparts(analysis_dir);  % .../omniroute-manuscript-analysis
            end
            root = cached_root;
        end

        function p = analysis()
            % Returns absolute path to the analysis directory.
            p = fullfile(path_utils.repo_root(), 'analysis');
        end

        function p = results()
            % Returns absolute path to the results directory.
            p = fullfile(path_utils.repo_root(), 'results');
        end

        function p = results_figures(subdir)
            % Returns absolute path to results/figures/<subdir> and creates it if needed.
            base = fullfile(path_utils.repo_root(), 'results', 'figures');
            if nargin < 1
                p = base;
                if ~exist(p, 'dir'), mkdir(p); end
                return;
            end
            p = fullfile(base, subdir);
            if ~exist(p, 'dir'), mkdir(p); end
        end

        function p = results_tables(subdir)
            % Returns absolute path to results/tables/<subdir> and creates it if needed.
            base = fullfile(path_utils.repo_root(), 'results', 'tables');
            if nargin < 1
                p = base;
                if ~exist(p, 'dir'), mkdir(p); end
                return;
            end
            p = fullfile(base, subdir);
            if ~exist(p, 'dir'), mkdir(p); end
        end

        function p = results_log(fi_name)
            % Returns absolute path to results/logs/<fi_name>.txt and ensures the log directory exists.
            if nargin < 1 || isempty(fi_name)
                error('path_utils:MissingFileName', 'fi_name must be provided.');
            end

            base = fullfile(path_utils.repo_root(), 'results', 'logs');
            if ~exist(base, 'dir')
                mkdir(base);
            end

            p = fullfile(base, [fi_name, '.txt']);

            % Delete existing log file if present
            if exist(p, 'file')
                delete(p);
            end
        end

        function p = data_repo()
            % Returns absolute path to the repo-local data directory.
            p = fullfile(path_utils.repo_root(), 'data');
        end

        function p = dataset_root()
            % Returns external dataset root defined by dataset_root in .env.
            env = path_utils.env();
            if isfield(env, 'DATASET_ROOT')
                p = env.DATASET_ROOT;
            else
                error('path_utils:MissingEnvKey', 'dataset_root not defined in .env.');
            end
        end

        function env_struct = env()
            % Parses the .env file at the repo root into a struct of key/value pairs.
            persistent cached_env
            if isempty(cached_env)
                env_path = fullfile(path_utils.repo_root(), '.env');
                env_struct = struct();
                if exist(env_path, 'file')
                    txt = fileread(env_path);
                    lines = regexp(txt, '\r\n|\n|\r', 'split');
                    for i = 1:numel(lines)
                        line = strtrim(lines{i});
                        if isempty(line) || startsWith(line, '#')
                            continue;
                        end
                        eq_idx = strfind(line, '=');
                        if isempty(eq_idx)
                            continue;
                        end
                        key   = strtrim(line(1:eq_idx(1)-1));
                        value = strtrim(line(eq_idx(1)+1:end));
                        if isempty(key)
                            continue;
                        end
                        field = matlab.lang.makeValidName(key);
                        env_struct.(field) = value;
                    end
                end
                cached_env = env_struct;
            else
                env_struct = cached_env;
            end
        end

        function clean_results_dirs()
            % CLEAN_RESULTS_DIRS
            % Remove all files and subfolders inside results/{figures,logs,tables}
            % while preserving the top-level .gitkeep file.

            root = path_utils.repo_root();

            dirs = {
                fullfile(root, 'results', 'figures')
                fullfile(root, 'results', 'logs')
                fullfile(root, 'results', 'tables')
                };

            for i = 1:numel(dirs)
                d = dirs{i};
                if ~isfolder(d), continue; end

                items = dir(d);
                for k = 1:numel(items)
                    name = items(k).name;

                    % Skip ".", "..", and ".gitkeep"
                    if any(strcmp(name, {'.', '..', '.gitkeep'}))
                        continue;
                    end

                    fp = fullfile(d, name);

                    if items(k).isdir
                        rmdir(fp, 's');
                    else
                        delete(fp);
                    end
                end
            end
        end

    end

end
