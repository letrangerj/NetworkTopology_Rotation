% MATLAB script to load, inspect, and export .mat files
% Strategy:
% - syn.mat: Convert to HDF5 (too large for v7 format)
% - rec.mat, pairing_result.mat: Save as MATLAB v7 format for R.matlab package
% - Syn_PWA_inGroup.mat, Rec_MSA_inGroup.mat: Convert to HDF5 using proven method
% Safe ASCII-only version. Run with:
%   matlab -batch "run('scripts/phase_01_explore_validate/inspect_export_mat.m')"

% Compute project root relative to this script when run from
% scripts/phase_01_explore_validate. Expected usage:
%   cd scripts/phase_01_explore_validate
%   matlab -batch "inspect_export_mat"
% Using two levels up from pwd to reach the project root.
project_root = fullfile(pwd, '..', '..');
raw_dir = fullfile(project_root, 'data', 'raw', 'Burkholderiaceae');
interim_dir = fullfile(project_root, 'data', 'interim');

% Ensure interim directory exists
if ~exist(interim_dir, 'dir')
    mkdir(interim_dir);
end

% List of .mat files to process
mat_files = {'pairing_result.mat', 'Syn_PWA_inGroup.mat', 'Rec_MSA_inGroup.mat'};

% Files to save as MATLAB v7 format
v7_files = {'pairing_result.mat'};

% Files to convert to HDF5 format
h5_files = {'Syn_PWA_inGroup.mat', 'Rec_MSA_inGroup.mat'};

% Iterate through files
for i = 1:numel(mat_files)
    mat_file = mat_files{i};
    mat_path = fullfile(raw_dir, mat_file);

    if ~isfile(mat_path)
        fprintf('File not found: %s. Skipping.\n', mat_path);
        continue;
    end

    fprintf('--- Processing %s ---\n', mat_file);

    % Load inside try/catch to capture load errors
    try
        S = load(mat_path);
    catch ME
        fprintf('Error loading %s: %s\n', mat_file, ME.message);
        continue;
    end

    % Inspect variables loaded
    vars = fieldnames(S);
    fprintf('Variables in %s: %s\n', mat_file, strjoin(vars', ', '));
    for v = 1:numel(vars)
        varname = vars{v};
        val = S.(varname);
        try
            fprintf('  %s: class=%s, size=%s\n', varname, class(val), mat2str(size(val)));
        catch
            fprintf('  %s: class=%s, size=%s\n', varname, class(val), mat2str(size(val)));
        end
    end

    % File-specific checks and summaries

    if isfield(S, 'final_pairs')
        try
            fp = S.final_pairs;
            fprintf('final_pairs size=%s, num_pairs=%d\n', mat2str(size(fp)), size(fp,1));
        catch
            fprintf('  Unable to inspect final_pairs\n');
        end
    end

    if isfield(S, 'cluster_info')
        try
            ci = S.cluster_info;
            fprintf('cluster_info fields: %s\n', strjoin(fieldnames(ci)', ', '));
            if isfield(ci, 'distm')
                fprintf('cluster_info.distm size=%s\n', mat2str(size(ci.distm)));
            end
            if isfield(ci, 'dataname')
                fprintf('cluster_info.dataname length=%d\n', numel(ci.dataname));
            end
        catch
            fprintf('  Unable to inspect cluster_info in %s\n', mat_file);
        end
    end

    % Export based on file type
    if ismember(mat_file, v7_files)
        % Save as MATLAB v7 format for R.matlab package
        output_file = strrep(mat_file, '.mat', '_v7.mat');
        output_path = fullfile(interim_dir, output_file);
        try
            % Remove existing file
            if isfile(output_path)
                delete(output_path);
            end

            % Save as MATLAB v7 format
            save(output_path, '-struct', 'S', '-v7');
            fprintf('Exported MATLAB v7: %s\n', output_path);
        catch ME
            fprintf('Failed to export %s to MATLAB v7: %s\n', mat_file, ME.message);
        end

    elseif ismember(mat_file, h5_files)
        % Convert to HDF5
        h5_file = strrep(mat_file, '.mat', '.h5');
        h5_path = fullfile(interim_dir, h5_file);
        try
            % Remove existing HDF5 file
            if isfile(h5_path)
                delete(h5_path);
            end

            if contains(mat_file, 'syn.mat')
                % Handle syn.mat - export struct fields directly
                if isfield(S, 'syn')
                    syn = S.syn;
                    fields = fieldnames(syn);
                    for j = 1:length(fields)
                        fieldName = fields{j};
                        fieldData = syn.(fieldName);
                        datasetPath = ['/syn/' fieldName];

                        if iscell(fieldData)
                            % Handle cell arrays - convert to strings then uint8
                            if all(cellfun(@(x) ischar(x) || isstring(x), fieldData))
                                % String-only cells
                                stringData = string(fieldData);
                                uint8Data = uint8(char(stringData));
                                h5create(h5_path, datasetPath, size(uint8Data), 'Datatype', 'uint8');
                                h5write(h5_path, datasetPath, uint8Data);
                            else
                                % Mixed cells - save as MAT blob
                                tmp = getByteStreamFromArray(fieldData);
                                blobPath = [datasetPath '_matblob'];
                                h5create(h5_path, blobPath, size(tmp), 'Datatype', 'uint8');
                                h5write(h5_path, blobPath, tmp);
                                fprintf('    Saved %s as MAT blob\n', fieldName);
                            end
                        elseif ischar(fieldData)
                            % Handle strings: convert to uint8
                            uint8Data = uint8(fieldData);
                            h5create(h5_path, datasetPath, size(uint8Data), 'Datatype', 'uint8');
                            h5write(h5_path, datasetPath, uint8Data);
                        else
                            % Handle numeric/other data directly
                            h5create(h5_path, datasetPath, size(fieldData));
                            h5write(h5_path, datasetPath, fieldData);
                        end
                        fprintf('    Exported syn.%s\n', fieldName);
                    end
                    fprintf('Exported syn.mat to HDF5: %s\n', h5_path);
                else
                    fprintf('No syn struct found in syn.mat, skipping HDF5 conversion\n', mat_file);
                end

            elseif contains(mat_file, 'rec.mat')
                % Handle rec.mat - export struct fields directly
                if isfield(S, 'rec')
                    rec = S.rec;
                    fields = fieldnames(rec);
                    for j = 1:length(fields)
                        fieldName = fields{j};
                        fieldData = rec.(fieldName);
                        datasetPath = ['/rec/' fieldName];

                        if iscell(fieldData)
                            % Handle cell arrays - convert to strings then uint8
                            if all(cellfun(@(x) ischar(x) || isstring(x), fieldData))
                                % String-only cells
                                stringData = string(fieldData);
                                uint8Data = uint8(char(stringData));
                                h5create(h5_path, datasetPath, size(uint8Data), 'Datatype', 'uint8');
                                h5write(h5_path, datasetPath, uint8Data);
                            else
                                % Mixed cells - save as MAT blob
                                tmp = getByteStreamFromArray(fieldData);
                                blobPath = [datasetPath '_matblob'];
                                h5create(h5_path, blobPath, size(tmp), 'Datatype', 'uint8');
                                h5write(h5_path, blobPath, tmp);
                                fprintf('    Saved %s as MAT blob\n', fieldName);
                            end
                        elseif ischar(fieldData)
                            % Handle strings: convert to uint8
                            uint8Data = uint8(fieldData);
                            h5create(h5_path, datasetPath, size(uint8Data), 'Datatype', 'uint8');
                            h5write(h5_path, datasetPath, uint8Data);
                        else
                            % Handle numeric/other data directly
                            h5create(h5_path, datasetPath, size(fieldData));
                            h5write(h5_path, datasetPath, fieldData);
                        end
                        fprintf('    Exported rec.%s\n', fieldName);
                    end
                    fprintf('Exported rec.mat to HDF5: %s\n', h5_path);
                else
                    fprintf('No rec struct found in rec.mat, skipping HDF5 conversion\n', mat_file);
                end

            else
                % Process cluster_info structure (for PWA/MSA files)
                if isfield(S, 'cluster_info')
                    cluster_info = S.cluster_info;

                    % Convert dataname using proven method: cellstr -> string -> char -> uint8
                    if isfield(cluster_info, 'dataname')
                        dataname = string(cluster_info.dataname);
                        uint8_array = uint8(char(dataname));

                        h5create(h5_path, '/dataname', size(uint8_array), 'Datatype', 'uint8');
                        h5write(h5_path, '/dataname', uint8_array);
                        fprintf('  Exported dataname as uint8 array\n');
                    end

                    % Save distm directly
                    if isfield(cluster_info, 'distm')
                        h5create(h5_path, '/distm', size(cluster_info.distm));
                        h5write(h5_path, '/distm', cluster_info.distm);
                        fprintf('  Exported distm as numeric array\n');
                    end

                    fprintf('Exported HDF5: %s\n', h5_path);
                else
                    fprintf('No cluster_info found in %s, skipping HDF5 conversion\n', mat_file);
                end
            end
        catch ME
            fprintf('Failed to export %s to HDF5: %s\n', mat_file, ME.message);
        end
    end

    fprintf('Finished %s\n\n', mat_file);
end

fprintf('All files processed.\n');
fprintf('MATLAB v7 files saved for: %s\n', strjoin(v7_files, ', '));
fprintf('HDF5 files saved for: %s\n', strjoin(h5_files, ', '));
fprintf('Note: All files exported to HDF5 format for consistency\n');
