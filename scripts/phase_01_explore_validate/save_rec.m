% Script to save MATLAB struct 'rec' to HDF5 using helper function saveStructToHDF5

% Load the MATLAB struct 'rec' from .mat file

load('/home/wulong/RotationWorks/PhylogenicTree/data/raw/Burkholderiaceae/rec.mat'); % Load the MATLAB struct 'rec'
rec_path = '/home/wulong/RotationWorks/PhylogenicTree/data/interim/rec.h5';
% Save the struct 'rec' to HDF5 format
saveRecToHDF5(rec_path, rec);

% Function to save struct to HDF5 with proper type handling
function saveRecToHDF5(filename, rec)
    % Delete existing file
    if exist(filename, 'file')
        delete(filename);
    end

    fprintf('Saving rec struct to %s...\n', filename);

    % 1. count - scalar value
    h5create(filename, '/rec/count', 1);
    h5write(filename, '/rec/count', rec.count);
    fprintf('  Saved count: %d\n', rec.count);

    % 2. foldername - cell array of strings
    foldername_chars = char(rec.foldername);
    foldername_uint8 = uint8(foldername_chars);
    h5create(filename, '/rec/foldername', size(foldername_uint8), 'Datatype', 'uint8');
    h5write(filename, '/rec/foldername', foldername_uint8);
    fprintf('  Saved foldername: %d entries\n', length(rec.foldername));

    % 3. fragmentname - cell array of strings
    fragmentname_chars = char(rec.fragmentname);
    fragmentname_uint8 = uint8(fragmentname_chars);
    h5create(filename, '/rec/fragmentname', size(fragmentname_uint8), 'Datatype', 'uint8');
    h5write(filename, '/rec/fragmentname', fragmentname_uint8);
    fprintf('  Saved fragmentname: %d entries\n', length(rec.fragmentname));

    % 4. boarders - numeric matrix
    h5create(filename, '/rec/boarders', size(rec.boarders));
    h5write(filename, '/rec/boarders', rec.boarders);
    fprintf('  Saved boarders: %dx%d matrix\n', size(rec.boarders));

    % 5. seqs - cell array of strings (DNA/protein sequences)
    seqs_chars = char(rec.seqs);
    seqs_uint8 = uint8(seqs_chars);
    h5create(filename, '/rec/seqs', size(seqs_uint8), 'Datatype', 'uint8');
    h5write(filename, '/rec/seqs', seqs_uint8);
    fprintf('  Saved seqs: %d sequences\n', length(rec.seqs));

    % 6. tag - cell array of strings
    tag_chars = char(rec.tag);
    tag_uint8 = uint8(tag_chars);
    h5create(filename, '/rec/tag', size(tag_uint8), 'Datatype', 'uint8');
    h5write(filename, '/rec/tag', tag_uint8);
    fprintf('  Saved tag: %d entries\n', length(rec.tag));

    % 7. recname - cell array of strings
    recname_chars = char(rec.recname);
    recname_uint8 = uint8(recname_chars);
    h5create(filename, '/rec/recname', size(recname_uint8), 'Datatype', 'uint8');
    h5write(filename, '/rec/recname', recname_uint8);
    fprintf('  Saved recname: %d entries\n', length(rec.recname));

    % 8. domloc - numeric matrix
    h5create(filename, '/rec/domloc', size(rec.domloc));
    h5write(filename, '/rec/domloc', rec.domloc);
    fprintf('  Saved domloc: %dx%d matrix\n', size(rec.domloc));

    % 9. domseq - cell array of strings
    domseq_chars = char(rec.domseq);
    domseq_uint8 = uint8(domseq_chars);
    h5create(filename, '/rec/domseq', size(domseq_uint8), 'Datatype', 'uint8');
    h5write(filename, '/rec/domseq', domseq_uint8);
    fprintf('  Saved domseq: %d sequences\n', length(rec.domseq));

    % 10. group - numeric vector
    h5create(filename, '/rec/group', size(rec.group));
    h5write(filename, '/rec/group', rec.group);
    fprintf('  Saved group: %d entries\n', length(rec.group));

    fprintf('Successfully saved all fields!\n');
end
