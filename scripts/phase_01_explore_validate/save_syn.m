% Usage:
load('/home/wulong/RotationWorks/PhylogenicTree/data/raw/Burkholderiaceae/syn.mat');
saveSynToHDF5('/home/wulong/RotationWorks/PhylogenicTree/data/interim/syn.h5', syn);

function saveSynToHDF5(filename, syn)
    % Delete existing file
    if exist(filename, 'file')
        delete(filename);
    end

    fprintf('Saving syn struct to %s...\n', filename);

    % Helper function to safely convert cell array to char matrix
    function char_mat = safeCellToChar(cellArray)
        % Replace empty cells with space
        cellArray = cellfun(@(x) char(x), cellArray, 'UniformOutput', false);
        for i = 1:length(cellArray)
            if isempty(cellArray{i})
                cellArray{i} = ' ';
            end
        end
        char_mat = char(cellArray);
    end

    % 1. clusterblast - cell array of strings
    clusterblast_chars = safeCellToChar(syn.clusterblast);
    clusterblast_uint8 = uint8(clusterblast_chars);
    h5create(filename, '/syn/clusterblast', size(clusterblast_uint8), 'Datatype', 'uint8');
    h5write(filename, '/syn/clusterblast', clusterblast_uint8);
    fprintf('  Saved clusterblast: %d entries\n', length(syn.clusterblast));

    % 2. strainName - cell array of strings
    strainName_chars = safeCellToChar(syn.strainName);
    strainName_uint8 = uint8(strainName_chars);
    h5create(filename, '/syn/strainName', size(strainName_uint8), 'Datatype', 'uint8');
    h5write(filename, '/syn/strainName', strainName_uint8);
    fprintf('  Saved strainName: %d entries\n', length(syn.strainName));

    % 3. regionName - cell array of strings
    regionName_chars = safeCellToChar(syn.regionName);
    regionName_uint8 = uint8(regionName_chars);
    h5create(filename, '/syn/regionName', size(regionName_uint8), 'Datatype', 'uint8');
    h5write(filename, '/syn/regionName', regionName_uint8);
    fprintf('  Saved regionName: %d entries\n', length(syn.regionName));

    % 4. assemblyDefinition - cell array of strings
    assemblyDefinition_chars = safeCellToChar(syn.assemblyDefinition);
    assemblyDefinition_uint8 = uint8(assemblyDefinition_chars);
    h5create(filename, '/syn/assemblyDefinition', size(assemblyDefinition_uint8), 'Datatype', 'uint8');
    h5write(filename, '/syn/assemblyDefinition', assemblyDefinition_uint8);
    fprintf('  Saved assemblyDefinition: %d entries\n', length(syn.assemblyDefinition));

    % 5. location - cell array of numeric arrays
    % Convert to matrix if all are 1x2
    try
        location_mat = cell2mat(syn.location);
        h5create(filename, '/syn/location', size(location_mat));
        h5write(filename, '/syn/location', location_mat);
        fprintf('  Saved location: %dx%d matrix\n', size(location_mat));
    catch
        % If cell2mat fails, save as JSON
        fprintf('  Location has irregular structure, saving as JSON...\n');
        location_json = jsonencode(syn.location);
        location_uint8 = uint8(location_json);
        h5create(filename, '/syn/location_json', [1, length(location_uint8)], 'Datatype', 'uint8');
        h5write(filename, '/syn/location_json', location_uint8);
        fprintf('  Saved location: %d entries (as JSON)\n', length(syn.location));
    end

    % 6. wholeCDS - nested cell array (save as JSON)
    fprintf('  Saving wholeCDS as JSON (this may take a while)...\n');
    wholeCDS_json = jsonencode(syn.wholeCDS);
    wholeCDS_uint8 = uint8(wholeCDS_json);
    h5create(filename, '/syn/wholeCDS_json', [1, length(wholeCDS_uint8)], 'Datatype', 'uint8');
    h5write(filename, '/syn/wholeCDS_json', wholeCDS_uint8);
    fprintf('  Saved wholeCDS: %d entries (as JSON)\n', length(syn.wholeCDS));

    % 7. biosynCDS - nested cell array (save as JSON)
    fprintf('  Saving biosynCDS as JSON (this may take a while)...\n');
    biosynCDS_json = jsonencode(syn.biosynCDS);
    biosynCDS_uint8 = uint8(biosynCDS_json);
    h5create(filename, '/syn/biosynCDS_json', [1, length(biosynCDS_uint8)], 'Datatype', 'uint8');
    h5write(filename, '/syn/biosynCDS_json', biosynCDS_uint8);
    fprintf('  Saved biosynCDS: %d entries (as JSON)\n', length(syn.biosynCDS));

    % 8. group - numeric vector
    h5create(filename, '/syn/group', size(syn.group));
    h5write(filename, '/syn/group', syn.group);
    fprintf('  Saved group: %d entries\n', length(syn.group));

    % 9. regionIdentifier - cell array of strings
    regionIdentifier_chars = safeCellToChar(syn.regionIdentifier);
    regionIdentifier_uint8 = uint8(regionIdentifier_chars);
    h5create(filename, '/syn/regionIdentifier', size(regionIdentifier_uint8), 'Datatype', 'uint8');
    h5write(filename, '/syn/regionIdentifier', regionIdentifier_uint8);
    fprintf('  Saved regionIdentifier: %d entries\n', length(syn.regionIdentifier));

    fprintf('Successfully saved all fields!\n');
end
