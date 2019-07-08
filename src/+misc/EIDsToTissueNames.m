function tissueNames = EIDsToTissueNames(tissueEIDs)
    ROADMAP_NAMES_CSV_PATH = '../data/peaks/help/full_tissue_names.csv';
    assert(isfile(ROADMAP_NAMES_CSV_PATH));
    fid = fopen(ROADMAP_NAMES_CSV_PATH);
    csvData = textscan(fid, '%s%s', 'delimiter',',');
    fclose(fid);
    namesDict = containers.Map(csvData{1}, csvData{2});
    tissueNames = tissueEIDs;
    for i = 1:length(tissueNames)
        if any(strcmp(tissueEIDs{i}, namesDict.keys))
            EID = tissueEIDs{i};
            tissueNames{i} = namesDict(EID);
            fprintf('Tissue name found: %s -> %s\n', EID, tissueNames{i});
        end
    end
end

