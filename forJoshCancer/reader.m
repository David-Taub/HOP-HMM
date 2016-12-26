function dataset = reader()
    filenames = {...
                 {'healthy_train', false}, ...
                 {'healthy_test', false}, ...
                 {'cancer_train', true}, ...
                 {'cancer_test', true}, ...
                };
    dirpath = '/cs/cbio/david/projects/CompGenetics/forJoshCancer/';
    dataset.tissueType = {};
    dataset.samples = {};
    dataset.peopleNum = [];
    dataset.isCancer = [];
    isCancer = [];
    tissueType = {};
    samples = [];
    for i = 1:length(filenames)
        tic
        filename = filenames{i}{1}
        isCancerFile = filenames{i}{2};

        ssCsvPath = strcat(dirpath, 'ss_', filename, '.csv');
        CSVPath = strcat(dirpath, filename, '.csv');

        % original data
        T = readtable(CSVPath);
        data = table2array(T(2:end,2:end));
        % sitesNames = table2array(T(2:end, 1));
        % samplesNames = table2array(T(1, 2:end));

        % tissue type
        ss = readtable(ssCsvPath);
        tissueType = cat(1, tissueType, ss(:,10));
        samples = cat(2, samples, data);
        isCancer(end+1:end+size(ss,1)) = isCancerFile;
        toc
    end
    keyboard
    [dataset.tissueType, ~, indices] = unique(tissueType);
    l = length(dataset.tissueType);
    
    for c = [true, false]
        for i = 1:l
            dataset.peopleNum(i+ c * l) = sum(indices == i && isCancer == c);
            dataset.isCancer(i+ c * l) = c;
            dataset.samples{i+ c * l} = samples(:, indices == i && isCancer == c);
        end
    end
end
