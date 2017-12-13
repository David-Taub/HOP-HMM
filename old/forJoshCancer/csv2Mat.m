% cd /cs/cbio/david/projects/CompGenetics/forJoshCancer
% dirname = '/cs/cbio/david/projects/CompGenetics/forJoshCancer/';
% csv2Mat(dirname,'cancer_test'); csv2Mat(dirname,'cancer_train'); csv2Mat(dirname,'health_test'); csv2Mat(dirname,'health_train');

function csv2Mat(dirpath, filename)
    fullCSVPath = strcat(dirpath, filename, '.csv')
    fullMatPath = strcat(dirpath, filename, '.mat')

    T = readtable(fullCSVPath);
    data = table2array(T(2:end,2:end));
    sitesNames = table2array(T(2:end, 1));
    samplesNames = table2array(T(1, 2:end));

    save(fullMatPath, 'data', 'sitesNames', 'samplesNames', '-v7.3');
end

