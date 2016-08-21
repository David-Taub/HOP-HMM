function k = reader()
    suffix = '-H3K27ac.peaks.mat';
    homePath = '/cs/cbio/tommy/Enhancers/Data/GSE29184_Bing_Ren';
    fileList = dir(homePath);
    filtered = regexp({fileList.name} ,'(.+.H3K27ac.peaks.mat$)','match');
    filtered = [filtered{:}];
    cellLines = strrep(filtered, suffix, '')
    N = length(cellLines);
    k = zeros(N);
    for i = 1 : N
        i
        f = fullfile(homePath, strcat(cellLines{i}, suffix))
        load(f, 'S');
        S2=[S{:}];
        heightTreshold = 0.1;
        % singular = sum(cat(1,S2.overlap), 2) == 1;
        distal = ismember({S2.class},{'Distal'}).';
        high = ([S2.height] > S{ceil(length(S) * heightTreshold)}.height).';
        % mask = singular & distal & high;
        mask = distal & high;
        S3 = [S{mask}];
        k(i, :) = sum(cat(1, S3.overlap), 1);
        imagesc(k)
        drawnow
    end
    k;
    imagesc(k)
    title('Normalized')
    figure
    imagesc(bsxfun(@times, k, 1 ./ max(k, [], 2)));
    title('Real')
end