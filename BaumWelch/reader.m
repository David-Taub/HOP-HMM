function [peaks] = reader()
    suffix = '-H3K27ac.peaks.mat';
    homePath = '/cs/cbio/tommy/Enhancers/Data/GSE29184_Bing_Ren';
    fileList = dir(homePath);
    filtered = regexp({fileList.name} ,'(.+.H3K27ac.peaks.mat$)','match');
    filtered = [filtered{:}];
    cellLines = strrep(filtered, suffix, '');
    N = length(cellLines);
    chrs = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrM', 'chrX', 'chrY'};

    peaks = [];
    for i = 1 : N
        i
        f = fullfile(homePath, strcat(cellLines{i}, suffix));
        load(f, 'S');
        S2=[S{:}];
        heightTreshold = 0.5;
        % singular = sum(cat(1,S2.overlap), 2) == 1;
        distal = ismember({S2.class},{'Distal'}).';
        high = ([S2.height] > S{ceil(length(S) * heightTreshold)}.height).';
        % mask = singular & distal & high;
        mask = distal & high;
        S3 = [S{mask}];
        M = length(S3);
        [~,chrInd] = ismember({S3.chr}, chrs);
        peaksAdd(1:M, 1:i-1) = 0;
        peaksAdd(1:M, i) = 1;
        peaksAdd(1:M, i+1:N) = 0;
        peaksAdd(1:M, N+1) = [S3.from] + (10^9 .* chrInd);
        peaksAdd(1:M, N+2) = [S3.to] + (10^9 .* chrInd);
        peaks = cat(1, peaks, peaksAdd);
    end
end
