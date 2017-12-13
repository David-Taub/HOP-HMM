
function josh(origM)
    close all;
    M = origM;
    smalls = 0 : 0.01 : 0.3;
    bigs = 0 : 0.01 : 0.5;
    covs = zeros(length(smalls), length(bigs));
    sits = zeros(length(smalls), length(bigs));
    
    for si = 1 : length(smalls)
        parfor bi = 1 : length(bigs)
            [si, bi]
            if smalls(si) > bigs(bi)
                continue
            end
            [covered, sites] = getSites(origM, smalls(si), bigs(bi));
            covs(si, bi) = length(covered);
            sits(si, bi) = length(sites);
        end
        subplot(1,2,1); imagesc(bigs, smalls, covs);
        title('Covered');
        xlabel('Big Gap TH');
        ylabel('Small Gap TH');
        colorbar;
        subplot(1,2,2); imagesc(bigs, smalls, sits);
        title('Sites');
        xlabel('Big Gap TH');
        ylabel('Small Gap TH');
        colorbar;
        drawnow;
    end
end
% sites are the indices of CpG sites we use to cover the tissues
function [covered, sites] = getSites(origM, smallGapTh, bigGapTh)
    
    pos1Peaks = [];
    neg1Peaks = [];
    covered = [];
    while true;
        % iterations over the matrix
        C = length(covered);
        % covered is the indices of the tissues we already "got"
        [covered, pos1PeakCur, neg1PeakCur] = ...
            get1Peaks(origM, covered, smallGapTh, bigGapTh);
        pos1Peaks = [pos1Peaks; pos1PeakCur];
        neg1Peaks = [neg1Peaks; neg1PeakCur];
        if C == length(covered) || length(covered) >= size(origM, 2) - 1
            % no new tissues were added
            break
        end
    end
    sites = [pos1Peaks; neg1Peaks];
end

function [covered, pos1Peaks, neg1Peaks] = get1Peaks(origM, alreadyCovered, smallGapTh, bigGapTh)
    M = origM;
    % removed covered columns
    M(:, alreadyCovered) = [];

    [N,m] = size(M);
    
    % N x m
    % find the interesting "peaks" lines
    MSorted = sort(M, 2, 'descend');

    % boolean vector, 1 means the line is a peak, 0 otherwise. 
    % N x 1
    pos1Peaks = find(MSorted(:, 1) - MSorted(:, 2) > bigGapTh & ...
                    MSorted(:, 2) - MSorted(:, m) < smallGapTh);
    neg1Peaks = find(MSorted(:, m-1) - MSorted(:, m) > bigGapTh & ...
                    MSorted(:, 1) - MSorted(:, m-1) < smallGapTh);
    
    
    % n x 1, where n is the number of peaks
    % get the tissue types of the peaks (negative and positive)
    [~, pos1PeaksType] = max(origM(pos1Peaks, :), [], 2);
    [~, neg1PeaksType] = min(origM(neg1Peaks, :), [], 2);
    covered = unique([alreadyCovered; pos1PeaksType; neg1PeaksType]);
end

   

function filtered = filterM(M, vals, a)
    N = length(vals);
    s = sort(vals, 'descend');
    threshold = s(ceil(a * N));
    filtered = M(vals > threshold, :);
end

function readerCSV()
    data = importdata('betas_means.csv'); 
    origM = data.data;
    
    cellTypes = data.textdata(1, 2:end);

    sites = data.textdata(2:end, 1);

    clear data
    save('M.mat');
end

