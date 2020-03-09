
function selectedPWMs = PWMsFeatureSelect(mergedPeaksMin, wantedK)
    fprintf('Selecting PWMs\n');


    dbstop if error
    % PWM - k x n x J
    % lengths - k x 1
    % names - k x 1
    [PWMs, lengths, names] = misc.PWMs();
    m = size(mergedPeaksMin.overlaps, 2);
    L = size(mergedPeaksMin.seqs, 2);
    [k, n, J] = size(PWMs);
    outFilepath = sprintf('../data/precomputation/SelectedPWMs_L%dwk%dm%dtissues%s.mat', L,...
                          wantedK, m, sprintf('%d', mergedPeaksMin.tissueList));

    if isfile(outFilepath)
        fprintf('Found selected PWMS cache file at %s\n', outFilepath);
        load(outFilepath);
        return;
    end
    [X, Y] = preprocess(mergedPeaksMin);
    Xs1H = matUtils.mat23Dmat(X, n);


    selectedPWMs = [];
    selectedPWMs = getExpressedPWMs(mergedPeaksMin, PWMs, names, selectedPWMs);
    fprintf('Selected PWMs by strength: %s\n', sprintf('%d ', selectedPWMs));

    selectedPWMs = getStrongPWMs(Xs1H, PWMs, lengths, wantedK, selectedPWMs);
    fprintf('Selected PWMs by strength: %s\n', sprintf('%d ', selectedPWMs));

    selectedPWMs = getDistinguishablePWMs(Xs1H, Y, PWMs, lengths, wantedK, selectedPWMs);
    fprintf('Added PWMs by AUC-ROC between sequences: %s\n', sprintf('%d ', selectedPWMs));

    save(outFilepath, 'selectedPWMs');
    fprintf('Saved feature selected PWMs in %s\n', outFilepath);
end

function selectedPWMs = getExpressedPWMs(mergedPeaksMin, PWMs, names, selectedPWMs)
    THREHSOLD = 0.4;
    MAX_REVERSE_MAP = 3;
    mat = readtable('..\data\Jaspar\raw\TF_expressions.txt');
    TFNames = mat.gene_id;
    PWMsNames = upper([names{:}])';
    reverseIds = {};
    for i = 1:length(TFNames)
        ids = find(cellfun(@isempty, strfind(PWMsNames, TFNames{i})) == 0);
        if ~isempty(ids)
            if length(ids) > MAX_REVERSE_MAP
                ids = ids(1:MAX_REVERSE_MAP);
            end
            reverseIds{i} = ids';
        end
    end
    for EID = mergedPeaksMin.tissueEIDs
        EidMask = strcmp(EID, mat.Properties.VariableNames);
        if any(EidMask)
            tfPresence = mat(:, EidMask);
            [tfPresenceSorted, tfPresenceIds]  = sort(tfPresence{:, :}, 'desc');
            tfPresenceIdsBest = tfPresenceIds(tfPresenceSorted > THREHSOLD);
            expressedPWMsIds = [reverseIds{tfPresenceIdsBest}];
            selectedPWMs = [selectedPWMs, expressedPWMsIds];
        end
    end
    selectedPWMs
end

function selectedPWMs = getDistinguishablePWMs(Xs1H, Y, PWMs, lengths, wantedK, selectedPWMs)
    aucRocs = misc.oneVsAllAucRoc(Xs1H, Y, PWMs, lengths);
    [aucRocsSorted, aucRocsSortedInd] = sort(aucRocs, 2, 'descend');

    i = 1;
    while length(selectedPWMs) < wantedK
        selectedPWMs = [selectedPWMs, unique(aucRocsSortedInd(1:i))];
        i = i + 1;
    end
end

function selectedPWMs = getStrongPWMs(Xs1H, PWMs, lengths, wantedK, selectedPWMs)
    EXPECTED_NUM_OF_PEAKS_IN_SEQ = 3;
    PICK_STRONGEST_RATIO = 0.3;
    [k, n, ~] = size(PWMs);
    rates = zeros(1, k);
    for pwmId = 1:k
        PWMLogLike = misc.PWMLogLikelihood(PWMs, lengths, Xs1H, pwmId);
        % mask x EXPECTED_NUM_OF_PEAKS_IN_SEQ
        rates(pwmId) = mean(mean(misc.maxN(PWMLogLike, 2, EXPECTED_NUM_OF_PEAKS_IN_SEQ), 2), 1);
        rates(pwmId) = rates(pwmId) - log(1 / n) * lengths(pwmId);
        fprintf('Checked PWM %d - log likelihood: %.2f (%.2f), length: %d\n', pwmId, rates(pwmId), max(rates), lengths(pwmId));
    end
    [rates, inds] = sort(rates, 'descend');
    selectedPWMs = [selectedPWMs, inds(1: floor(wantedK * PICK_STRONGEST_RATIO))];
end

function [X, Y] = preprocess(mergedPeaksMin)
    X = mergedPeaksMin.seqs;
    Y = (mergedPeaksMin.overlaps > 0) * (1:size(mergedPeaksMin.overlaps, 2))';
    [overlaps, seqInd] = sortrows(mergedPeaksMin.overlaps);
    X = X(seqInd, :);
    Y = Y(seqInd);
end
