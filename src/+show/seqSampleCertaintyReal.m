% sample sequences, and draw for each colorful plots with what the posterior
% probability was compared to the correct state per letter
function seqSampleCertaintyReal(params, theta, dataset, outpathBase)
    fprintf('Showing real sequences and their epigenetics\n');
    trackNames = {'H3K27ac', 'DNase'};

    % tissues x 2
    bedGraphs = misc.readAllBedGraphs(dataset.tissueEIDs, trackNames);
    for tissueId = 1:length(dataset.tissueEIDs)
        seqInd = 1;
        for i = 1:10
            outpath = sprintf('%sEid%s.%d.jpg', outpathBase, dataset.tissueEIDs{tissueId}, i);
            tmp = plotSequence(params, dataset, theta, bedGraphs, seqInd, tissueId, false, outpath);
            seqInd = plotSequence(params, dataset, theta, bedGraphs, seqInd, tissueId, true, outpath);
        end
    end
    deeptoolsPlot(params, theta, dataset);
end


function seqInd = plotSequence(params, dataset, theta, bedGraphs, seqInd, tissueId, doSwitch, outpath)
    [N, L] = size(dataset.X);
    LOW_BAR_HIEGHT = 0.1;
    cMap = lines(params.m);
    [tracks, seqInd] = getSeqWithTracks(dataset, bedGraphs, seqInd + 1, tissueId);
    if seqInd == -1
        return
    end

    % gamma - 1 x m x L
    % psi - 1 x m x k x L
    [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, dataset.X(seqInd, :), ...
        dataset.pcPWMp(seqInd, :, :));
    % 1 x m x L
    posterior = calcPosterior(params, gamma, psi);
    YEst = misc.viterbi(params, theta, dataset.X(seqInd, :), dataset.pcPWMp(seqInd, :, :));
    if doSwitch
        tmp = posterior(1, 1, :);
        posterior(1, 1, :) = posterior(1, 2, :);
        posterior(1, 2, :) = tmp;
        tfs = YEst(:, :, 2);
        bgs = YEst(:, :, 1);
        bgs(bgs == 1) = -1;
        bgs(bgs == 2) = 1;
        bgs(bgs == -1) = 2;
        YEst = cat(3, bgs, tfs);
    end


    H3K27acTrack = tracks(:, :, 1);
    DNaseTrack = tracks(:, :, 2);

    pwmVal = params.m + 1;


    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    hold on;


    % LOW PLOT
    subplot(3, 1, 3);
    plotProbabilityMap(params, permute(posterior(1, :, :), [2, 3, 1]), YEst(1, :, :), cMap, ...
                       dataset.starts(seqInd), pwmVal, dataset.chrs{seqInd}, LOW_BAR_HIEGHT);

    % LEGEND
    ax = gca;
    ax.YDir = 'normal';
    legendStrings1 = strcat({'Enhancer Type '}, num2str([1:params.enhancerAmount]'));
    legendStrings2 = strcat({'Background '}, num2str([1:params.backgroundAmount]'));
    legendStrings = {legendStrings1{:}, legendStrings2{:}};
    legendStrings{pwmVal} = 'TFBS';
    legend(legendStrings);

    % HIGH PLOT
    subplot(3, 1, 1);
    chromatinPlot(cMap, H3K27acTrack, L);
    title(sprintf('Viterbi Classification of sequence %s:%d-%d (%s peak)', ...
        dataset.chrs{seqInd}, dataset.starts(seqInd), ...
        dataset.starts(seqInd) + L, dataset.tissueEIDs{dataset.Y(seqInd)}));
    % MIDDLE PLOT
    subplot(3, 1, 2);
    dnasePlot(params, dataset, seqInd, cMap, DNaseTrack, L, LOW_BAR_HIEGHT);

    saveas(gcf, outpath);
end


% tracks - tissues x L x tracks
function [tracks, seqInd] = getSeqWithTracks(dataset, bedGraphs, startInd, tissueId)
    N = size(dataset.X, 1);
    tissueSeqsInds = find(dataset.Y == tissueId);
    for seqInd = tissueSeqsInds(tissueSeqsInds >= startInd)'
        % fprintf('Trying sequence %d\n', seqInd)
        tracks = getTracks(dataset, bedGraphs, seqInd);
        if all(all(any(tracks > 0, 2), 1), 3)
            fprintf('Found seq %d\n', seqInd);
            return;
        end
    end
    fprintf('No sequence found with enough data in its track');
    seqInd = -1;
end


% tracks - tissues x L x tracks
function tracks = getTracks(dataset, bedGraphs, seqInd)
    L = size(dataset.X, 2);
    chr = dataset.chrs{seqInd};
    to = dataset.starts(seqInd) + L;
    from = dataset.starts(seqInd);
    tracks = zeros(size(bedGraphs, 1), L, size(bedGraphs, 2));
    for trackInd = 1:size(bedGraphs, 2)
        for i = 1:size(bedGraphs, 1)
            fprintf('looking for track %d of sequence %d %s[%d:%d], for tissue index %d\n', ...
                    trackInd, seqInd, chr, from, to, i);
            tracks(i, :, trackInd) = getTrack(bedGraphs{i, trackInd}, chr, from, to);
            if all(tracks(i, :, trackInd) == 0)
                return;
            end
        end
    end
end


function rotateYLabel()
    ylh = get(gca,'ylabel');
    gyl = get(ylh);
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation', 0, 'Position', ylp, 'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'right');
end


% posterior - N x m x L
function posterior = calcPosterior(params, gamma, psi)
    N = size(gamma, 1);
    posterior = gamma;
    for l = 1:params.k
        for t = 1:params.lengths(l)
            subStatePost = permute(cat(4, -inf(N, params.m, 1, t), psi(:, :, l, 1:end-t)), [1,2,4,3]);
            posterior = matUtils.logAdd(posterior, subStatePost);
        end
    end
    posterior = exp(posterior);
end


function ret = getTrack(bedGraph, trackChr, trackFrom, trackTo)
    MIN_SAMPLES_IN_SEQ = 100;
    mask = strcmp(bedGraph.chrs, trackChr) & (bedGraph.tos >= trackFrom) & (bedGraph.froms <= trackTo);
    % dilation of bin vector
    mask = [mask(2:end); false] | mask | [false; mask(1:end-1)];
    % linear interpa
    foundSamples = sum(mask);
    fprintf('Found %d sample points in (%s:%d-%d)\n', foundSamples, trackChr, trackFrom, trackTo);
    if foundSamples < MIN_SAMPLES_IN_SEQ
        ret = zeros(trackTo - trackFrom, 1);
        fprintf('%d < %d, not enough samples in sequence', foundSamples, MIN_SAMPLES_IN_SEQ);
        return
    end
    knownPoints = double([(bedGraph.tos(mask) + bedGraph.froms(mask)) / 2]');
    knownVals = bedGraph.vals(mask)';
    [knownPoints, inds] = unique(knownPoints);
    knownVals = knownVals(inds);
    wantedPoints = double([trackFrom : trackTo - 1]);
    ret = interp1(knownPoints, knownVals, wantedPoints);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% YEst - 1 x L x 2
% probMap - m x L
function plotProbabilityMap(params, probMap, YEst, cMap, start, pwmVal, chr, lowBarHieght)
    BARS_PLOT_DARKNESS_FACTOR = 0.85;
    PWM_COLOR = [0, 0, 0];

    cMapWithError = [cMap;  PWM_COLOR];
    % probMap = permute(probMap, [2, 3, 1]);
    L = size(YEst, 2);
    hold on;
    % if params.m >= size(probMap, 1) >= params.enhancerAmount
    % m x L
    YEstOneHot = matUtils.vec2mat(YEst(:, :, 1), params.m);
    selectedProb = probMap .* YEstOneHot(1:size(probMap, 1), :);
    % m + 1 x L
    barHandle = bar(start: start + L - 1, selectedProb', 1, 'stacked', 'FaceColor', 'flat');
    for j = 1:size(selectedProb, 1)
        barHandle(j).CData = cMap(j, :) * BARS_PLOT_DARKNESS_FACTOR;
    end
    % end

    plotHandle = plot(start: start + L - 1, probMap', 'LineWidth', 1.5);
    cellCMap = {};
    for i = 1:size(probMap, 1)
        cellCMap{i, 1} = cMap(i, :);
    end
    set(plotHandle, {'Color'}, cellCMap);
    rotateYLabel();


    % Comma separated ticks
    xlim([start, start + L - 1]);
    ax = gca;
    ax.XAxis.Exponent = 0;
    xtickformat('%,d')

    YEstColored = YEst;
    YEstColored(:, YEst(:, :, 2) > 0, 1) = pwmVal;
    imagesc([start - .5, start + L - 1.5], [-3 * lowBarHieght / 4, -lowBarHieght / 4], ...
            [YEstColored(:, :, 1); YEstColored(:, :, 1)], [1, pwmVal]);
    ylim([-lowBarHieght, 1]);
    colormap(gca, cMapWithError);
    text(L + 1, 0.5, 'Posterior Prob.', 'FontSize', 10);
    text(L + 1, -lowBarHieght / 2, 'Viterbi Est.', 'FontSize', 10);
    xlabel(sprintf('Position in %s', chr));
    ylabel('P(y_{t}|X)');
    hold off;
end

function chromatinPlot(cMap, H3K27acTrack, L)
    hold on;
    set(gca,'xtick',[]);
    text(L + 1, 0.5, 'H3K27ac', 'FontSize', 10);
    plotHandle = plot(1:L, H3K27acTrack', 'LineWidth', 1.5);
    cellCMap = {};
    for i = 1:size(H3K27acTrack, 1)
        cellCMap{i, 1} = cMap(i, :);
    end
    set(plotHandle, {'Color'}, cellCMap);
    xlim([1, L]);
    ylim([0, max(H3K27acTrack(:))]);
    ylabel('-log_{10}(p-value)');
    hold off;
end


function dnasePlot(params, dataset, seqInd, cMap, DNaseTrack, L, lowBarHieght)
    hold on;
    lowBarHieghtRel = max(DNaseTrack(:)) * lowBarHieght;
    set(gca,'xtick',[]);
    text(L + 1, 0.5, 'DNase', 'FontSize', 10);
    plotHandle = plot(1:L, DNaseTrack', 'LineWidth', 1.5);
    cellCMap = {};
    for i = 1:size(DNaseTrack, 1)
        cellCMap{i, 1} = cMap(i, :);
    end
    set(plotHandle, {'Color'}, cellCMap);
    xlim([1, L]);
    ylim([-lowBarHieghtRel, max(DNaseTrack(:))]);
    ylabel('-log_{10}(p-value)');
    plotChromHMM(dataset.tissueEIDs(1: 2), dataset.chrs{seqInd}, ...
        dataset.starts(seqInd), L, lowBarHieghtRel);
    hold off;
end


function plotChromHMM(eids, chr, from, L, lowBarHieghtRel)
    CHROMHMM_COLORMAP = [255, 0, 0;
                         255, 69, 0;
                         50, 205, 50;
                         0, 128, 0;
                         0, 100, 0;
                         194, 225, 5;
                         255, 255, 0;
                         102, 205, 170;
                         138, 145, 208;
                         205, 92, 92;
                         233, 150, 122;
                         189, 183, 107;
                         128, 128, 128;
                         192, 192, 192;
                         255, 255, 255];
    MARKS_NAMES = {'Active TSS',
                   'Flanking Active TSS Orange',
                   "Transcr. at gene 5' and 3'",
                   'Strong transcription',
                   'Weak transcription',
                   'Genic enhancers',
                   'Enhancers',
                   'ZNF genes & repeats Medium',
                   'Heterochromatin',
                   'Bivalent/Poised TSS',
                   'Flanking Bivalent TSS/Enh',
                   'Bivalent Enhancer',
                   'Repressed PolyComb',
                   'Weak Repressed PolyComb',
                   'Quiescent/Low'};

    to = from + L;

    chromHMMmarks = zeros(length(eids), L);
    for i = 1:length(eids)
        bedFilePath = sprintf('../data/ChromHMM/%s_15_coreMarks_mnemonics.bed', eids{i});
        fid = fopen(bedFilePath);
        bedData = textscan(fid, '%s%d%d%d_%s', 'delimiter','\t');
        chrs = bedData{1};
        froms = bedData{2};
        tos = bedData{3};
        marks = bedData{4};
        markStrs = bedData{5};
        fclose(fid);

        mask = strcmp(chrs, chr) & (froms <= to) & (tos >= from);
        subFroms = froms(mask);
        subFroms(1) = from;
        subTos = tos(mask);
        subTos(end) = to;
        subMarks = marks(mask);
        diffs = subTos - subFroms;
        chromHMMmarks(i, :) = repelem(subMarks, diffs);
    end

    imagesc([.5, L - .5], [-3 * lowBarHieghtRel / 4, -lowBarHieghtRel / 4], chromHMMmarks - 1, [0, 15]);
    colormap(gca, CHROMHMM_COLORMAP ./ 255);
    c = colorbar;
    set(c, 'Position', get(c, 'Position') + [0.05, 0, 0, 0]);
    c.Ticks = .5: 1: 15.5;
    c.TickLabels = MARKS_NAMES;
end


% Deeptools-like plot
function deeptoolsPlot(params, theta, dataset)
    [N, L] = size(dataset.X);
    subsampleMask = 1: N < 500;
    [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, dataset.X(subsampleMask, :), ...
        dataset.pcPWMp(subsampleMask, :, :));
    % N x m x L
    posterior = calcPosterior(params, gamma, psi);
    isEnhancer = max(posterior(:, 1:2, :), [], 2);
    isEnhancer = permute(isEnhancer, [1,3,2]);
    [~, seqsInds] = sort(isEnhancer(:, floor(L / 2)), 'desc');
    isEnhancer = isEnhancer(seqsInds, :);
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    subplot(1, 2, 1);
    imagesc(isEnhancer);
    colormap(gca, 'Gray');
    subplot(1, 2, 2);
    imagesc(dataset.Y(seqsInds));
end

