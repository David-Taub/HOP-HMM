% sample sequences, and draw for each colorful plots with what the posterior
% probability was compared to the correct state per letter
function seqSampleCertainty5Real(params, theta, dataset, outpathBase)
    fprintf('Showing real sequences and their epigenetics\n');
    trackNames = {'H3K27ac', 'DNase'};

    % tissues x 2
    bedGraphs = misc.readAllBedGraphs(dataset.tissueEIDs, trackNames);
    for tissueId = 1:length(dataset.tissueEIDs)
        seqInd = 1;
        for i = 1:40
            close all;
            outpath = sprintf('%s.%svs%s.%dS0.jpg', outpathBase, dataset.tissueEIDs{tissueId}, dataset.tissueEIDs{3-tissueId}, i);
            tmp = plotSequence(params, dataset, theta, bedGraphs, seqInd, tissueId, false, outpath);
            outpath = sprintf('%s.%svs%s.%dS1.jpg', outpathBase, dataset.tissueEIDs{tissueId}, dataset.tissueEIDs{3-tissueId}, i);
            seqInd = plotSequence(params, dataset, theta, bedGraphs, seqInd, tissueId, true, outpath);
            if seqInd == -1
                return
            end
        end
    end
end


function seqInd = plotSequence(params, dataset, theta, bedGraphs, seqInd, tissueId, doSwitch, outpath)
    [N, L] = size(dataset.X);
    LOW_BAR_HEIGHT = 0.1;
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
    % 1 x L x 2
    YEstViterbi = misc.viterbi(params, theta, dataset.X(seqInd, :), dataset.pcPWMp(seqInd, :, :));
    % 1 x 1 x L
    [~, YEstMax] = max(posterior, [], 2);

    % 1 x m x k+1 x L
    a = cat(3, permute(gamma, [1, 2, 4, 3]), psi);
    % 1 x m x 1 x L
    [v, indsTfs] = max(a, [], 3);
    % 1 x 1 x 1 x L
    [~, inds] = max(v, [], 2);
    % 1 x L x 2
    YEstMax = cat(3, inds(:)', indsTfs(sub2ind(size(indsTfs), ones(1,L), inds(:)', ones(1,L), 1:L)) - 1);

    % 1 x m x 1 x L
    if doSwitch
        tmp = posterior(1, 1, :);
        posterior(1, 1, :) = posterior(1, 2, :);
        posterior(1, 2, :) = tmp;
        tfs = YEstViterbi(:, :, 2);
        bgs = YEstViterbi(:, :, 1);
        bgs(bgs == 1) = -1;
        bgs(bgs == 2) = 1;
        bgs(bgs == -1) = 2;
        YEstViterbi = cat(3, bgs, tfs);

        tfs = YEstMax(:, :, 2);
        bgs = YEstMax(:, :, 1);
        bgs(bgs == 1) = -1;
        bgs(bgs == 2) = 1;
        bgs(bgs == -1) = 2;
        YEstMax = cat(3, bgs, tfs);
    end



    H3K27acTrack = tracks(:, :, 1);
    DNaseTrack = tracks(:, :, 2);

    pwmVal = params.m + 1;
    xRange = [dataset.starts(seqInd): dataset.starts(seqInd) + L - 1];

    figure('units', 'pixels', 'Position', [0 0 1500 1000]);
    % 1
    subplot(5, 1, 1);

    areaHandle = area(xRange, H3K27acTrack(1, :), 'LineWidth', 1.5);
    set(gca,'XTickLabel',[]);
    areaHandle(1).FaceColor = cMap(1, :);
    xlim([xRange(1), xRange(end)]);
    ylim([0, max(H3K27acTrack(:))]);


    % TF names
    tfs = YEstViterbi(1, :, 2);
    tfs = diff(tfs);
    tfs = tfs(tfs > 0);
    tfsNames = params.names(tfs);
    tfsNames = [tfsNames{:}];
    if length(tfsNames) > 0
        stfs = sprintf('%s, ', tfsNames{:});
    else
        stfs = '';
    end

    % Title
    title(sprintf('%s vs %s [%s:%d-%d] [%s]', ...
        dataset.tissueEIDs{1 + doSwitch}, dataset.tissueEIDs{2 - doSwitch}, ...
        dataset.chrs{seqInd}, dataset.starts(seqInd), dataset.starts(seqInd) + L, stfs));

    % 2
    subplot(5, 1, 2);
    hold on;
    areaHandle = area(xRange, DNaseTrack(1, :), 'LineWidth', 1.5);
    set(gca,'XTickLabel',[]);
    areaHandle(1).FaceColor = cMap(1, :);
    xlim([xRange(1), xRange(end)]);
    ylim([-LOW_BAR_HEIGHT * max(DNaseTrack(:)), max(DNaseTrack(:))]);
    plotChromHMM(dataset.tissueEIDs{1}, dataset.chrs{seqInd}, ...
        dataset.starts(seqInd), L, LOW_BAR_HEIGHT * max(DNaseTrack(:)));
    hold off;


    % 3
    subplot(5, 1, 3);
    set(gca,'XTickLabel',[]);
    plotProbabilityMap(params, permute(posterior(1, :, :), [2, 3, 1]), YEstViterbi(1, :, :), YEstMax(1, :, :), cMap, ...
                       dataset.starts(seqInd), pwmVal, dataset.chrs{seqInd}, LOW_BAR_HEIGHT);

    % 4
    subplot(5, 1, 4);
    areaHandle = area(xRange, H3K27acTrack(2, :), 'LineWidth', 1.5);
    set(gca,'XTickLabel',[]);
    areaHandle(1).FaceColor = cMap(2, :);
    xlim([xRange(1), xRange(end)]);
    ylim([0, max(H3K27acTrack(:))]);

    % 5
    subplot(5, 1, 5);
    hold on;
    areaHandle = area(xRange, DNaseTrack(2, :), 'LineWidth', 1.5);
    areaHandle(1).FaceColor = cMap(2, :);
    xlim([xRange(1), xRange(end)]);
    ylim([-LOW_BAR_HEIGHT * max(DNaseTrack(:)), max(DNaseTrack(:))]);
    plotChromHMM(dataset.tissueEIDs{2}, dataset.chrs{seqInd}, ...
        dataset.starts(seqInd), L, LOW_BAR_HEIGHT * max(DNaseTrack(:)));
    hold off;

    % X ticks
    ax = gca;
    ax.XAxis.Exponent = 0;
    xtickformat('%,d')

    % % LEGEND
    % ax = gca;
    % ax.YDir = 'normal';
    % legendStrings1 = strcat({'Enhancer Type '}, num2str([1:params.enhancerAmount]'));
    % legendStrings2 = strcat({'Background '}, num2str([1:params.backgroundAmount]'));
    % legendStrings = {legendStrings1{:}, legendStrings2{:}};
    % legendStrings{pwmVal} = 'TFBS';
    % legend(legendStrings);

    % Save
    saveas(gcf, outpath);
end


% tracks - tissues x L x tracks
function [tracks, seqInd] = getSeqWithTracks(dataset, bedGraphs, startInd, tissueId)
    [N, L] = size(dataset.X);
    tissueSeqsInds = find(dataset.Y == tissueId);
    for seqInd = tissueSeqsInds(tissueSeqsInds >= startInd)'
        % fprintf('Trying sequence %d\n', seqInd)
        tracks = getTracks(dataset, bedGraphs, seqInd);
        if all(all(any(tracks > 0, 2), 1), 3)
            fprintf('Sequence %d - %s[%s:%s]\n', seqInd, dataset.chrs{seqInd}, ...
                dataset.starts(seqInd), dataset.starts(seqInd) + L);
            return;
        end
    end
    fprintf('No sequence found with enough data in its track\n');
    seqInd = -1;
    tracks = [];
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
            % fprintf('looking for track %d of sequence %d %s[%d:%d], for tissue index %d\n', ...
            %         trackInd, seqInd, chr, from, to, i);
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
    MAX_SAMPLE_RATE = 1 / 5;

    mask = strcmp(bedGraph.chrs, trackChr) & (bedGraph.tos >= trackFrom) & (bedGraph.froms <= trackTo);
    % dilation of bin vector
    mask = [mask(2:end); false] | mask | [false; mask(1:end-1)];
    % linear interpa
    % fprintf('Found %d sample points in (%s:%d-%d)\n', foundSamples, trackChr, trackFrom, trackTo);
    knownPoints = sort([bedGraph.froms(mask)' + 1, bedGraph.tos(mask)']);
    knownPoints(1) = trackFrom;
    knownPoints(end) = trackTo;

    if length(knownPoints) * MAX_SAMPLE_RATE < 1 | ...
            max(diff(knownPoints)) > (trackTo - trackFrom) * MAX_SAMPLE_RATE
        ret = zeros(trackTo - trackFrom, 1);
        fprintf('samples: %d - max distance: %d\n', length(knownPoints), max(diff(knownPoints)))
        return
    end
    knownVals = repelem(bedGraph.vals(mask)', 2);
    [knownPoints, inds] = unique(knownPoints);
    knownVals = knownVals(inds);
    wantedPoints = double([trackFrom : trackTo - 1]);
    ret = interp1(double(knownPoints), double(knownVals), wantedPoints);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% YEst - 1 x L x 2
% probMap - m x L
function plotProbabilityMap(params, probMap, YEstViterbi, YEstMax, cMap, start, pwmVal, chr, lowBarHeight)
    BARS_PLOT_DARKNESS_FACTOR = 0.85;
    PWM_COLOR = [0, 0, 0];
    cMapWithError = [cMap;  PWM_COLOR];
    L = size(YEstViterbi, 2);

    hold on;
    % Fill by Viterbi (darker color)

    % m x L
    YEstOneHot = matUtils.vec2mat(YEstViterbi(:, :, 1), params.m);
    selectedProb = probMap .* YEstOneHot(1:size(probMap, 1), :);
    % m + 1 x L
    barHandle = bar(start: start + L - 1, selectedProb', 1, 'stacked', 'FaceColor', 'flat');
    for j = 1:size(selectedProb, 1)
        barHandle(j).CData = cMap(j, :) * BARS_PLOT_DARKNESS_FACTOR;
    end

    % Plot Line
    plotHandle = plot(start: start + L - 1, probMap', 'LineWidth', 1.5);
    cellCMap = {};
    for i = 1:size(probMap, 1)
        cellCMap{i, 1} = cMap(i, :);
    end
    set(plotHandle, {'Color'}, cellCMap);

    % Viterbi bar
    % 1 x L
    YEstViterbi(:, YEstViterbi(:, :, 2) > 0, 1) = pwmVal;
    YEstMax(:, YEstMax(:, :, 2) > 0, 1) = pwmVal;
    % imagesc([start - .5, start + L - 1.5], [-3 * lowBarHeight / 4, -lowBarHeight / 4], ...
    %         [YEstViterbi(:, :, 1); YEstMax(:, :, 1)], [1, pwmVal]);
    imagesc([start - .5, start + L - 1.5], [-3 * lowBarHeight / 4, -lowBarHeight / 4], ...
            YEstViterbi(:, :, 1), [1, pwmVal]);
    ylim([-lowBarHeight, 1]);
    colormap(gca, cMapWithError);

    hold off;
end


function plotChromHMM(eid, chr, from, L, LOW_BAR_HEIGHT)
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
    % read ChromHMM bed
    to = from + L;
    bedFilePath = sprintf('../data/ChromHMM/%s_15_coreMarks_mnemonics.bed', eid);
    fid = fopen(bedFilePath);
    bedData = textscan(fid, '%s%d%d%d_%s', 'delimiter','\t');
    chrs = bedData{1};
    froms = bedData{2};
    tos = bedData{3};
    marks = bedData{4};
    markStrs = bedData{5};
    fclose(fid);

    % Find location
    mask = strcmp(chrs, chr) & (froms <= to) & (tos >= from);
    subFroms = froms(mask);
    subFroms(1) = from;
    subTos = tos(mask);
    subTos(end) = to;
    subMarks = marks(mask);
    diffs = subTos - subFroms;
    chromHMMmarks = repelem(subMarks, diffs)';
    % show image
    imagesc([from, to - 1], [-3 * LOW_BAR_HEIGHT / 4, -LOW_BAR_HEIGHT / 4], chromHMMmarks - 1, [0, 15]);

    % coloring
    colormap(gca, CHROMHMM_COLORMAP ./ 255);
    c = colorbar;
    set(c, 'Position', get(c, 'Position') + [0.05, 0, 0, 0]);
    c.Ticks = .5: 1: 15.5;
    c.TickLabels = MARKS_NAMES;
end


