% sample sequences, and draw for each colorful plots with what the posterior
% probability was compared to the correct state per letter
function seqSampleCertaintyReal(params, theta, dataset, outpath, tissueEIDs)
    [N, L] = size(dataset.X);
    % gamma - N x m x L
    % psi - N x m x k x L
    [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, dataset.X, dataset.pcPWMp);
    % N x m x L
    posterior = calcPosterior(params, gamma, psi);
    cMap = lines(params.m);

    PWM_COLOR = [0, 0, 0];
    cMapWithError = [cMap;  PWM_COLOR];
    seqInd = 1;
    trackNames = {'H3K27ac', 'DNase'};
    bedGraphs = misc.readAllBedGraphs(tissueEIDs, trackNames);
    YEst = misc.viterbi(params, theta, dataset.X, dataset.pcPWMp);
    for i = 1:10
        [tracks, seqInd] = getSeqWithTracks(dataset, bedGraphs, seqInd + 1);
        H3K27acTrack = tracks(:, :, 1);
        DNaseTrack = tracks(:, :, 2);
        outpathI = sprintf('%s.%s.jpg', outpath, i);
        figure('units', 'pixels', 'Position', [0 0 1000 1000]);
        % N x L x 2
        pwm_val = params.m + 1;
        LOW_BAR_HIEGHT = 0.1;

        subplot(3, 1, 3);
        hold on;
        ylim([-LOW_BAR_HIEGHT, 1]);
        plotProbMap(params, permute(posterior(seqInd, :, :), [2, 3, 1]), YEst(seqInd, :, 1), cMap);
        ylim([-LOW_BAR_HIEGHT, 1]);
        YEstColored = YEst;
        YEstColored(seqInd, YEst(seqInd, :, 2) > 0, 1) = pwm_val;
        imagesc([.5, L-.5], [-3 * LOW_BAR_HIEGHT / 4, -LOW_BAR_HIEGHT / 4], ...
                [YEstColored(seqInd, :, 1); YEstColored(seqInd,:,1)], [1, pwm_val]);
        ylim([-LOW_BAR_HIEGHT, 1]);
        colormap(cMapWithError);
        text(L + 1, 0.5, 'Posterior Prob.', 'FontSize', 10);
        text(L + 1, -LOW_BAR_HIEGHT / 2, 'Viterbi Est.', 'FontSize', 10);
        ylabel(['Seq ', num2str(seqInd)]);
        xlabel('Position in Sequence');
        ylabel('P(y_{t}|X)');
        ax = gca;
        ax.YDir = 'normal';
        legendStrings1 = strcat({'Enhancer Type '}, num2str([1:params.m - params.backgroundAmount]'));
        legendStrings2 = strcat({'Background '}, num2str([1:params.backgroundAmount]'));
        legendStrings = {legendStrings1{:}, legendStrings2{:}};
        legendStrings{pwm_val} = 'TFBS';
        legend(legendStrings);

        subplot(3, 1, 1);
        set(gca,'xtick',[]);
        title('Posterior & Viterbi Estimation');
        text(L + 1, 0.5, 'H3K27ac', 'FontSize', 10);
        hold on;
        ylim([0, max(H3K27acTrack(:))]);
        plotProbMap(params, H3K27acTrack, YEst(seqInd, :, 1), cMap);
        ylim([0, max(H3K27acTrack(:))]);
        ylabel('-log_{10}(p-value)');

        subplot(3, 1, 2);
        set(gca,'xtick',[]);
        text(L + 1, 0.5, 'DNase', 'FontSize', 10);
        hold on;
        ylim([0, max(DNaseTrack(:))]);
        plotProbMap(params, DNaseTrack, YEst(seqInd, :, 1), cMap);
        ylim([0, max(DNaseTrack(:))]);
        ylabel('-log_{10}(p-value)');
        saveas(gcf, outpath);
    end
end


function [tracks, seqInd] = getSeqWithTracks(dataset, bedGraphs, startInd)
    N = size(dataset.X, 1);
    for seqInd = startInd:N
        % fprintf('Trying sequence %d\n', seqInd)
        tracks = getTracks(dataset, bedGraphs, seqInd);
        if all(all(any(tracks > 0, 2), 1), 3)
            fprintf('Found seq %d\n', seqInd);
            return;
        end
    end
end


function tracks = getTracks(dataset, bedGraphs, seqInd)
    L = size(dataset.X, 2);
    chr = dataset.chrs{seqInd};
    to = dataset.starts(seqInd) + L;
    from = dataset.starts(seqInd);
    tracks = zeros(size(bedGraphs, 1), L, size(bedGraphs, 2));
    for trackInd = 1:size(bedGraphs, 2)
        for i = 1:size(bedGraphs, 1)
            fprintf('looking for track index %d, for tissue index %d\n', trackInd, i)
            tracks(i, :, trackInd) = getTrack(bedGraphs{i, trackInd}, chr, from, to);
            if all(tracks(i, :, trackInd) == 0)
                fprintf('Seq has not enough samples\n', trackInd)
                return;
            end
        end
    end
end


% YEst - 1 x L
% probMap - m x L
function plotProbMap(params, probMap, YEst, cMap)
    BARS_PLOT_DARKNESS_FACTOR = 0.85;
    % probMap = permute(probMap, [2, 3, 1]);
    L = size(YEst, 2);
    hold on;
    if size(probMap, 1) == params.m
        % m x L
        YEstOneHot = matUtils.vec2mat(YEst, params.m);
        selectedProb = probMap .* YEstOneHot(1:size(probMap, 1), :);
        % m + 1 x L
        b = bar(1:L, selectedProb', 1, 'stacked', 'FaceColor','flat');
        for j = 1:size(selectedProb, 1)
            b(j).CData = cMap(j, :) * BARS_PLOT_DARKNESS_FACTOR;
        end
    end

    p = plot(1:L, probMap', 'LineWidth', 1.5);
    cellCMap = {};
    for i = 1:size(probMap, 1)
        cellCMap{i, 1} = cMap(i, :);
    end
    set(p, {'Color'}, cellCMap);

    % m x L
    xlim([1, L]);
    rotateYLabel();
end


function rotateYLabel()
    ylh = get(gca,'ylabel');
    gyl = get(ylh);
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation', 0, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
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
    % linear interp
    foundSamples = sum(mask);
    % fprintf('Found %d sample points in (%s:%d-%d)\n', foundSamples, trackChr, trackFrom, trackTo);
    if foundSamples < MIN_SAMPLES_IN_SEQ
        ret = zeros(trackTo - trackFrom, 1);
        return
    end
    knownPoints = double([(bedGraph.tos(mask) + bedGraph.froms(mask)) / 2]');
    knownVals = bedGraph.vals(mask)';
    wantedPoints = double([trackFrom : trackTo - 1]);
    ret = interp1(knownPoints, knownVals, wantedPoints);
end

