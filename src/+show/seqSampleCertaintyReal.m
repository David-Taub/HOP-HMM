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
    seqInd = randsample(N, 1);
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    % N x L x 2
    YEst = misc.viterbi(params, theta, dataset.X, dataset.pcPWMp);
    pwm_val = params.m + 1;
    LOW_BAR_HIEGHT = 0.1;

    subplot(3, 1, 3);
    ylim([-LOW_BAR_HIEGHT, 1]);
    YEstColored = YEst;
    YEstColored(seqInd, YEst(seqInd, :, 2) > 0, 1) = pwm_val;
    imagesc([.5, L-.5], [-3 * LOW_BAR_HIEGHT / 4, -LOW_BAR_HIEGHT / 4], ...
            [YEstColored(seqInd, :, 1); YEstColored(seqInd,:,1)], [1, pwm_val]);
    colormap(cMapWithError);
    text(L + 1, 0.5, 'Posterior Probability', 'FontSize', 10);
    text(L + 1, -LOW_BAR_HIEGHT / 2, 'Viterbi Estimation', 'FontSize', 10);
    ylabel(['Seq ', num2str(seqInd)]);
    xlabel('Position in Sequence');

    ax = gca;
    ax.YDir = 'normal';
    plotProbMap(params, posterior(seqInd, :, :), YEst(seqInd, :, 1), cMap);

    H3K27ac = readAllTissuesBedGraphs(dataset.chrs{seqInd}, dataset.starts(seqInd),...
                                      dataset.starts(seqInd) + L, tissueEIDs, 'H3K27ac');
    DNase = readAllTissuesBedGraphs(dataset.chrs{seqInd}, dataset.starts(seqInd),...
                                      dataset.starts(seqInd) + L, tissueEIDs, 'DNase');
    subplot(3, 1, 1);
    set(gca,'xtick',[]);
    title('Posterior & Viterbi Estimation');
    text(L + 1, 0.5, 'H3K27ac', 'FontSize', 10);
    plotProbMap(params, H3K27ac(:, :), YEst(seqInd, :, 1), cMap);

    subplot(3, 1, 2);
    set(gca,'xtick',[]);
    text(L + 1, 0.5, 'DNase', 'FontSize', 10);
    plotProbMap(params, DNase(:, :), YEst(seqInd, :, 1), cMap);


    legendStrings1 = strcat({'Enhancer Type '}, num2str([1:params.m - params.backgroundAmount]'));
    legendStrings2 = strcat({'Background '}, num2str([1:params.backgroundAmount]'));
    legendStrings = {legendStrings1{:}, legendStrings2{:}};
    legendStrings{pwm_val} = 'TFBS';
    legend(legendStrings);
    saveas(gcf, outpath);
end


% tracks - m x L
function tracks = readAllTissuesBedGraphs(chr, from, to, tissueEIDs, trackName)
    L = to - from;
    m = length(tissueEIDs);
    tracks = zeros(m, L)
    for i = 1:m
        EID = tissueEIDs{i};
        bedGraph = readBedGraph(EID, trackName);
        tracks(i, :) = getTrack(bedGraph, chr, from, to);
    end
end


% YEst - 1 x L
% probMap - 1 x m x L
function plotProbMap(params, probMap, YEst, cMap)
    BARS_PLOT_DARKNESS_FACTOR = 0.85;
    probMap = permute(probMap, [2, 3, 1]);
    L = size(YEst, 2);
    hold on;
    % m x L
    YEstOneHot = matUtils.vec2mat(YEst, params.m);
    p = plot(1:L, probMap', 'LineWidth', 1.5);
    cellCMap = {};
    for i = 1:params.m
        cellCMap{i, 1} = cMap(i, :);
    end

    set(p, {'Color'}, cellCMap);
    % m x L
    selectedProb = probMap .* YEstOneHot;
    % m + 1 x L
    b = bar(1:L, selectedProb', 1, 'stacked', 'FaceColor','flat');
    for j = 1:size(selectedProb, 1)
        b(j).CData = cMap(j, :) * BARS_PLOT_DARKNESS_FACTOR;
    end
    ylim([0, 1]);
    xlim([1, L]);
    hold off;
    yticks([0:0.2:1]);
    rotateYLabel();
end


function rotateYLabel()
    ylh = get(gca,'ylabel');
    gyl = get(ylh);
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation', 0, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
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


% row example:
% chr1    5113983 5113984 2.03288
function bedGraph = readBedGraph(EID, trackName)
    bedGraphFilePath = sprintf('../data/peaks/processed_bedgraphs/%s-%s.enh.bedgraph', EID, trackName);
    fprintf('Loading bed graph from %s\n', bedGraphFilePath);
    fid = fopen(bedGraphFilePath);

    bedGraphData = textscan(fid, '%s%d%d%f', 'delimiter','\t');
    bedGraph.chrs = bedGraphData{1};
    bedGraph.froms = bedGraphData{2};
    bedGraph.tos = bedGraphData{3};
    bedGraph.vals = bedGraphData{4};
    fclose(fid);
end


function ret = getTrack(bedGraph, trackChr, trackFrom, trackTo)
    mask = strcmp(bedGraph.chrs, trackChr) & (bedGraph.tos > trackFrom) & (bedGraph.froms < trackTo);
    % linear interp
    ret = interp1((bedGraph.tos(mask) + bedGraph.froms(mask)) / 2, bedGraph.vals(mask), trackFrom : trackTo);
    fprintf('got track %s: %d-%d\n', trackChr, trackFrom, trackTo);
end

