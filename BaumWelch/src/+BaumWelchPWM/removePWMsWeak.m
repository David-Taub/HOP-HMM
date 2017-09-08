% PWMs - k x n x J
function [PWMs, lengths, names] = removePWMsWeak(PWMs, lengths, names)
    PART_TO_REMOVE = 0.97;
    % k x J
    vars = var(exp(PWMs), 0, 2);
    % k x 1
    varSums = sum(vars, 3);
    % k x 1
    % varMeans = varSums ./ lengths';
    varMeans = varSums;

    [~, inds] = sort(varMeans, 1);
    strongMap = inds(round(length(inds)*PART_TO_REMOVE):end);

    figure
    i = inds(round(length(inds)*PART_TO_REMOVE));
    imagesc(permute(exp(PWMs(i,:,:)), [2,3,1]));

    PWMs = PWMs(strongMap, :, :);
    lengths = lengths(strongMap);
    names = names(strongMap);
end
