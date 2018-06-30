% PWMs - k x n x J
function [PWMs, lengths, names] = removePWMsWeak(PWMs, lengths, names, partsToRemove)
    fprintf('Removing "weak" background-like PWMs\n')
    % k x J
    vars = var(exp(PWMs), 0, 2);
    % k x 1
    varSums = sum(vars, 3);
    % k x 1
    % varMeans = varSums ./ lengths';
    varMeans = varSums;

    [~, inds] = sort(varMeans, 1);
    strongMap = inds(round(length(inds)*partsToRemove):end);

    figure
    i = inds(ceil(length(inds)*partsToRemove));
    imagesc(permute(exp(PWMs(i,:,:)), [2,3,1]));
    title('Weakest PWM selected')

    PWMs = PWMs(strongMap, :, :);
    lengths = lengths(strongMap);
    names = names(strongMap);
end
