% PWMs - k x J x n
function [PWMs, lengths, names] = removedPWMsDuplicates(PWMs, lengths, names, partToRemove)

    close all;
    fprintf('Removing Duplicates\n')
    k = size(PWMs, 1);
    dist = zeros(k, k);
    for i=1:k
        for j=1:i-1
            PWM1 = permute(PWMs(i, :, :), [2, 3, 1]);
            PWM2 = permute(PWMs(j, :, :), [2, 3, 1]);
            Similarity = conv2(PWM1, rot90(PWM2, 2));
            Similarity = max(Similarity(4,:));
            Similarity = Similarity ./ min(lengths(i), lengths(j));
            dist(i, j) = Similarity;
        end
    end
    % k x 1
    mostSimilar = max(dist, [], 2);
    % k x 1
    [mostSimilarSorted, inds] = sort(mostSimilar, 1);
    uniqueMask = true(k, 1);
    uniqueMask(inds(end-ceil(k*partToRemove):end)) = false;

    dd = dist(uniqueMask, :);
    dd = dd(:, uniqueMask);
    [vals, ind1] = max(dd, [], 1);
    [~, j] = max(vals, [], 2);
    i = ind1(j);
    figure
    subplot(1,2,1);imagesc(permute(exp(PWMs(i,:,:)), [2,3,1]))
    title('most similar PWM selected')
    subplot(1,2,2);imagesc(permute(exp(PWMs(j,:,:)), [2,3,1]))
    title('most similar PWM selected')
    PWMs = PWMs(uniqueMask, :, :);
    lengths = lengths(uniqueMask);
    names = names(uniqueMask);
end