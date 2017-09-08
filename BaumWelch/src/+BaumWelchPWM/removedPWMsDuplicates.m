% PWMs - k x J x n
function [PWMs, lengths, names] = removedPWMsDuplicates(PWMs, lengths, names)
    % PART_TO_REMOVE = 1 / 15;
    PART_TO_REMOVE = 1 / 75;
    close all;
    k = size(PWMs, 1);
    dist = zeros(k, k);
    for i=1:k
        for j=i+1:k
            a = permute(PWMs(i, :, :), [2,3,1]);
            b = permute(PWMs(j, :, :), [2,3,1]);
            c = conv2(a, rot90(b, 2));
            dist(i, j) = max(c(4,:));
            dist(i, j) = dist(i, j) / min(lengths(i),lengths(j));
        end
    end
    d = sort(dist(:));
    thresh = d(end - round(k*k*PART_TO_REMOVE));
    uniqueMask = all(dist <= thresh, 2);

    dd = dist(uniqueMask, :);
    dd = dd(:, uniqueMask);
    [vals, ind1] = max(dd, [], 1);
    [~, j] = max(vals, [], 2);
    i = ind1(j);
    figure
    subplot(1,2,1);imagesc(permute(exp(PWMs(i,:,:)), [2,3,1]))
    subplot(1,2,2);imagesc(permute(exp(PWMs(j,:,:)), [2,3,1]))

    PWMs = PWMs(uniqueMask, :, :);
    lengths = lengths(uniqueMask);
    names = names(uniqueMask);
end