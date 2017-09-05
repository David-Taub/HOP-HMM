
function [PWMs, lengths, names] = removedPWMsDuplicates(PWMs, lengths, names)
    % [PWMs, lengths, names] = BaumWelchPWM.PWMs();
    REMOVE_RATIO = 1 / 20;
    close all;
    k = size(PWMs, 1);
    dist = zeros(k, k);
    for i=1:k
        for j=i+1:k
            a = permute(PWMs(i, :, :), [2,3,1]);
            b = permute(PWMs(j, :, :), [2,3,1]);
            c = conv2(a, rot90(b, 2));
            dist(i, j) = max(c(:));
            dist(i, j) = dist(i, j) / min(lengths(i),lengths(j));
        end
        i
    end
    % dist = dist - diag(diag(dist));

    % [distOrdered, order] = matUtils.clustMatRows(dist);
    % distOrdered = distOrdered(:, order);

    d = sort(dist(:));
    thresh = d(end - round(k*k*REMOVE_RATIO));
    % distC = dist;
    % distC() = nan;

    % while any(isnan(distC(:)))
    %     for i=1:size(distC, 1)
    %         if any(isnan(distC(i, :)),2)
    %             distC(i,:) = [];
    %             distC(:,i) = [];
    %             lengths(i) = [];
    %             PWMs(i, :, :) = [];
    %             names{i} = [];
    %             break;
    %         end
    %     end
    % end
    similarMask = any(dist > thresh, 2);
    PWMs = PWMs(~similarMask, :, :);
    lengths = lengths(~similarMask);
    names = names(~similarMask);

    % figure;
    % imagesc(dist); colorbar;

    % figure;
    % imagesc(distC); colorbar;
end