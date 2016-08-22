
% assumes M - N x N
function peaks = makeBlockWise(peaks)
    N = 19;
    [~, ord] = sort(peaks(:, N+1));
    peaks(:, N+3) = (peaks(:, N+1) + peaks(:, N+2)) ./ 2;
    peaks = peaks(ord, :);
    i = 1;
    while i < size(peaks, 1)
        fprintf('\r%d / %d (%d)', i, size(peaks, 1), sum(peaks(i, 1:N), 2));
        if sum(peaks(i, 1:N), 2) > N
            fprintf('!\n');
        end
        center = (peaks(i+1, N+1) + peaks(i+1, N+2)) / 2;
        % [peaks(i, N+1) ,center , peaks(i, N+2)]
        if peaks(i, N+1) < center & center < peaks(i, N+2)
            peaks(i, N+1) = min(peaks(i, N+1), peaks(i+1, N+1));
            peaks(i, N+2) = max(peaks(i, N+2), peaks(i+1, N+2));
            peaks(i, 1:N) = peaks(i, 1:N) + peaks(i+1, 1:N);
            peaks(i + 1, :) = [];
        else
            i = i + 1;
        end
    end
end
function M = makeBlockWise2(M)
    N = length(M);
    reps = [];
    for i = 1 : 100
        subplot(1,2,1); imagesc(M);
        MW = M;
        MW = MW - triu(MW);
        MK = M';
        MW(MW>M') = MK(MW>M');
        MW(reps, :) = 0;
        [val, ind] = max(MW(:));
        [sr, c] = ind2sub([N, N], ind);
        tr = find(MW(:, c) > 0 & MW(:, c) < val, 1);
        [sr, tr]
        reps(end + 1) = tr;
        if sr <= tr
            continue;
        end
        M = repRowCol(M, tr, sr);

        subplot(1,2,2); imagesc(M);
        drawnow
        pause;
    end
    
end
% assumes ind1 < ind2
function M = repRowCol(M, ind1, ind2)
    M = M(:,[1:ind1-1, ind2, ind1+1:ind2-1, ind1, ind2+1:end]);
    M = M([1:ind1-1, ind2, ind1+1:ind2-1, ind1, ind2+1:end], :);
end
