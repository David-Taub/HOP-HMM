function genPosNegSeqs(N1, N2, L)

    [PWMs, lengths, names] = misc.PWMs();
    % [k, n, J] = size(PWMs);
    peakedPWM = 130;
    names(peakedPWM)
    overlaps = zeros(N1 + N2, 2);
    overlaps(1:N1, 1) = 1;
    overlaps(N1+1:end, 2) = 1;

    seqs = mnrnd(1, [0.2573, 0.2433, 0.2426, 0.2568], (N1 + N2) * L) * [1,2,3,4]';
    seqs = reshape(seqs, [N1 + N2, L]);
    for i = 1:N1
        orientation = rand(1) > 0.5;
        % startIndex = randi(L - lengths(peakedPWM)-1);
        startIndex = round(L/2);
        for j = 1:lengths(peakedPWM)
            if orientation
                seqs(i, startIndex + j) = mnrnd(1, PWMs(peakedPWM, :, j)) * [1,2,3,4]';
            else
                seqs(i, startIndex - j) = 5-mnrnd(1, PWMs(peakedPWM, :, j)) * [1,2,3,4]';
            end
        end
    end
    lengths = ones(N1+N2, 1) * L;
    save('data/peaks/randGenerated/mergedPeaksMinimized.mat', 'seqs', 'overlaps', 'lengths');
end