
% assumes ranges are with same lengths
% lows / highs - 1 x N
function ranges = genRanges(lows, highs)
    n = cumsum([1;highs(:) - lows(:)]);
    z = ones(n(end)-1,1);
    z(n(1:end-1)) = [lows(1),lows(2:end)-highs(1:end-1)];
    rangesV = cumsum(z);
    ranges = reshape(rangesV, [highs(1) - lows(1), length(lows)]).';
end
