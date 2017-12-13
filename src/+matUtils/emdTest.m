
function [dist] = emdTest(x, y)
    nbins = 60;
    [xh, xEdges] = histcounts(x, nbins, 'Normalization', 'probability');
    [yh, yEdges] = histcounts(y, nbins, 'Normalization', 'probability');
    xCenters = (xEdges(1:end-1) + xEdges(2:end)) / 2;
    yCenters = (yEdges(1:end-1) + yEdges(2:end)) / 2;
    [~, dist] = matUtils.emd(xh', yh', xCenters', yCenters', @matUtils.gdf);
end