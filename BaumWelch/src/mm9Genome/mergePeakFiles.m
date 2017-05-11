% merge all peaks.mat from /cs/cbio/tommy/Enhancers/Data/GSE29184_Bing_Ren
% to a single mat file
function mergePeakFiles()
    totalpeaks
    peaksBasePath = '../data/peaks';
    peakFiles = dir(peaksBasePath);
    for i = 1:length(peakFiles)
        if peakFiles(i).isdir
            continue
        end
        peaks = load(fullfile(peaksBasePath, peakFiles(i).name));
        S2 = [peaks.S{:}]
    end
end