% bedGraphs - tissues x tracks
function bedGraphs = readAllBedGraphs(tissueEIDs, trackNames)
    tissueWithTracksCount = length(tissueEIDs);
    for trackInd = 1:length(trackNames)
        for i = 1:tissueWithTracksCount
            EID = tissueEIDs{i};
            if startsWith(EID, 'E')
                bedGraphs{i, trackInd} = readBedGraph(EID, trackNames{trackInd});
            end
        end
    end
end


% row example:
% chr1    5113983 5113984 2.03288
function bedGraph = readBedGraph(EID, trackName)
    bedGraphFilePath = sprintf('../data/peaks/processed_bedgraphs/%s-%s.enh.u.bedgraph', EID, trackName);
    % bedGraphFilePath = sprintf('../data/peaks/raw_bedgraphs/%s-%s.pval.bedgraph', EID, trackName);
    fprintf('Loading bed graph from %s\n', bedGraphFilePath);
    fid = fopen(bedGraphFilePath);

    bedGraphData = textscan(fid, '%s%d%d%f', 'delimiter','\t');
    bedGraph.chrs = bedGraphData{1};
    bedGraph.froms = bedGraphData{2};
    bedGraph.tos = bedGraphData{3};
    bedGraph.vals = bedGraphData{4};
    bedGraph.centers = (bedGraph.froms + bedGraph.tos) ./ 2;
    fclose(fid);
end
