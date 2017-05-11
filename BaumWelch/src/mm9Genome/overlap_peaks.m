% function [SS2,names] = overlap_peaks()
function peakFilePaths = getPeakFilePaths()
    baseDir = '../data/Tommy/GSE29184_Bing_Ren';
    peakFilePaths={};
    dirs = dir();
    for i=1:length(dirs),
        dirName = dirs(i).name;
        if ~isdir(dirName)
            continue;
        end
        peakMatFiles=dir(fullfile(baseDir, dirName, '*k27ac*.peak.mat'));
        for j=1:length(peakMatFiles)
            peakFilePaths{end+1}=fullfile(baseDir, dirName, peakMatFiles(j).name);
        end;
    end;
end
clear dirs i j d f


names = {};
for j=1:length(peakFilePaths),
    % save peaks
    A=load(peakFilePaths{j});
    sequences{j}=A.S;
    % save names
    ii=strfind(peakFilePaths{j},'/')-1; names{j}=peakFilePaths{j}(1:ii);
end

for j=1:length(peakFilePaths),
    tic;
    S=sequences{j};
    ov=NaN*ones(length(peakFilePaths),length(S),'single');

    % over versus all
    parfor i=1:length(peakFilePaths),
	ov(i,:) = overlap(S,sequences{i},50);
    end

    % save
    for i=1:length(S),
	S{i}.overlap=ov(:,i)';
    end
    SS2{j}=S;
end

for j=1:length(peakFilePaths),
    S = SS2{j};
    % S = annotate_peaks(SS2{j});
    mname = [names{j} '-H3K27ac.peaks.mat'];
    save(mname, 'S', 'names');
end
