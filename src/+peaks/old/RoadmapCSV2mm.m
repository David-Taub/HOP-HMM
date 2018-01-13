% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
% [seqs, overlaps] = peaks.RoadmapCSV2mm();
function RoadmapCSV2mm()
    HG19_FASTA_DIR = 'data/peaks/raw/roadmap/HG19/uncompressed';
    fprintf('Fasta 2 mm\n');
    % save in dict opened hg19 fasta files as memory mapped files
    dict = containers.Map;
    fastaFiles = dir([HG19_FASTA_DIR, '/*.fa']);
    for i = 1:length(fastaFiles)
        if not(fastaFiles(i).isdir)
            fastaFilePath = fullfile(HG19_FASTA_DIR,fastaFiles(i).name)
            mmFilePath = fasta2MM(fastaFilePath);
        end
    end
end

function mmFilename = fasta2MM(FASTAfilename)
    fidIn = fopen(FASTAfilename,'r');
    fgetl(fidIn);
    [dirPath, filename, ~] = fileparts(FASTAfilename);
    mmFilename = [fullfile(dirPath, filename) '.mm'];
    fidOut = fopen(mmFilename,'w');
    newLine = sprintf('\n');
    blockSize = 2^20;
    while ~feof(fidIn)
        % Read in the data
        charData = fread(fidIn,blockSize,'*char')';
        % Remove new lines
        charData = strrep(charData,newLine,'');
        % Convert to integers
        intData = nt2int(charData);
        % Write to the new file
        fwrite(fidOut,intData,'uint8');
    end
    fclose(fidIn);
    fclose(fidOut);
end
