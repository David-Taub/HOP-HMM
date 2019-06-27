% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
% dict = peaks.fasta2mem()
% seqs{i} = dict(chr).Data(from:to);
function dict = fasta2mem()
    HG19_FASTA_DIR = '../data/peaks/raw_genome';
    fprintf('Fasta -> mm\n');
    % save in dict opened hg19 fasta files as memory mapped files
    dict = containers.Map;
    fastaFiles = dir([HG19_FASTA_DIR, '/*.fa']);
    for i = 1:length(fastaFiles)
        if not(fastaFiles(i).isdir)
            fastaFilePath = fullfile(HG19_FASTA_DIR, fastaFiles(i).name);
            [~, chr, ~] = fileparts(fastaFilePath);
            mmFilePath = [fullfile(HG19_FASTA_DIR, chr) '.mm'];
            if not(exist(mmFilePath, 'file') == 2)
                fprintf('Converting file %s -> %s\n', fastaFilePath, mmFilePath);
                singleFasta2mm(fastaFilePath, mmFilePath);
            end
            dict(chr) = memmapfile(mmFilePath, 'format', 'uint8');
        end
    end
end

function singleFasta2mm(fastaFilePath, mmFilePath)
    fidIn = fopen(fastaFilePath,'r');
    fgetl(fidIn);
    fidOut = fopen(mmFilePath,'w');
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