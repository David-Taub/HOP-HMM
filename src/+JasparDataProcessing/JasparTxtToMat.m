function JasparTxtToMat()
    dbstop if error
    txtFilePath = '../data/Jaspar/raw/JASPAR_CORE_nonredundant_pfm_vertebrates.txt';
    % txtFilePath = 'data/Jaspar/raw/JASPAR_CORE_redundant_pfm_vertebrates.txt';
    % txtFilePath = 'data/Jaspar/raw/JASPAR_CORE_individual_pfm_vertebrates.txt';
    % txtFilePath = 'data/Jaspar/raw/JASPAR_CNE.txt';
    % txtFilePath = 'data/Jaspar/raw/JASPAR_OLD_PMWs.txt';

    % k x J x n
    [PWMs, lengths, names] = parseTxt(txtFilePath);

    % high values remove many, low values remove few
    % duplicatesToRemove = 0.10; strengthToRemove = 0.30; % 519 -> 26

    % n x J x k -> k x J x n
    % PWMs = permute(PWMs, [3, 2, 1]);
    % lengths = sum(sum(PWMs, 3) > 0, 2)';
    emptyMap = repmat(sum(PWMs, 3) == 0, [1,1,4]);
    % psaudo count
    PWMs = PWMs + 1;
    PWMs = bsxfun(@times, PWMs, 1 ./ sum(PWMs, 3));
    PWMs(emptyMap) = 0;
    % k x J x N -> k x n x J
    PWMs = permute(PWMs, [1, 3, 2]);

    out_filepath = '../data/Jaspar/PWMs.mat';
    save(out_filepath, 'PWMs', 'lengths', 'names');
    fprintf('Saved PWMs in %s\n', out_filepath)
end


% PWM - k x J x n
function [PWM, lengths, names] = parseTxt(txtFilePath)
    fileID = fopen(txtFilePath,'r');
    PWM = [];
    names = {};
    lengths = [];
    pwmId = 1;
    while not(feof(fileID))
        name = textscan(fileID,'>%*s%s',1,'Delimiter',{'\n', ' '}, 'MultipleDelimsAsOne', 1);
        names{pwmId} = name{1};
        for i=1:4
            textscan(fileID,'%*C [','Delimiter',{' ', '\n'}, 'MultipleDelimsAsOne', 1);
            counts = textscan(fileID,'%f','Delimiter',{' ', '\n'}, 'MultipleDelimsAsOne', 1);
            textscan(fileID, ' ]', 'Delimiter', {' ', '\n'}, 'MultipleDelimsAsOne', 1);
            counts = counts{1}';
            if length(counts)-size(PWM, 2) > 0
                PWM = cat(2, PWM, zeros(pwmId-1, length(counts)-size(PWM, 2), 4));
            elseif length(counts)-size(PWM, 2) < 0
                counts = cat(2, counts, zeros(1, size(PWM, 2)-length(counts)));
            end
            PWM(pwmId, :, i) = counts;
        end
        lengths(pwmId) = sum(sum(PWM(pwmId,:,:), 3) > 0, 2);
        pwmId = pwmId + 1;
    end
end