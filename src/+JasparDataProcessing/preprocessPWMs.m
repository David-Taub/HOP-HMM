% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
% JasparDataProcessing.preprocessPWMs()
function preprocessPWMs()
	dbstop if error
    delete(fullfile('..', 'data', 'precomputation', 'pcPWMp.mat'));
    txtFilePath = '../data/Jaspar/raw/JASPAR_CORE_nonredundant_pfm_vertebrates.txt';
    % txtFilePath = 'data/Jaspar/raw/JASPAR_CORE_redundant_pfm_vertebrates.txt';
    % txtFilePath = 'data/Jaspar/raw/JASPAR_CORE_individual_pfm_vertebrates.txt';
    % txtFilePath = 'data/Jaspar/raw/JASPAR_CNE.txt';
    % txtFilePath = 'data/Jaspar/raw/JASPAR_OLD_PMWs.txt';

	% k x J x n
    [PWM, lengths, names] = parseTxt(txtFilePath);

	% high values remove many, low values remove few
	DUPLICATES_TO_REMOVE = 0.10; STRENGTH_TO_REMOVE = 0.30; % 519 -> 26

	% n x J x k -> k x J x n
	% PWM = permute(PWM, [3, 2, 1]);
	% lengths = sum(sum(PWM, 3) > 0, 2)';
	emptyMap = repmat(sum(PWM, 3) == 0, [1,1,4]);
	PWM = PWM + 1;
	PWM = bsxfun(@times, PWM, 1 ./ sum(PWM, 3));
	PWM(emptyMap) = 0;
	% k x J x N -> k x n x J
	PWM = permute(PWM, [1, 3, 2]);

	length(lengths)
	[PWM, lengths, names] = JasparDataProcessing.removedPWMsDuplicates(PWM, lengths, names, DUPLICATES_TO_REMOVE);
	length(lengths)
	[PWM, lengths, names] = JasparDataProcessing.removePWMsWeak(PWM, lengths, names, STRENGTH_TO_REMOVE);
	length(lengths)
	save('../data/Jaspar/PWMs.mat', 'PWM', 'lengths', 'names');
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
			textscan(fileID,' ]','Delimiter',{' ', '\n'}, 'MultipleDelimsAsOne', 1);
			counts = counts{1}';
			if length(counts)-size(PWM, 2) > 0
				PWM = cat(2, PWM, zeros(pwmId-1, length(counts)-size(PWM, 2), 4));
			elseif length(counts)-size(PWM, 2) < 0
				counts = cat(2, counts, zeros(1, size(PWM, 2)-length(counts)));
			end
			PWM(pwmId, :, i) = counts;
		end
		lengths(pwmId) = sum(sum(PWM(pwmId,:,:), 3) > 0, 2)
		pwmId = pwmId + 1;
	end
end