% BaumWelchPWM.preprocessPWMs()
function preprocessPWMs()
    delete(fullfile('data', 'precomputation', 'pcPWMp.mat'));
	S = load('data/PWMsRaw.mat');
	PWM = S.PWMs;
	names = S.names;
	% DUPLICATES_FACTOR = 0.97; STRENGTH_FACTOR = 0.50; % 519 -> 26
	DUPLICATES_FACTOR = 0.1; STRENGTH_FACTOR = 0.1; % 519 -> 446
	% DUPLICATES_FACTOR = 0.97; STRENGTH_FACTOR = 0.97; % 519 -> 2
	% DUPLICATES_FACTOR = 0.97; STRENGTH_FACTOR = 0.95; % 519 -> 3
	% DUPLICATES_FACTOR = 0.97; STRENGTH_FACTOR = 0.91; % 519 -> 5
	% n x J x k -> k x J x n
	PWM = permute(PWM, [3, 2, 1]);
	lengths = sum(sum(PWM, 3) > 0, 2)';
	emptyMap = repmat(sum(PWM, 3) == 0, [1,1,4]);
	PWM = PWM + 1;
	PWM = bsxfun(@times, PWM, 1 ./ sum(PWM, 3));
	PWM(emptyMap) = 0;
	% k x J x N -> k x n x J
	PWM = permute(PWM, [1, 3, 2]);

	length(lengths)
	[PWM, lengths, names] = BaumWelchPWM.removedPWMsDuplicates(PWM, lengths, names, DUPLICATES_FACTOR);
	length(lengths)
	[PWM, lengths, names] = BaumWelchPWM.removePWMsWeak(PWM, lengths, names, STRENGTH_FACTOR);
	length(lengths)
	save('data/PWMs.mat', 'PWM', 'lengths', 'names');
end