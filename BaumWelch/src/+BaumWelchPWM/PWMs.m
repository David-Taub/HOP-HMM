% k x n x J
function [PWM, lengths, names] = PWMs()
	S = load('data/PWMs.mat');
	% k x n x J
	PWM = S.PWM;
	names = S.names;
	% k x n x J -> k x 1
	lengths = sum(sum(PWM, 2) > 0, 3)';
	emptyMap = repmat(sum(PWM, 2) == 0, [1,4,1]);
	PWM = PWM + 1;
	% k x n x J -> k x J x n
	PWM = permute(PWM, [1, 3, 2]);
	PWM = bsxfun(@times, PWM, 1 ./ sum(PWM, 3));
	PWM = permute(PWM, [1, 3, 2]);
	PWM(emptyMap) = 0;
	% k x J x N -> k x n x J
end