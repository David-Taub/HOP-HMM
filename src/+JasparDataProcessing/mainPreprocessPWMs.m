
function selectedPWMs = mainPreprocessPWMs(duplicatesToRemove, strengthToRemove, mergedPeaksMin, tissueList, k)
	[PWMs, lengths, names] = misc.PWMs();
	selectedPWMs = PWMsFeatureSelect(mergedPeaksMin, tissueList, PWMs, lengths);
	PWMs = PWMs(selectedPWMs, :, :);
	lengths = lengths(selectedPWMs);
	names = names(selectedPWMs);
	length(lengths)
	if duplicatesToRemove > 0
		[PWMs, lengths, names] = JasparDataProcessing.removedPWMsDuplicates(PWMs, lengths, names, duplicatesToRemove);
		length(lengths)
	end
	if strengthToRemove > 0
		[PWMs, lengths, names] = JasparDataProcessing.removePWMsWeak(PWMs, lengths, names, strengthToRemove);
		length(lengths)
	end
	k = min(k, size(lengths));
	PWMs = PWMs(1:k, :, :);
	lengths = lengths(1:k);
	names = names(1:k);
	out_filepath = '../data/precomputation/SelectedPWMs.mat';
	save(out_filepath, 'PWMs', 'lengths', 'names');
	fprintf('Saved PWMs in %s\n', out_filepath)
end
