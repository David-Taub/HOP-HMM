cd /cs/stud/boogalla/projects/CompGenetics/mm9Genome

fprintf('Loading Genome\n');
load('/cs/cbio/tommy/Enhancers/Data/genome_mm9.mat');
fprintf('Loading histone peaks\n');
load('/cs/stud/boogalla/Work/data/mat/peaks_raw.mat');
fprintf('Loading background sequences\n');
negSeqs = readSeq('/cs/cbio/tommy/Enhancers/Data/NEnhancers.seq', 500);
fprintf('Sorting background sequences by CG-content\n');
negSeqs = sortBaseContent(negSeqs);

% work on learner:
fprintf('Saving positive enhancers sequences\n');
reader(T, genome, negSeqs, false);
fprintf('Loading positive enhancers sequences\n');
load('/cs/stud/boogalla/projects/CompGenetics/mm9Genome/data/posSeqs.mat')

fprintf('Learning...\n');
learn(posSeqs, negSeqs, overlaps); 
% work on best dataset:
% while 1; reader(T, genome, negSeqs, true); end;
