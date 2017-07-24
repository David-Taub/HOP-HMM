cd /cs/stud/boogalla/projects/CompGenetics/mm9Genome
addpath('/cs/cbio/david/projects/CompGenetics/BaumWelch')

fprintf('Loading Genome\n');
load('/cs/cbio/tommy/Enhancers/Data/genome_mm9.mat');

fprintf('Loading histone peaks\n');
load('/cs/stud/boogalla/Work/data/mat/peaks_raw.mat');

fprintf('Loading background sequences\n');
seqsLength = 500;
negSeqs = readSeq('/cs/cbio/tommy/Enhancers/Data/NEnhancers.seq', seqsLength);

fprintf('Sorting background sequences by CG-content\n');
negSeqs = sortBaseContent(negSeqs);

% work on learner:
    % fprintf('Saving positive enhancers sequences\n');
    reader(T, genome, negSeqs, false, seqsLength);
    % fprintf('Loading positive enhancers sequences\n');
    load('/cs/stud/boogalla/projects/CompGenetics/mm9Genome/data/posSeqs.mat')
    size(posSeqs)
    % fprintf('Learning...\n');
    close all; learn(posSeqs, negSeqs, overlaps);
% work on best dataset:
while 1; reader(T, genome, negSeqs, true, seqsLength); end;
reader(T, genome, negSeqs, true, seqsLength);



addpath('C:\Users\booga\Dropbox\bio\projects\CompGenetics\mm9Genome')
addpath('C:\Users\booga\Dropbox\bio\projects\CompGenetics\BaumWelch')
load('C:\Users\booga\Dropbox\bio\projects\CompGenetics\mm9Genome\data\posSeqs.mat')

% [startT, T, E, ~, ~] = EM(posSeqs, m, 4, maxIter, tEpsilon, order);
[startT, T, E, ~, ~] = EM(posSeqs, 5, 4, 400, 0.05, 3);