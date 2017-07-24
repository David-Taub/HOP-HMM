% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
% outputPath = fullfile('data', 'outMainPWMCor.mat');
% load(outputPath)
% tempMain(overlaps, maxPeaks)
function tempMain(overlaps, maxPeaks)
    close all
    % checkRegression(overlaps, maxPeaks);
    [aucs, regAuc] = aucMat(overlaps, maxPeaks);
    save('/cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src/data/temp', 'aucs', 'regAuc')
end


function loss = trainPredict(XTrain, YTrain, XTest, YTest)
    mdl = fitcecoc(XTrain, YTrain);
    YTestEst = predict(mdl, XTest);
    accuricy = sum(YTestEst == YTest) / length(YTest);
    fprintf('%d x %d - %.2f \n', size(XTrain, 1), size(XTrain, 2), accuricy);
    loss = 1 - accuricy;
end

function [aucs, regAuc] = aucMat(overlaps, maxPeaks)
    testRatio = 0.25;
    T = 15;

    % chosen = sequentialfs(@trainPredict, X(sampleMask, :), Y(sampleMask))
    % % trainPredict(, XTest(:, chosen), YTest)
    % mdl = fitcecoc(XTrain(:, chosen), YTrain);
    % [YTestEst, score, cost] = predict(mdl, XTest(:, chosen));

    [N, r] = size(overlaps);
    [N2, k] = size(maxPeaks);
    assert(N2 == N);
    testMask = rand(N,1) < testRatio;
    X = maxPeaks;
    aucs = zeros(r, r, k);
    regAuc = zeros(r, r);
    figure
    for i = 1:r
        for j = i+1:r
            maskIJ = sum(overlaps(:, [i, j]) > 0, 2) == 1;
            overlapsTrain = overlaps((~testMask) & maskIJ, [i,j]);
            overlapsTest = overlaps(testMask & maskIJ, [i,j]);
            YTrain = (overlapsTrain(:, 1) == 0) + 1;
            YTest = (overlapsTest(:, 1) == 0) + 1;
            for t = 1:k
                XTrainT = X((~testMask) & maskIJ, t);
                auc = matUtils.getAucRoc(XTrainT(YTrain == 1), XTrainT(YTrain == 2), false);
                aucs(i, j, t) = auc; aucs(j, i, t) = auc;
                fprintf('%d,%d: %d - %.2f\n', i, j, t, auc);
                if mod(t, 40) == 0
                    [aucsSorted, ind] = sort(aucs, 3, 'descend');
                    subplot(1,2,2);imagesc(aucsSorted(:,:,1)); colorbar;
                    fprintf('%d,%d: %d - %.2f Best so far: %.2f (%d)\n', i, j, t, auc, aucsSorted(i,j,1), ind(i,j,1));
                    drawnow
                end
            end
            [~, ind] = sort(aucs, 3, 'descend');
            XTest = X(testMask & maskIJ, ind(i,j,1:T));
            XTrain = X((~testMask) & maskIJ, ind(i,j,1:T));
            [B,~,~] = mnrfit(XTrain, YTrain);
            pihat = mnrval(B, XTest);
            pos = pihat(YTest == 1, 1);
            neg = pihat(YTest == 2, 1);

            regAuc(i,j) = matUtils.getAucRoc(pos, neg, true);
            regAuc(j,i) = regAuc(i,j);
            drawnow
        end
    end
    figure
    imagesc(regAuc(:,:))
    addTicks();
    title(['AUC ROC of tissue specific enhancers by sum of 10\\5\\3\\1 best PSSM in each sequence of ',int2str(T),' best PWM'])
end

function addTicks()
    tissue_names = names();
    f = gca;
    f.XTick = 1:length(tissue_names);
    f.YTick = 1:length(tissue_names);
    f.XTickLabel = tissue_names;
    f.YTickLabel = tissue_names;
    f.XTickLabelRotation = 45;
end
function out = names()
    out = ...
    {'E003 H1 Cells',...
     'E004 H1 BMP4 Derived Mesendoderm Cultured Cells',...
     'E005 H1 BMP4 Derived Trophoblast Cultured Cells',...
     'E006 H1 Derived Mesenchymal Stem Cells',...
     'E007 H1 Derived Neuronal Progenitor Cultured Cells',...
     'E008 H9 Cells',...
     'E017 IMR90 fetal lung fibroblasts Cell Line',...
     'E021 iPS DF 6.9 Cells',...
     'E022 iPS DF 19.11 Cells',...
     'E029 Primary monocytes from peripheral blood',...
     'E032 Primary B cells from peripheral blood',...
     'E034 Primary T cells from peripheral blood',...
     'E046 Primary Natural Killer cells from peripheral blood',...
     'E050 Primary hematopoietic stem cells G-CSF-mobilized Female',...
     'E055 Foreskin Fibroblast Primary Cells skin01',...
     'E056 Foreskin Fibroblast Primary Cells skin02',...
     'E059 Foreskin Melanocyte Primary Cells skin01',...
     'E080 Fetal Adrenal Gland',...
     'E084 Fetal Intestine Large',...
     'E085 Fetal Intestine Small',...
     'E089 Fetal Muscle Trunk',...
     'E090 Fetal Muscle Leg',...
     'E091 Placenta',...
     'E092 Fetal Stomach',...
     'E093 Fetal Thymus',...
     'E094 Gastric',...
     'E097 Ovary',...
     'E098 Pancreas',...
     'E100 Psoas Muscle',...
     'E109 Small Intestine',...
     'E114 A549 EtOH 0.02pct Lung Carcinoma Cell Line',...
     'E116 GM12878 Lymphoblastoid Cells',...
     'E117 HeLa-S3 Cervical Carcinoma Cell Line',...
     'E118 HepG2 Hepatocellular Carcinoma Cell Line',...
     'E119 HMEC Mammary Epithelial Primary Cells',...
     'E120 HSMM Skeletal Muscle Myoblasts Cells',...
     'E121 HSMM cell derived Skeletal Muscle Myotubes Cells',...
     'E122 HUVEC Umbilical Vein Endothelial Primary Cells',...
     'E123 K562 Leukemia Cells',...
     'E124 Monocytes-CD14+ RO01746 Primary Cells',...
     'E125 NH-A Astrocytes Primary Cells',...
     'E126 NHDF-Ad Adult Dermal Fibroblast Primary Cells',...
     'E127 NHEK-Epidermal Keratinocyte Primary Cells',...
     'E128 NHLF Lung Fibroblast Primary Cells'};
end

function checkRegression(overlaps, maxPeaks)
    testRatio = 0.1;
    [N, r] = size(overlaps);
    Y = overlaps > 0;
    Y = Y * [1:r]';
    X = maxPeaks;
    testMask = rand(N,1) < testRatio;
    XTest = X(testMask, :);
    XTrain = X(~testMask, :);
    YTest = Y(testMask, :);
    YTrain = Y(~testMask, :);
    fprintf('feature selection:\n ');

    sampleMask = rand(N,1) < 0.5;
    sum(sampleMask)
    chosen = sequentialfs(@trainPredict, X(sampleMask, :), Y(sampleMask))
    % trainPredict(, XTest(:, chosen), YTest)
    mdl = fitcecoc(XTrain(:, chosen), YTrain);
    [YTestEst, score, cost] = predict(mdl, XTest(:, chosen));


    figure
    subplot(1,2,1);imagesc(overlaps(testMask, :)); colorbar;
    xlabel('Cell Type'); ylabel('Sequences');
    title('Overlaps (height of H3k27ac peak)');
    subplot(1,2,2);imagesc(score); colorbar;
    xlabel('Cell Type'); ylabel('Sequences');
    title('Estimation');


end