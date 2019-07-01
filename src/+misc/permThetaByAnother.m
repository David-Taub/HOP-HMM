% use the Hungarian Algorithm implementation downloaded from MathWorks.com
function theta = permThetaByAnother(params, thetaOrig, thetaEst)
    vectorizedOrig = misc.thetaToMat(params, thetaOrig, false);
    vectorizedEst = misc.thetaToMat(params, thetaEst, false);
    % perm = matUtils.repermuteMat(vectorizedOrig, vectorizedEst);
    perm = misc.munkres(pdist2(vectorizedOrig, vectorizedEst))';

    theta = permTheta(thetaEst, perm);
end


function theta = permTheta(theta, perm)
    theta.T = theta.T(perm, :);
    theta.T = theta.T(:, perm);
    theta.startT = theta.startT(perm);
    theta.G = theta.G(perm, :);
    theta.E(:, :) = theta.E(perm, :);
    % for i = 1:length(perm)
    %     theta.E(i, :) = theta.E(perm(i), :);
    % end
end
