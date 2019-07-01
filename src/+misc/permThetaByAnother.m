% use the Hungarian Algorithm implementation downloaded from MathWorks.com
function theta = permThetaByAnother(params, thetaOrig, thetaEst)
    vectorizedOrig = misc.thetaToMat(params, thetaOrig, false);
    vectorizedEst = misc.thetaToMat(params, thetaEst, false);
    % perm = matUtils.repermuteMat(vectorizedOrig, vectorizedEst);
    perm = misc.munkres(pdist2(vectorizedOrig, vectorizedEst))';
    theta = misc.permTheta(thetaEst, perm);
end

