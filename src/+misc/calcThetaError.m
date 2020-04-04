
function rmse = calcThetaError(params, thetaOrig, thetaEst)
    thetaOrigMat = misc.thetaToMat(params, thetaOrig, false);
    thetaEstMat = misc.thetaToMat(params, thetaEst, false);
    rmse = sqrt(mean((exp(thetaOrigMat(:)) - exp(thetaEstMat(:))) .^ 2, 1));
end