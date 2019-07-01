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
