
% psi - N x m x k x L
% pX - N x 1
% gamma - N x m x L
function drawStatus(theta, params, gamma)%, alpha, beta, gamma, pX, xi)
    figure
    YsEst = mean(gamma, 3);
    YsEst = YsEst(:, 1);
    subplot(1,3,1);
    title('Gamma')
    for i = 1:params.m
        hold on;
        plot(permute(gamma(1, i, :), [3,1,2]));
    end
    legend(num2str([1:params.m]'));
    % [~, YsEst] = max(alpha, [], 2);
    % subplot(1,3,2);imagesc(permute(YsEst, [1,3,2])); colorbar;
    % % title('alpha')
    % [~, YsEst] = max(beta(:,:,1:end), [], 2);
    % subplot(1,3,3);imagesc(permute(YsEst, [1,3,2])); colorbar;
    % title('beta')

    subplot(1,3,2);
    for i = 1:params.m
        hold on;
        plot(exp(theta.E(i,:)));
    end
    ylim([0,1]);
    legend(num2str([1:params.m]'))
    title('E');

    subplot(1,3,3);
    for i = 1:params.m
        hold on;
        plot(exp(theta.G(i, :)));
    end

    ylim([0,max(exp(theta.G(:)))]);
    legend(num2str([1:params.m]'))
    title('G');
    drawnow
end
