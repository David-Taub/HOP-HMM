function showTheta(theta)
    % fprintf('Showing Theta\n')
    figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(1,4,1);
    imagesc(exp(theta.startT)); colorbar; title('startT')
    caxis([0,1])
    subplot(1,4,2);
    imagesc(exp(theta.T)); colorbar; title('T')
    subplot(1,4,3);
    imagesc(exp(theta.G)); colorbar; title('G')
    subplot(1,4,4);
    imagesc(exp(theta.E(:, :)))'; colorbar; title('E')
    drawnow
end