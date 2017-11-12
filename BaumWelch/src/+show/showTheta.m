function showTheta(theta)
    fprintf('Showing Theta\n')
    figure
    subplot(1,4,1);
    imagesc(exp(theta.startT)); colorbar; title('startT')
    subplot(1,4,2);
    imagesc(exp(theta.T) - diag(diag(exp(theta.T)))); colorbar; title('T')
    subplot(1,4,3);
    imagesc(exp(theta.G)); colorbar; title('G')
    subplot(1,4,4);
    plot(exp(theta.E(:))); colorbar; title('E')
    drawnow
end