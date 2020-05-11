function showTheta(theta)
    % fprintf('Showing Theta\n')
    figure('units', 'pixels', 'Position', [0 0 1500 1000]);

    subplot(1,4,1);
    imagesc(exp(theta.startT)); colorbar; title('startT');
    addText(exp(theta.startT));
    % caxis([0,1])

    subplot(1,4,2);
    imagesc(exp(theta.T)); colorbar; title('T');
    addText(exp(theta.T));

    subplot(1,4,3);
    imagesc(exp(theta.G)); colorbar; title('G');

    subplot(1,4,4);
    imagesc(exp(theta.E(:, :)))'; colorbar; title('E');
    drawnow
end


function addText(A)
    hold on;
    A(A > 0.5) = 1 - A(A > 0.5);
    A = round(1 ./ A);
    x = repmat(1:size(A,2), size(A,1), 1); % generate x-coordinates
    y = repmat([1:size(A,1)]', 1, size(A,2)); % generate y-coordinates
    t = num2cell(A); % extact values into cells
    t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
    text(x(:), y(:), t, 'HorizontalAlignment', 'Center')
    hold off;
end