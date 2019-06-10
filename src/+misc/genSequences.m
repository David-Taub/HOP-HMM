
function [X, Y] = genSequences(theta, params, N, L)
    fprintf('Generating sequences... \n');
    params.order = length(size(theta.E)) - 1;
    Y = zeros(N, L, 2);
    X = zeros(N, L);
    E = exp(theta.E);
    T = exp(theta.T);
    G = exp(theta.G);
    startT = exp(theta.startT);
    for j = 1:N
        % first letter
        Etemp = matUtils.sumDim(E, 2 : params.order);
        Etemp = bsxfun(@times, Etemp, 1 ./ sum(Etemp, length(size(Etemp))));
        Y(j, 1, 1) = sampleMultiVarDist(startT);
        Y(j, 1, 2) = 0;
        X(j, 1) = sampleMultiVarDist(Etemp(Y(j, 1), :)');
        fprintf('Sequence index %d / %d: %d ', j, N, Y(j, 1, 1));
        t = 2;
        while t <= L
            yt = Y(j, t-1, 1);
            state = sampleMultiVarDist([T(yt, :), G(yt, :)]');

            if params.m < state
                % PWM step
                motif = state - params.m;
                % fprintf('[%d]', motif)
                for i=1:params.lengths(motif)
                    if mod(t, 20) == 0
                        fprintf('*');
                    end
                    X(j, t) = sampleMultiVarDist(params.PWMs(motif, :, i)');
                    % fprintf('%d',X(j, t))
                    Y(j, t, 1) = yt;
                    Y(j, t, 2) = motif;
                    assert(X(j, t) > 0);
                    t = t + 1;
                    if t > L
                        % sequence too long
                        break
                    end
                end
                if t <= L
                    if mod(t, 20) == 0
                        fprintf('%d', state);
                    end
                    X(j, t) = emitBaseState(X, params, E, t, yt, j);
                    Y(j, t, 1) = yt;
                    Y(j, t, 2) = 0;
                    t = t + 1;
                end
            else
                if mod(t, 20) == 0
                    fprintf('%d', state);
                end
                ytNext = state;
                X(j, t) = emitBaseState(X, params, E, t, ytNext, j);
                if ytNext ~= yt
                    fprintf('(%d>%d)', yt, ytNext);
                end
                Y(j, t, 1) = ytNext;
                Y(j, t, 2) = 0;
                t = t + 1;
            end
        end
        fprintf('\n');
    end
end

% emit one letter from base state
function ret = emitBaseState(X, params, E, t, yt, j)
    % fprintf('.');
    % regular E step
    if t >= params.order
        Etemp = E;
    else
        % first letters
        Etemp = matUtils.sumDim(E, 2 : 1 + params.order - t);
        Etemp = bsxfun(@times, Etemp, 1 ./ sum(Etemp, length(size(Etemp))));
    end
    % 1 x order -1
    Xlast = X(j, max(t-params.order+1,1) : t-1);
    % params.order x params.n
    subscripts = [repmat([yt ; Xlast'], [1, params.n]); 1:params.n];
    indices = matUtils.matSub2ind(size(Etemp), subscripts);
    ret = sampleMultiVarDist(Etemp(indices)');
end


% distribution - N x 1
function ret = sampleMultiVarDist(distribution)
    assert(size(distribution, 2) == 1);
    assert(abs(sum(distribution,1) - 1) < 0.01);
    ret = find(mnrnd(1, distribution));
end

