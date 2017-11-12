
function [X, Y] = genSequencesJ(theta, params)
    fprintf('Generating sequences... \n');
    params.order = length(size(theta.E)) - 1;
    Y = zeros(params.N, params.L, 2);
    X = zeros(params.N, params.L);
    E = exp(theta.E);
    T = exp(theta.T);
    G = exp(theta.G);
    PWMs = params.PWMs;
    lengths = params.lengths;
    startT = exp(theta.startT);
    for j = 1:params.N
        % first letter
        Etemp = matUtils.sumDim(E, 2 : params.order);
        Etemp = bsxfun(@times, Etemp, 1 ./ sum(Etemp, length(size(Etemp))));
        Y(j, 1, 1) = smrnd(startT);
        Y(j, 1, 2) = 0;
        X(j, 1) = smrnd(Etemp(Y(j, 1), :)');
        fprintf('%d ', Y(j, 1, 1));
        t = 2;
        while t <= params.L
            yt = Y(j, t-1, 1);
            state = smrnd([T(yt, :), G(yt, :)]');
            if params.m < state
                % PWM step
                motif = state - params.m;
                for i=1:lengths(motif)
                    X(j, t) = smrnd(PWMs(motif, :, i)');
                    fprintf('%d',motif)
                    % fprintf('%d',X(j, t))
                    Y(j, t, 1) = yt;
                    Y(j, t, 2) = motif;
                    assert(X(j, t) > 0);
                    t = t + 1;
                    if t > params.L
                        % sequence too long
                        break
                    end
                end
                if t <= params.L
                    X(j, t) = emitBaseState(X, params, E, t, yt, j);
                    Y(j, t, 1) = yt;
                    Y(j, t, 2) = 0;
                    t = t + 1;
                end
            else
                ytNext = state;
                X(j, t) = emitBaseState(X, params, E, t, ytNext, j);
                % if ytNext ~= yt
                %     fprintf('!')
                % end
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
    ret = smrnd(Etemp(indices)');
end


% P - N x 1
function ret = smrnd(P)
    assert(size(P, 2) == 1);
    assert(abs(sum(P,1) - 1) < 0.01);
    ret = find(mnrnd(1, P));
end

