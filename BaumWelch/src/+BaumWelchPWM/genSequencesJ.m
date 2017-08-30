
function [X, Y] = genSequencesJ(theta, params)
    fprintf('Generating sequences... \n');
    params.order = length(size(theta.E)) - 1;
    Y = zeros(params.N, params.L);
    X = zeros(params.N, params.L);
    E = exp(theta.E);
    T = exp(theta.T);
    G = exp(theta.G);
    PWMs = exp(theta.PWMs);
    lengths = theta.lengths;
    startT = exp(theta.startT);
    F = exp(theta.F);
    for j = 1:params.N
        % first letter
        Etemp = matUtils.sumDim(E, 2 : params.order);
        Etemp = bsxfun(@times, Etemp, 1 ./ sum(Etemp, length(size(Etemp))));
        Y(j, 1) = smrnd(startT);
        X(j, 1) = smrnd(Etemp(Y(j, 1), :)');
        t = 2;
        while t <= params.L
            yt = Y(j, t-1);
            if rand(1) < F(yt)
                % PWM step
                motif = smrnd(G(yt, :)');
                for i=1:lengths(motif)
                    fprintf('.')
                    X(j, t) = smrnd(PWMs(motif, :, i)');
                    Y(j, t) = yt;
                    assert(X(j, t) > 0);
                    t = t + 1;
                    if t > params.L
                        % sequence too long
                        break
                    end
                end
                if t <= params.L
                    X(j, t) = emitBaseState(X, params, E, t, yt, j);
                    Y(j, t) = yt;
                    t = t + 1;
                end
            else

                ytNext = smrnd(T(yt, :)');
                X(j, t) = emitBaseState(X, params, E, t, ytNext, j);
                Y(j, t) = ytNext;
                t = t + 1;
            end
        end
        fprintf('\n');
    end
end

% emit one letter from base state
function ret = emitBaseState(X, params, E, t, yt, j)
    fprintf('o');
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

