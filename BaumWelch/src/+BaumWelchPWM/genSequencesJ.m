
function [X, Y] = genSequencesJ(theta, params)
    fprintf('Generating sequences... \n');
    params.order = length(size(theta.E)) - 1;
    Y = zeros(params.N, params.L);
    X = zeros(params.N, params.L);
    for j = 1:params.N
        Ep = matUtils.sumDim(theta.E, 2 : params.order);
        Y(j, 1) = smrnd(theta.startT);
        X(j, 1) = smrnd(Ep(Y(j, 1), :)');
        t = 2;
        while t <= params.L
            tMode = smrnd(theta.T(Y(j, t-1), :)');
            Y(j, t) = tMode;
            if rand(1) < theta.F(tMode)
                % PWM step
                motif = smrnd(theta.M(tMode, :)');
                for i=1:theta.lengths(motif)
                    X(j, t) = smrnd(theta.PWMs(motif, :, i)');
                    Y(j, t) = tMode;
                    assert(X(j, t) > 0)
                    t = t + 1;
                    if t > params.L
                        break
                    end
                end
            else
                % regular theta.E step
                if t >= params.order
                    Etemp = theta.E;
                else
                    Etemp = matUtils.sumDim(theta.E, 2 : 1 + params.order - t);
                end
                Xlast = X(j, max(t-params.order+1,1) : t-1);
                % params.order x params.n*params.N
                subscripts = [repmat([tMode ; Xlast'], [1, params.n]); 1:params.n];
                indices = matUtils.matSub2ind(size(Etemp), subscripts);
                X(j, t) = smrnd(Etemp(indices)');
                assert(X(j, t) > 0)
                t = t + 1;
            end
        end
    end
end



% P - params.N x 1
function ret = smrnd(P)
    assert(size(P, 2) == 1);
    P = P ./ sum(P, 1);
    ret = find(mnrnd(1, P));
    assert(ret > 0);
end

