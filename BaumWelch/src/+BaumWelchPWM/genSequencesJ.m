
function [X, Y] = genSequencesJ(startT, T, E, M, F, L, N, PWMs, lengths)
    fprintf('Generating sequences... \n');
    order = length(size(E)) - 1;
    n = size(PWMs, 2);
    Y = zeros(N, L);
    X = zeros(N, L);
    for j = 1:N
        Ep = matUtils.sumDim(E, 2 : order);
        Y(j, 1) = smrnd(startT);
        X(j, 1) = smrnd(Ep(Y(j, 1), :)');
        t = 2;
        while t <= L
            tMode = smrnd(T(Y(j, t-1), :)');
            Y(j, t) = tMode;
            if rand(1) < F(tMode)
                % PWM step
                motif = smrnd(M(tMode, :)');
                for i=1:lengths(motif)
                    X(j, t) = smrnd(PWMs(motif, :, i)');
                    Y(j, t) = tMode;
                    assert(X(j, t) > 0)
                    t = t + 1;
                    if t > L
                        break
                    end
                end
            else
                % regular E step
                if t >= order
                    Etemp = E;
                else
                    Etemp = matUtils.sumDim(E, [2 : 1 + order - t]);
                end
                Xlast = X(j, max(t-order+1,1) : t-1);
                % order x n*N
                subscripts = [repmat([tMode ; Xlast'], [1, n]); 1:n];
                indices = matUtils.matSub2ind(size(Etemp), subscripts);
                X(j, t) = smrnd(Etemp(indices)');
                assert(X(j, t) > 0)
                t = t + 1;
            end
        end
    end
end



% P - N x 1
function ret = smrnd(P)
    assert(size(P, 2) == 1);
    P = P ./ sum(P, 1);
    ret = find(mnrnd(1, P));
    assert(ret > 0);
end
% populationM - a x b, where each column is a distribution of getting
function ret = mmnrand(populationM)
    sum(repmat(rand(1, N), [n,1]) > cumsum(Ep'), 1) + 1;
end

