function T = transitionRate(DT, X)
% FIXME Generator function

    % Markov states are triangle midpoints
    numStates = size(DT.ConnectivityList, 1);
    dim = size(DT.ConnectivityList, 2);
    T = zeros(numStates, numStates);

    K = DT.neighbors;

    % Pseudocounts (Eq. 11)
    % Uniform prior
    binCount = ones(numStates, 1);
    for ii = 1:numStates
        for jj = 1:dim
            n = K(ii, jj);
            if ~isnan(n)
                T(n, ii) = 1;
            end
        end
    end

    for ii = 1:size(X, 3)
        x = X(:,:,ii);
        for jj = 1:size(x, 2)-1

            startBin = pointLocation(DT, x(:, jj)');
            endBin   = pointLocation(DT, x(:, jj+1)');
            if isnan(startBin) || isnan(endBin)
                continue
            end

            neighborhood = K(startBin, :);

            if any(neighborhood == endBin)
                T(endBin, startBin) = T(endBin, startBin) + 1;
                binCount(startBin) = binCount(startBin) + 1;
            elseif startBin == endBin
                binCount(startBin) = binCount(startBin) + 1;
            end
        end
    end

    % Normalize columns
    T = T ./ binCount';

    % Enforce column sum to zero (generator matrix)
    T = T - diag(sum(T, 1));
end