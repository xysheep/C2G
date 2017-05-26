function m = mcl(m)
% test the explanations in stijn van dongens thesis.
%
% @author gregor :: arbylon . net

if nargin < 1
    % m contains T(G3 + I) as stochastic matrix
    load -ascii m.txt
end


p = 2;
minval = 0.00001;

%e = 1.;
emax = 0.0001;
origin_m = m;
while 1
    m = origin_m;
    e = 1.;
    while e > emax
        %fprintf('iteration %i before expansion:\n', i);
        %m

        %fprintf('iteration %i after expansion/before inflation:\n', i);
        m2 = expand(m);

        %fprintf('inflation:\n')
        [m, e] = inflate(m2, p, minval);

        %fprintf('residual energy: %f\n', e);

    end % while e
    if sum(isnan(m))==0
        break
    else
        minval = minval/10;
        emax = emax/10;
    end
end
end % mcl

% expand by multiplying m * m
% this preserves column (or row) normalisation
function m2 = expand(m)
    m2 = m * m;
end

% inflate by Hadamard potentiation
% and column re-normalisation
% prune elements of m that are below minval
function [m2, energy] = inflate(m, p, minval)
    % inflation
    m2 = m .^ p;
    % pruning
    m2(m2 < minval) = 0;
    % normalisation
    dinv = diag(1./sum(m2));
    m2 = m2 * dinv;
    % calculate residual energy
    maxs = max(m2);
    sqsums = sum(m2 .^ 2);
    energy = max(maxs - sqsums);
end



