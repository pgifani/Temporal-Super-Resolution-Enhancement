function [weights,used,noise2,errbars,basis] =SBL(PHI,t,noise2,termination,adaptive,optimal,scale)
%------------------------------------------------------------------
% Sparse Bayesian Learning
% Input for BCS:
%   PHI: projection matrix
%   t:   CS measurements
%   noise2: initial noise variance
%      If measurement noise exists and/or w is not truely sparse, 
%             then noise2 = std(t)^2/1e2 (suggested)
%      If no measurement noise and w is truely sparse,
%             then noise2 = std(t)^2/1e6 (suggested)
%      This term is in fact not updated in the implementation to allow 
%      the fast algorithm. For this reason, you are recommended to use
%      mt_CS.m, in which the noise variance is marginalized.
%   termination: threshold for stopping the algorithm (suggested value: 1e-8)
% Input for Adaptive CS:
%   adaptive: generate basis for adpative CS? (default: 0)
%   optimal: use the rigorous implementation of adaptive CS? (default: 1)
%   scale: diagonal loading parameter (default: 0.1)
% Output:
%   weights:  sparse weights
%   used:     the positions of sparse weights
%   noise2:   re-estimated noise variance
%   errbars:  one standard deviation around the sparse weights
%   basis:    if adaptive==1, then basis = the next projection vector
%
if nargin < 5
    adaptive = 0;
end
if nargin < 6
    optimal = 1;
end
if nargin < 7
    scale = 0.1;
end

% find initial alpha
[N,M] = size(PHI);
PHIt = PHI'*t;
PHI2 = sum(PHI.^2)';
ratio = (PHIt.^2)./PHI2;
[maxr,index] = max(ratio);
alpha = PHI2(index)/(maxr-noise2);
% compute initial miu, Sigma, S, Q
phi = PHI(:,index);
Hessian = alpha + phi'*phi/noise2;
Sigma = 1/Hessian;
miu = Sigma*PHIt(index)/noise2;
woodbury= PHI'*phi/noise2;
S = PHI2/noise2-Sigma*woodbury.^2;
Q = PHIt/noise2-Sigma*PHIt(index)/noise2*woodbury;
%
for count = 1:20000

    s = S; q = Q;
    s(index) = alpha.*S(index)./(alpha-S(index));
    q(index) = alpha.*Q(index)./(alpha-S(index));
    theta = q.^2-s;

    % choice the next alpha that maximizes marginal likelihood
    ml = -inf*ones(1,M);
    ig0 = find(theta>0);
    % index for re-estimate
    [ire_estimate,~,which] = intersect(ig0,index);
    if ~isempty(ire_estimate)
        Alpha = s(ire_estimate).^2./theta(ire_estimate);
        delta = (alpha(which)-Alpha)./(Alpha.*alpha(which));
        ml(ire_estimate) = Q(ire_estimate).^2.*delta./(S(ire_estimate).*delta+1)-log(1+S(ire_estimate).*delta);
    end
    % index for adding
    iadd = setdiff(ig0,ire_estimate);
    if ~isempty(iadd)
        ml(iadd) = (Q(iadd).^2-S(iadd))./S(iadd)+log(S(iadd)./(Q(iadd).^2));
    end
    is0 = setdiff([1:M],ig0);
    % index for deleting
    [idelete,~,which] = intersect(is0,index);
    if ~isempty(idelete)
        ml(idelete) = Q(idelete).^2./(S(idelete)-alpha(which))-log(1-S(idelete)./alpha(which));
    end

    [ML(count),idx] = max(ml);
    % check if terminates?
    if count > 2 & abs(ML(count)-ML(count-1)) < abs(ML(count)-ML(1))*termination
        break;
    end

    % update alphas
    which = find(index==idx);
    if theta(idx) > 0
        if ~isempty(which) % re-estimate
            Alpha = s(idx)^2/theta(idx);
            Sigmaii = Sigma(which,which); miui = miu(which); Sigmai = Sigma(:,which);
            delta = Alpha-alpha(which);
            ki = delta/(1+Sigmaii*delta);
            miu = miu-ki*miui*Sigmai;
            Sigma = Sigma-ki*Sigmai*Sigmai';
            comm = PHI'*(phi*Sigmai)/noise2;
            S = S + ki*comm.^2;
            Q = Q + ki*miui*comm;
            %
            alpha(which) = Alpha;
        else % adding
            Alpha = s(idx)^2/theta(idx);
            phii = PHI(:,idx); Sigmaii = 1/(Alpha+S(idx)); miui = Sigmaii*Q(idx);
            comm1 = Sigma*(phi'*phii)/noise2;
            ei = phii-phi*comm1;
            off = -Sigmaii*comm1;
            Sigma = [Sigma+Sigmaii*comm1*comm1', off; off', Sigmaii];
            miu = [miu-miui*comm1; miui];
            comm2 = PHI'*ei/noise2;
            S = S - Sigmaii*comm2.^2;
            Q = Q - miui*comm2;
            %
            index = [index;idx];
            alpha = [alpha;Alpha];
            phi = [phi,phii];
        end
    else
        if ~isempty(which) % deleting
            Sigmaii = Sigma(which,which); miui = miu(which); Sigmai = Sigma(:,which);
            Sigma = Sigma-Sigmai*Sigmai'/Sigmaii; Sigma(:,which) = []; Sigma(which,:) = [];
            miu  = miu-miui/Sigmaii*Sigmai; miu(which) = [];
            comm = PHI'*(phi*Sigmai)/noise2;
            S = S + comm.^2/Sigmaii;
            Q = Q + miui/Sigmaii*comm;
            %
            index(which) = [];
            alpha(which) = [];
            phi(:,which) = [];
        end
    end

end
weights	= miu;
used = index;
% re-estimated noise2
noise2 = sum((t-phi*miu).^2)/(N-length(index)+alpha'*diag(Sigma)); 
errbars = sqrt(diag(Sigma));

% generate a basis for adaptive CS?
if adaptive
    if optimal
        [V,D] = eig(Sigma);
        [~,idx] = max(diag(D));
        basis = V(:,idx)';
    else
        temp = phi'*phi/noise2;
        Sigma_inv = temp + scale*mean(diag(temp))*eye(length(used));
        [V,D] = eig(Sigma_inv);
        [~,idx] = min(diag(D));
        basis = V(:,idx)';
    end
end
