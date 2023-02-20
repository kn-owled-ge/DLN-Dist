function W = dlnrnd(muP,sigP,muN,sigN,rPN,sizeOut)
% Draws random variates from the DLN = exp(N1)-exp(N2) distribution with 
% [N1;N2]~MVN([muP;muN],[sigP^2 sigP*sigN*rPN; sigP*sigN*rPN sigN^2]).
% muP,sigP,muN,sigN,rPN must be scalars, with the ranges: 
% -inf < muP,muN < inf ; 0 <= sigP,sigN < inf ; abs(rPN) <= 1
%
% For theoretical derivation, see Parham (2022)

% Check parameters
if nargin < 5
    error('stats:DLN:BadInputs','Requires at least five input arguments.'); 
elseif ~isscalar(muP) || ~isscalar(muN) || ~isscalar(sigP) || ~isscalar(sigN) || ~isscalar(rPN) 
   error('stats:DLN:BadInputs','Parameter arguments must be scalars');
elseif sigP < 0 || sigN < 0 || abs(rPN) > 1
   error('stats:DLN:BadInputs','Required: sigP>=0, sigN>=0, abs(rPN)<=1');
end
if nargin==5
   sizeOut = [1,1];
end

% Generate sample
mu  = [muP;muN];
sig = [sigP^2 sigP*sigN*rPN; sigP*sigN*rPN sigN^2];
BVN = mvnrnd(mu,sig,prod(sizeOut));
DLN = exp(BVN(:,1))-exp(BVN(:,2));
W   = squeeze(reshape(DLN,sizeOut));
end
