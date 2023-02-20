function C = dlncdf(W,muP,sigP,muN,sigN,rPN)
% Returns the CDF of the DLN = exp(N1)-exp(N2) distribution with 
% [N1;N2]~MVN([muP;muN],[sigP^2 sigP*sigN*rPN; sigP*sigN*rPN sigN^2]),
% at the values of W[1xK]. muP,sigP,muN,sigN,rPN must be scalars, with the 
% ranges: -inf < muP,muN < inf ; 0 <= sigP,sigN < inf ; abs(rPN) <= 1
%
% Alternatively, the second argument (muP) can be a 1X5 vector containing
% the 5 paramteres.
%
% For theoretical derivation, see Parham (2022)

% Check parameters
if nargin==2
   if numel(muP)~=5
      error('stats:DLN:BadInputs','Second input must be a 5-vector.');
   end
   rPN  = muP(5);
   sigN = muP(4);
   muN  = muP(3);
   sigP = muP(2);
   muP  = muP(1);
elseif nargin~=6
    error('stats:DLN:BadInputs','Requires six input arguments.');
end
if ~isscalar(muP) || ~isscalar(muN) || ~isscalar(sigP) || ~isscalar(sigN) || ~isscalar(rPN) 
   error('stats:DLN:BadInputs','Bad arguments provided.');
end
if sigP < 0 || sigN < 0 || abs(rPN) > 1
   error('stats:DLN:BadInputs','Bad arguments provided.');
end

% Make sure orientation is predictable
TR = 0;
if size(W,1)~=1
   TR = 1;
   W  = W';
end

% Create interpolation points
M_L = min([W(~isinf(W)),-exp(muN+10*sigN)]);
M_H = max([W(~isinf(W)),+exp(muP+10*sigP)]);
W_1 = sinh(linspace(asinh(M_L),asinh(M_H),20));
f_1 = dlnpdf(W_1,muP,sigP,muN,sigN,rPN);
i_L = find(f_1>1e-15,1,'first');
if isempty(i_L)
   i_L=1;
end
i_H = find(f_1>1e-15,1,'last');
if isempty(i_H)
   i_H=numel(W_1);
end
i_L = max(1,i_L-1);
i_H = min(numel(W_1),i_H+1);
M_L = asinh(W_1(i_L));
M_H = asinh(W_1(i_H));
W_2 = hermitequad(300)';
W_2 = (W_2-W_2(1))./(W_2(end)-W_2(1));
W_2 = sinh(W_2.*(M_H-M_L) + M_L);
f_2 = dlnpdf(W_2,muP,sigP,muN,sigN,rPN);
F_2 = cumtrapz(W_2,f_2);
F_2 = F_2./max(F_2);

% Get interpolant
F_W = griddedInterpolant(asinh(W_2),F_2,'linear','nearest');

% Do calculation
C = F_W(asinh(W));

% Flip vector if required
if TR==1
   C = C';
end
end



% Generate nodes and weights for Hermite-Gauss quadrature.
% Note that x is a column vector and w is a row vector.
function [x, w] = hermitequad(n)
   u      = sqrt([1 : n-1]/2);
   [V, L] = eig(diag(u,1) + diag(u,-1));
   [x, i] = sort(diag(L));
   Vtop   = V(1, :);
   Vtop   = Vtop(i);
   w      = sqrt(pi)*Vtop.^2;
end
