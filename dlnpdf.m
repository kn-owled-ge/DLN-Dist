function P = dlnpdf(W,muP,sigP,muN,sigN,rPN)
% Returns the PDF of the DLN = exp(N1)-exp(N2) distribution with 
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

% Create matrices and define f_W
MU = [muP;muN];
SIG = [sigP^2 sigP*sigN*rPN; sigP*sigN*rPN sigN^2];
f_W = @(w) f_DLN(w,MU,SIG);

% If more than 300 pts, set up interpolation
if (numel(W)>300)
   % Create interpolation points
   W_1 = hermitequad(300)';
   W_1 = (W_1-W_1(1))./(W_1(end)-W_1(1));
   W_1 = sinh(W_1.*(asinh(max(W))-asinh(min(W))) + asinh(min(W)));

   % Calculate value at interpolation points
   f_1 = zeros(size(W_1));
   parfor i=1:numel(W_1)
      f_1(i) = f_W(W_1(i));
   end

   % Get interpolant and calculate
   f_W = griddedInterpolant(asinh(W_1),f_1,'linear','nearest');
   P   = max(0,f_W(asinh(W)));
else
   % Just calculate directly
   P   = zeros(size(W));
   parfor i=1:numel(W)
      P(i) = max(0,f_W(W(i)));
   end
end

P = max(1e-30,P);

% Flip vector if required
if TR==1
   P = P';
end
end



% pdf of DLN, W[1xK] MU[2,1] SIG[2,2]
function P = f_DLN(W,MU,SIG)
   warnStruct = warning('query');
   warning('off','all');
   P = max(integral( @(X_) f_MVLN([X_;X_-W],MU,SIG),max(0,W),Inf,'RelTol',1e-8,'AbsTol',1e-10),0);
   warning(warnStruct);
end



% pdf of MVLN, Y[NxK] MU[N,1] SIG[N,N]
function P = f_MVLN(Y,MU,SIG)
   N      = size(Y,1);
   Z      = sum(Y<1e-50);
   Y      = max(Y,1e-50);
   dY     = log(Y) - MU;
   Jacob  = prod(Y).^(-1);
   dG     = sum(dY.*(SIG\dY));
   P      = (2*pi)^(-N/2).*det(SIG)^(-1/2).*Jacob.*exp(-(1/2).*dG);
   P(Z>0) = 0;
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

