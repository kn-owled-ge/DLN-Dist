function C = mvdlncdf(W,MU,SIG,muDLN,sigDLN,rhoDLN)
% P = MVDLNCDF(W,MU,SIG,muDLN,sigDLN,rhoDLN) returns the cdf of the 
% multi-variate DLN with location MU and scale SIG, generated by the 
% SymDLN(muDLN,sigDLN,rhoDLN). W[NxK] is N-dimensional DLN with K
% observations.
%
% For theoretical derivation, see Parham (2022)

% Check parameters
N = size(W,1);
if nargin==1
   MU     = zeros(N,1);
   SIG    = eye(N);
   muDLN  = 0;
   sigDLN = 1;
   rhoDLN = 0;
elseif nargin==3
   muDLN  = 0;
   sigDLN = 1;
   rhoDLN = 0;
elseif nargin~=6
    error('stats:DLN:BadInputs','Requires 1, 3, or 6 input arguments.');
end
if (N~=2)
   error('stats:DLN:BadInputs','MVDLNCDF is currently implemented only for N=2');
end
if (any(size(MU)~=[N 1]) || any(size(SIG)~=[N N]))
   error('stats:DLN:BadInputs','Expected size(MU)=[N,1] and size(SIG)=[N,N]');
end
if ~isscalar(muDLN) || ~isscalar(sigDLN) || ~isscalar(rhoDLN) 
   error('stats:DLN:BadInputs','Expected size(muDLN,sigDLN,rhoDLN)=[1,1]');
end
if sigDLN <= 0 || abs(rhoDLN) >= 1
   error('stats:DLN:BadInputs','Expected sigDLN>0, abs(rDLN)<1');
end

% Keep persistence
persistent Ist
if (isempty(Ist) || any(abs([Ist.muDLN Ist.sigDLN Ist.rhoDLN Ist.N]-[muDLN sigDLN rhoDLN N])>1E-3)) || any(abs(Ist.MU-MU)>1e-3) || any(any(abs(Ist.SIG-SIG)>1e-3))
   % Save current state
   Ist.muDLN  = muDLN;
   Ist.sigDLN = sigDLN;
   Ist.rhoDLN = rhoDLN;
   Ist.N      = N;
   Ist.MU     = MU;
   Ist.SIG    = SIG;

   % Create interpolation points
   S05 = Riccati(SIG);
   M_1 = muDLN+15*sigDLN;
   W_1 = sinh(-M_1:1:M_1);
   W_2 = [W_1 ; W_1];
   W_2 = MU + S05*W_2;
   f_2 = mvdlnpdf(W_2,MU,SIG,muDLN,sigDLN,rhoDLN);
   i_L = find(f_2>0,1,'first');
   if isempty(i_L)
      i_L=1;
   end
   i_H = find(f_2>0,1,'last');
   if isempty(i_H)
      i_H=numel(W_2);
   end
   M_L = asinh(W_1(i_L));
   M_H = asinh(W_1(i_H));
   W_3 = hermitequad(300)';
   W_3 = (W_3-W_3(1))./(W_3(end)-W_3(1));
   W_3 = sinh(W_3.*(M_H-M_L) + M_L);
   W_4 = [W_3 ; W_3];
   W_4 = MU + S05*W_4;

   % Mesh and calculate CDF
   [X_4,Y_4] = meshgrid(W_4(1,:),W_4(2,:));
   vX_4      = reshape(X_4,1,[]);
   vY_4      = reshape(Y_4,1,[]);
   vf_4      = mvdlnpdf([vX_4;vY_4],MU,SIG,muDLN,sigDLN,rhoDLN);
   f_4       = reshape(vf_4,size(X_4));
   F_4       = cumtrapz(W_4(1,:),cumtrapz(W_4(2,:),f_4,2));

   % Make sure extrapolation yields 0 and 1 at edges
   F_4 = F_4 - F_4(1,1);
   F_4 = F_4./F_4(end,end);

   % Get interpolant
   Ist.FW = griddedInterpolant({asinh(W_4(1,:));asinh(W_4(2,:))},F_4,'makima','nearest');
end

% Calculate
C     = Ist.FW(asinh(W)')';
end



% Riccati square-root of Sig
function Sig05 = Riccati(Sig)
   % Sig05 = Sig05' ; Sig05^2=Sig 
   [V,D] = eig(Sig);
   Sig05 = V*D.^0.5*V;
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
