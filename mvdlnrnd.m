function W = mvdlnrnd(N,K,MU,SIG,muDLN,sigDLN,rhoDLN)
% MVDLNRND multi-variate DLN random number generator.
% W = MVDLNRND(N,K,MU,SIG,muDLN,sigDLN,rhoDLN) draws K obs W[N,K] from the 
% N-dimensional location-scale DLN with location MU, scale SIG, and 
% generator SymDLN(muDLN,sigDLN,rhoDLN).
%
% muDLN,sigDLN,rhoDLN must be scalars, with the ranges: 
% -inf < muDLN < inf ; 0 <= sigDLN < inf ; abs(rhoDLN) <= 1
%
% For theoretical derivation, see Parham (2022)

% Check parameters
if nargin==1
   K      = 1;
   MU     = zeros(N,1);
   SIG    = eye(N);
   muDLN  = 0;
   sigDLN = 1;
   rhoDLN = 0;
elseif nargin==2
   MU     = zeros(N,1);
   SIG    = eye(N);
   muDLN  = 0;
   sigDLN = 1;
   rhoDLN = 0;
elseif nargin==4
   muDLN  = 0;
   sigDLN = 1;
   rhoDLN = 0;
elseif nargin~=7
    error('stats:DLN:BadInputs','Requires 1, 2, 4, or 7 input arguments.');
end
if (N<1)
   error('stats:DLN:BadInputs','Expected N>=1');
end
if (any(size(MU)~=[N 1]) || any(size(SIG)~=[N N]))
   error('stats:DLN:BadInputs','Expected size(MU)=[N 1] and size(SIG)=[N,N]');
end
if ~isscalar(muDLN) || ~isscalar(sigDLN) || ~isscalar(rhoDLN) 
   error('stats:DLN:BadInputs','Expected size(muDLN,sigDLN,rhoDLN)=[1,1]');
end
if sigDLN <= 0 || abs(rhoDLN) >= 1
   error('stats:DLN:BadInputs','Expected sigDLN>0, abs(rDLN)<1');
end

% Keep persistence
persistent Ist
if (isempty(Ist) || any(abs([Ist.muDLN Ist.sigDLN Ist.rhoDLN Ist.N]-[muDLN sigDLN rhoDLN N])>1E-3))
   % Save current state
   Ist.muDLN  = muDLN;
   Ist.sigDLN = sigDLN;
   Ist.rhoDLN = rhoDLN;
   Ist.N      = N;

   % Create normalization integral
   r_I = sinh(1:1:30);
   v_I = r_I.^(N-1).*dlnpdf(r_I,muDLN,sigDLN,muDLN,sigDLN,rhoDLN);
   i_I = find(v_I==0,1,'first');
   if isempty(i_I)
      i_I=numel(v_I);
   end
   Int = integral(@(r_) r_.^(N-1).*dlnpdf(r_,muDLN,sigDLN,muDLN,sigDLN,rhoDLN),0,r_I(i_I),'RelTol',1e-8,'AbsTol',1e-10);

   % Get F_R(r)
   r_F = sinh([0:0.01:1 , 1.05:0.05:asinh(r_I(i_I))]);
   c_F = arrayfun(@(r_H) integral(@(r_) dlnpdf(r_,muDLN,sigDLN,muDLN,sigDLN,rhoDLN).*r_.^(N-1)./Int,0,r_H,'RelTol',1e-8,'AbsTol',1e-10),r_F);
   i_F = min(find(c_F>1-1E-10,1,'first'),find(c_F(2:end)<c_F(1:end-1),1,'first'));
   r_F = r_F(1:i_F);
   c_F = c_F(1:i_F);

   % Get invF_R(r) interpolant
   Ist.invFR = griddedInterpolant(c_F,asinh(r_F),'makima','nearest');
end

% Generate standardized Z variates
H = random('Uniform',0,1,[1,K]);
R = sinh(Ist.invFR(H));
U = mvnrnd(zeros(N,1),eye(N),K)';
U = (U./vecnorm(U));
Z = R.*U;

% Make location-scale
S05 = Riccati(SIG);
W   = MU + S05*Z;

end



% Riccati square-root of Sig
function Sig05 = Riccati(Sig)
   % Sig05 = Sig05' ; Sig05^2=Sig 
   [V,D] = eig(Sig);
   Sig05 = V*D.^0.5*V;
end
