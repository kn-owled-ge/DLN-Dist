function [m1,m2,m3,m4,m5] = dlnmom(muP,sigP,muN,sigN,rPN)
% Return the first five central moments of the DLN.
% muP,sigP,muN,sigN,rPN must be scalars, with ranges:
% -inf < muP,muN < inf ; 0 <= sigP,sigN < inf ; abs(rPN) <= 1
%
% Alternatively, the first argument (muP) can be a 1X5 vector containing
% the 5 paramteres.
%
% Note that returned kurtosis is regular rather than excess kurtosis.
%
% For theoretical derivation, see Parham (2022)
%
% See also: DLNPDF, DLNCDF, DLNINV, DLNFIT, DLNPAR, DLNRND

% Check parameters
if nargin==1
   if numel(muP)~=5
      error('stats:DLN:BadInputs','First input must be a 5-vector.');
   end
   rPN  = muP(5);
   sigN = muP(4);
   muN  = muP(3);
   sigP = muP(2);
   muP  = muP(1);
elseif nargin~=5
    error('stats:DLN:BadInputs','Requires five input arguments.');
end
if ~isscalar(muP) || ~isscalar(muN) || ~isscalar(sigP) || ~isscalar(sigN) || ~isscalar(rPN) 
   m1 = NaN;
   return;
elseif sigP < 0 || sigN < 0 || abs(rPN) > 1
   m1 = NaN;
   return;
end

% MU, SIG of the BVN and MVLN dists
MU_BVN = [muP,muN]';
SIG_BVN = [sigP^2 sigP*sigN*rPN; sigP*sigN*rPN sigN^2];
MU_MVLN = exp(MU_BVN+0.5*diag(SIG_BVN));
SIG_MVLN = exp(log(MU_MVLN*[1,1])+log((MU_MVLN*[1,1])')).*(exp(SIG_BVN)-1);

% The first two moments are easy.
m1 = MU_MVLN(1)-MU_MVLN(2);
m2 = SIG_MVLN(1,1)+SIG_MVLN(2,2)-2*SIG_MVLN(1,2);

% Pascal is the first few lines of the negative Pascal triangle
Pascal = [+1 0 0 0 0 0 ; +1 -1 0 0 0 0 ; +1 -2 +1 0 0 0 ; +1 -3 +3 -1 0 0 ; +1 -4 +6 -4 +1 0 ; +1 -5 +10 -10 +5 -1];

% E is a helper matrix with E[i,j] = Ex[Y_i^(j-1)];
E = zeros(2,6);
for j=1:6
   E(1,j) = exp((j-1)*MU_BVN(1) + 0.5*(j-1)^2*SIG_BVN(1,1));
   E(2,j) = exp((j-1)*MU_BVN(2) + 0.5*(j-1)^2*SIG_BVN(2,2));
end

% Get the S function-based skewness
S  = 0;
for s=1:4
   S1 = zeros(4,4);
   t = (4-s+1);
   for u=1:t
      for v=1:s
         S1(u,v) =         Pascal(t,u)*E(1,t-u+1)*MU_MVLN(1)^(u-1);
         S1(u,v) = S1(u,v)*Pascal(s,v)*E(2,s-v+1)*MU_MVLN(2)^(v-1);
         S1(u,v) = S1(u,v)*exp((t-u)*(s-v)*SIG_BVN(1,2));
      end
   end
   S = S + Pascal(4,s)*sum(sum(S1));
end
m3 = S/m2^(3/2);

% Get the K function-based kurtosis
K  = 0;
for k=1:5
   K1 = zeros(5,5);
   t = (5-k+1);
   for u=1:t
      for v=1:k
         K1(u,v) =         Pascal(t,u)*E(1,t-u+1)*MU_MVLN(1)^(u-1);
         K1(u,v) = K1(u,v)*Pascal(k,v)*E(2,k-v+1)*MU_MVLN(2)^(v-1);
         K1(u,v) = K1(u,v)*exp((t-u)*(k-v)*SIG_BVN(1,2));
      end
   end
   K = K + Pascal(5,k)*sum(sum(K1));
end
m4 = K/m2^(4/2);

% Get the A function-based tail-assymmetry
A  = 0;
for a=1:6
   A1 = zeros(6,6);
   t = (6-a+1);
   for u=1:t
      for v=1:a
         A1(u,v) =         Pascal(t,u)*E(1,t-u+1)*MU_MVLN(1)^(u-1);
         A1(u,v) = A1(u,v)*Pascal(a,v)*E(2,a-v+1)*MU_MVLN(2)^(v-1);
         A1(u,v) = A1(u,v)*exp((t-u)*(a-v)*SIG_BVN(1,2));
      end
   end
   A = A + Pascal(6,a)*sum(sum(A1));
end
m5 = A/m2^(5/2);

end
