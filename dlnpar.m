function Opar = dlnpar(Ipar,inv)
% Opar = DLNPAR(Ipar,inv) reparametrizes the DLN parameters in the 5-vector 
% Ipar.
%
% If inv== 1, it translates from M-S to A-E representation.
% If inv==-1, it translates from A-E to M-S representation.
% If inv== 0, it does nothing.
%
% For theoretical derivation, see Parham (2022)
%
% See also: DLNPDF, DLNCDF, DLNINV, DLNFIT, DLNMOM, DLNRND

% Check parameters
if nargin~=2
   error('stats:DLN:BadInputs','Requires two input arguments.');
elseif numel(Ipar)~=5
   error('stats:DLN:BadInputs','First input must be a 5-vector.');
end

% Which case?
switch inv
   case 0
      % Just return Ipar as identity
      Opar = Ipar;
      return;
   case 1
      % Ipar is M-S, return A-E
      Opar    = zeros(size(Ipar));
      Opar(1) = exp(Ipar(1)+0.5*Ipar(2)^2) - exp(Ipar(3)+0.5*Ipar(4)^2);
      Opar(2) = exp(Ipar(1)+0.5*Ipar(2)^2) + exp(Ipar(3)+0.5*Ipar(4)^2);
      Opar(3) = (exp(Ipar(2)^2)-1) - (exp(Ipar(4)^2)-1);
      Opar(4) = (exp(Ipar(2)^2)-1) + (exp(Ipar(4)^2)-1);
      Opar(5) = exp(Ipar(5)*Ipar(2)*Ipar(4))-1;
      Opar    = asinh(Opar);
      return;
   case -1
      % Ipar is A-E - return M-S
      Opar    = zeros(size(Ipar));
      Ipar    = sinh(Ipar);
      Opar(2) = sqrt(log((Ipar(4) + Ipar(3))/2 + 1));
      Opar(4) = sqrt(log((Ipar(4) - Ipar(3))/2 + 1));
      Opar(1) = log((Ipar(2) + Ipar(1))/2) - 0.5*Opar(2)^2;
      Opar(3) = log((Ipar(2) - Ipar(1))/2) - 0.5*Opar(4)^2;
      Opar(5) = log(1+Ipar(5))/(Opar(2)*Opar(4));
      return;
   otherwise
      error('stats:DLN:BadInputs','Unknown representation.');
end
end