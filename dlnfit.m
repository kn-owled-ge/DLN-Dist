function par = dlnfit(W)
% Returns estimates of the five parameters in a fit of the data W to the
% DLN distribution.
%
% PARAM is a 1X5 vector containing the muP, sigP, muN, sigN, rPN params.
%
% For theoretical derivation, see Parham (2022)


% Check parameters
if nargin==1
   % Require sufficient obs
   if numel(W(W<0))<1E2 || numel(W(W>0))<1E2
      error('stats:DLN:BadInputs','Requires at least 100 pos and neg input obs each.');
   end
else
   error('stats:DLN:BadInputs','Requires exactly one input argument.');
end

% Prepare optimization 
mom      = [median(log(W(W>0))),iqr(log(W(W>0)))/1.35,median(log(-W(W<0))),iqr(log(-W(W<0)))/1.35];
mom      = [ones(5,1)*mom , [-0.8 -0.3 0 0.3 0.8]'];
lbP      = [-9 0 -9 0 -1];
ubP      = [ 9 9  9 9  1];
dlnObj   = @(p) -sum(log(max(1e-30,dlnpdf(W,p))))/numel(W);
options  = optimoptions('fmincon');
options  = optimoptions(options, 'Display'       , 'none'           );  % 'iter'/'final'/'none'
options  = optimoptions(options, 'MaxFunEvals'   , 1000             );
options  = optimoptions(options, 'Algorithm'     , 'interior-point' );  % 'interior-point'/'active-set'/'sqp'
options  = optimoptions(options, 'TolX'          , 1E-3             );
options  = optimoptions(options, 'UseParallel'   , false				  );
options  = optimoptions(options, 'StepTolerance' , 1E-3             );
options  = optimoptions(options, 'DiffMinChange' , 1E-3             );
options  = optimoptions(options, 'TypicalX'      , mom(2,:)         );
ms       = MultiStart('UseParallel',true,'Display','none');
problem  = createOptimProblem('fmincon','objective',dlnObj,'x0',mom(2,:),'lb',lbP,'ub',ubP,'options',options);
points   = CustomStartPointSet(mom);

% Run MultiStart Optimization
[~,~,~,~,solutions] = run(ms,problem,points);
if numel(solutions) == 0
   par = NaN;
else
   par = solutions(1).X;
end
end