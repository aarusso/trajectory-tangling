%% TANGLEANALYSIS.m
% This analysis determines how 'tangled' a given dataset is. That is, how
% often the projections in high-D space cross over each other in different
% directions. The tangle index (q) of a pair of time points is given by the
% squared  difference in their derivative divided by the squared difference
% in their state. A constant (proportional to the variance of the data) is
% added to the denominator to prevent hyperbolic growth.
%
%
% Note that although tangling can be calculated for every pair of times,
% the value that is most informative is the maximum a given time point
% becomes tangled with any other time point. This value is the output 'Q'.
%
% SYNTAX
%   [ Q, out] = tangleAnalysis( Data, .001 ); % for data collected at 1kHz
%
% INPUTS:
%
% Required Inputs:
% * Data (struct) - a C-dimensional structure where C is the number of
% conditions. For a given condition, Data(c).A should hold the data (e.g.
% firing rates). Each column of A corresponds to a neuron, and each row to
% a timepoint. Data(c).times and Data(c).analyzeTimes are optional
% fields. size(Data(c).times,1) should equal size(Data(c).A,1). If
% Data(c).analyzeTimes is provided, only those entries that match between
% 'times' and 'analyzeTimes' will be used for the analysis. If analyzeTimes
% == [], all times will be used.
%
% * samplingRate (scalar) - a scalar indicating how much time in seconds
% passes with each sample in the data.
%
% Optional Inputs:
% Several parameters can be specified from the command line (see Usage
% Examples below)
% * numPCs (scalar, default: 8) - specifies how many PCs to use when
% calculating tangling
% * fracVar (scalar, default: .85) - if numPCs is empty, the number of PCs
% will be chosed based on the number needed to caputure as much variance as
% specified by fracVar.
% * cond2use (vector, default: 1:length(Data)) - specifies which conditions
% to use in calculating tangling
% * withinConditionsOnly (logical, default: false) - if desired, tangling
% will only be calculated for pairs of times that correspond to the same
% condition
% * timeStep (scalar, default: 20) - Because tangling is a pairwise
% computation, things can get quite slow if enough data is included. After
% calculating the derivitive, it can be efficient to calculate tangling
% only on every ith time sample where i is specified by timeStep. It is
% generally advisible to keep this number as low as possible without
% computation times getting unreasonably long (e.g. > ~30 seconds). Note
% that if this number is too high, aliasing could occur and tangling could
% appear lower than it actually is.
% * softenNorm (scalar, default: 5) - In the usual Churchlandian fashion,
% soft normalization is performed on the neural data such that each neuron
% is devided by its range (across all times and conditions) + some constant
% indicated by softenNorm. EMG data should be fully normalized (i.e.
% softenNorm set to 0).
% 
% OUTPUTS:
% * Q (vector) - the tangling index for each time point,  a T/timesStep x 1 vector
% * out (struct) - out is a struct containing supplementary information
% that may be useful for troubleshooting and/or analyzing the results more
% deeply. out contains the following fields:
%     - out.X: out.X is a matrix containing the principal componenets that
%     were analyzed after pre-processing. X is T x K where T is the
%     number of time points to be analyzes (corresponing to
%     Data.analyzeTimes or Data.times if .analyze times is not specified),
%     and K is the number of PCs used.
%     - out.X_dot: out.X_dot is the time-derivative of X
%     - out.conditionMask: is a T x 1 mask where each element specifies
%     which condition each time point came from
%     - out.times_t1: is a T/timesStep x 1 vector that indicates which time
%     point correspond to each element in Q
%     - out.times_t2: is a T/timesStep x 1 vector that indicates which time
%     point was maximally tangled with each element in times_t1.
%     - out.allCombos: is a T/timesStep choose 2 x 2 matrix that
%     corresponds to each pair of times that tangling was computed for
%     (i.e. before the maximum over t2 was taken)
%     - out.q_allPairs: is a T/timesStep choose 2 x 1 vector that
%     corresponds to the tangling value associated with each pair of points
%     in out.allCombos.
%  
%
% USAGE EXAMPLES
%   [ Q, out] = tangleAnalysis( Data, .02, 'numPCs', 10 ); % to calculate tangling with 10 PCs; data collected every 20 ms
%   [ Q, out] = tangleAnalysis( Data, .001, 'cond2use', [1 4:6 10] ,'withinConditionsOnly',true); % to calculate tangling only on conditions #1,4,5,6,10 and only within conditions
%
% AUTHORSHIP
% |Author: Abigail Russo,|
% |Email: russoaa@gmail.com,|
% |Dated: November 2017|

function [ Q, out] = tangleAnalysis(  Data, samplingRate, varargin )
params = inputParser;

% Set default options
params.addParameter('numPCs',  8, @(x) isscalar(x) | isempty(x));
params.addParameter('fracVar', [], @(x) x < 1 & x > 0); % use as many PCs as is necessary to capture this fraction of variance in the data
params.addParameter('cond2use', 1:length(Data), @isvector); % default is to use all conditions
params.addParameter('withinConditionsOnly', false, @islogical); % if true, tangling will only be computed between pairs of times that come from the same condition
params.addParameter('timeStep', floor(.02/samplingRate), @isscalar); % larger values will speed computation
params.addParameter('softenNorm', [], @(x) isvector(x) | isempty(x));

% Parse command line for parameter settings
params.parse(varargin{:});

numPCs = params.Results.numPCs;
fracVar = params.Results.fracVar;
cond2use = params.Results.cond2use;
withinConditionsOnly = params.Results.withinConditionsOnly;
timeStep = params.Results.timeStep;
softenNorm = params.Results.softenNorm;

%% Format data

% Ensure data is formatted correctly
if size(Data(1).A,1) < size(Data(1).A,2)
   warning('Data(c).A should be a t x n matrix containing the firing rates for condition c. Ensure that this is the case')
end

% Unwrap all data into a ct X n matrix
A = [];
conditionMask = [];
for cc = 1:length(Data)
   theseData = Data(cc).A;
   A = [A; theseData];
   conditionMask = [conditionMask; cc*ones(size(theseData,1),1)];
end

% soft-normalize firing rates in the usual Churchland way
if ~isempty(softenNorm)
   normFactors = range(A,1)+softenNorm;
   A = bsxfun(@times, A, 1./normFactors);
elseif any(range(A,1) > 1) && isempty(softenNorm)
   warning('A should be normalized or soft-normalized such that the range of each neuron <= 1')
end
A = bsxfun(@minus, A, mean(A,1)); % mean center

%% Get masks to analyze the times and conditions of interest
timeMask = false(size(conditionMask));

for cc = 1:length(Data)
   thisCondition = conditionMask == cc;
   if ~isfield(Data(cc), 'times')
      Data(cc).times = 1:size(Data(cc).A,1);
   elseif isempty(Data(cc).times)
      Data(cc).times = 1:size(Data(cc).A,1);
   end
   
   if ~isfield(Data(cc), 'analyzeTimes')
      Data(cc).analyzeTimes = Data(cc).times;
   elseif isempty(Data(cc).analyzeTimes)
      Data(cc).analyzeTimes = Data(cc).times;
   end
   theseTimes = ismember(Data(cc).times,  Data(cc).analyzeTimes);
   timeMask_thisCondition = timeMask(thisCondition);
   timeMask_thisCondition(theseTimes) = true;
   timeMask(thisCondition) = timeMask_thisCondition;
end
myMask = timeMask & ismember(conditionMask, cond2use);

%% Reduce the data, get the derivative being careful to respect condition boundaries
A_masked = A(myMask,:);

% reduce to numPCs or however many PCs capture fracVar of the variance
[PCs, ~, v] = pca(A_masked);
V = cumsum(v)./sum(v);
if ~isempty(fracVar)
   % determine how many PCs to use
   v = cumsum(var(A_masked*PCs))./sum(var(A_masked*PCs));
   idx = 1:size(v,2);
   numPCs = idx(v > fracVar);
   numPCs = numPCs(1);
end
if numPCs > size(PCs,2)
   topPCs = PCs;
else
   topPCs = PCs(:,1:numPCs);
end
X = A_masked * topPCs;

% get derivatives separately for each condition
X_dot = zeros(size(X));
newCondMask = conditionMask(myMask);
conds = unique(newCondMask);
for cc = conds(:)'
   thisCondMask = newCondMask == cc;
   
   theseData = X(thisCondMask,:);
   derivThisCond = diff(theseData);
   
   % make sure the derivative and the projection are the same size so that
   % the same mask can be used
   derivThisCond = [derivThisCond(1,:); derivThisCond]; % pad so that the same size
   
   X_dot(thisCondMask,:) = derivThisCond;
end

%% Make tangling scale and sampling rate invariant
X_dot = X_dot / samplingRate; % correct for sampling rate
alpha = .1; % this will be added to the denominator to prevent hyperbolic growth
constant = alpha * sum(var(X)); % make proportional to variance so metric is approximately scale invarient

%% Make sure choice of timeStep is appropriate

%%% get time pairs to be compared
t = 1:timeStep:sum(myMask); % if this vector is large (greater than ~1000) this will take very long
timeEst = 0.5729*exp(size(t,2)*0.002873);

if size(t,2) < 100 
   warning('Few time points are being analyzed. Consider using more conditions or time points to get higher resolution')
end

if samplingRate*timeStep >= .05 % sampling every 50 ms
   warning('Temporal resolution is poor. Consider increasing the parameter timeStep to get higher resolution and prevent aliasing.')
end


if timeEst < 60
   fprintf('Estimated computation time: %1.0f seconds\n', timeEst)
elseif timeEst < 300 % 5 minues
   fprintf('Estimated computation time: %1.0f seconds\n', timeEst)
   warning('Computation time is long. Consider increasing the parameter timeStep')
elseif timeEst < 1800 % 30 minutes
   fprintf('Estimated computation time: %1.0f minutes\n', timeEst/60)
   warning('Computation time is very long!! Increase the parameter timeStep, decrease number of conditions, or shorten analysis window.')
else
   fprintf('Estimated computation time: %1.1f hours\n', timeEst/3600)
   warning('Computation time is very long!! Increase the parameter timeStep, decrease number of conditions, or shorten analysis window.')
end

%% COMPUTE TANGLING
allCombos = nchoosek(t,2);

% if desired, only consider time pairs that come from the same condition
if withinConditionsOnly
   
   useThisPair = false(size(allCombos,1),1);
   for currPair = 1:size(allCombos,1)
      t1 = allCombos(currPair,1);
      t2 = allCombos(currPair,2);
      if newCondMask(t1) == newCondMask(t2)
         useThisPair(currPair) = true;
      end
   end
   allCombos = allCombos(useThisPair,:);
end

%%% for each pair of times, calculate the tangling index (q)

% initialize
q = zeros(size(allCombos,1),1);

% caluclate q for all pairs
for currPair = 1:size(allCombos,1)
   t1 = allCombos(currPair,1);
   t2 = allCombos(currPair,2);
   
   % calculate the difference in derivative
   derivDiff =  norm(X_dot(t1,:) - X_dot(t2,:))^2;
   % calculate the difference in state
   stateDiff =  norm(X(t1,:) - X(t2,:))^2;
   
   % calculate the tangle index for every pair, q
   q(currPair) = derivDiff / (stateDiff + constant);
end

%%% get the max tangling index for each time point (Q)

% initialize
time_t1 = unique(allCombos(:));
Q = zeros(size(time_t1,1),1);
time_t2_Max = zeros(size(time_t1,1),1);

% caluclate Q
loopCount = 1;
for t = time_t1'
   time_t2 = allCombos(any(allCombos == t, 2),:);
   q_ti = q(any(allCombos == t,2));
   [Q(loopCount), idx] = max(q_ti);
   time_t2 = time_t2(time_t2 ~= t);
   time_t2_Max(loopCount) = time_t2(idx);
   loopCount = loopCount+1;
end

%% Save some data that can be useful for plotting and understanding tangling
out.X = X;
out.X_dot = X_dot;
out.conditionMask = newCondMask;
out.times_t1 = time_t1;
out.times_t2 = time_t2_Max;
out.allCombos = allCombos;
out.q_allPairs = q;

end