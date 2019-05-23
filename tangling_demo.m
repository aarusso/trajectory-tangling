%% TANGLING DEMO
%
% Data formatting
% For tangling and divergence analyses, data should be formatted as for
% jPCA. Each data structure (e.g. D_m1, loacted in M1_sampleData.mat)
% contains C elements corresponding to the number fo conditions (here, 2).
% The t x n matrix 'A' contains the trial-averaged firing rates as a
% function of time for each neuron. All other fields are optional. The
% 'times' field indicates the time course of each sample in A and the
% 'analyzeTimes' field indicates which of these times should be analyzed
% (e.g. by jPCA, tangling, or divergence).
% 
% Sample data contents 
% The data included in this demo were recorded in the Churchland lab at
% Columbia University. Single-unit neural data were recorded sequentially
% from a macaque monkey during the performance of a hand-pedaling task.
% Monkeys grasped a hand crank and cycled through a virtual environment for
% a number of prescribed cycles as indicated by a visual cue. The data
% included here correspond to two seven-cycle conditions. EMG data were
% recorded with intramuscular electrodes and are rectified and filtered to
% obtain an envelop of the muscle response. Neural and EMG were filtered
% with a Gaussian kernel with a standard deviation of 25ms and trial
% averaged. For more information about data processing and the task, see
% Russo et al Neuron 2018. 
%
%
% AUTHORSHIP
% |Author: Abigail Russo,|
% |Email: russoaa@gmail.com,|
% |Dated: May 2019|
%% Plot example behavioral kinamtics
% The sample data contain two conditions: seven-cycle for forward pedaling
% starting at the bottom and seven cycle for backward pedaling starting at
% the bottom. Plot the kinematic data (hand position and velocity, torque,
% and wrist angle) from 8 example sessions to get a feel for the behavior.
% Data have been trial averaged and aligned as described in Russo et al
% 2018. The analysis time window is 100ms before movement onset to 100ms
% after movement offset.

load kinematic_sampleData


T = length(p(1).A);
ticks = sort([500 1500:500:T T-500]);
tickLabels = arrayfun(@(n) {num2str(n)}, (ticks/1000) - 1.5); % time from movement onset in seconds
tickLabels{1} = 'target on';
tickLabels{end-1} = 'stop';

figure; 
subplot(421); hold on;
h = plot(p(1).A,'color',[.2 .8 .2]);
j = plot(p(2).A,'color',[.8 .2 .2]);
ylabel('world position')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])
analyzeTimes = p(1).analyzeTimes;
k = plot([analyzeTimes(1) analyzeTimes(end)],[0 0],'k');
plot([analyzeTimes(1) analyzeTimes(1)],[-1 1],'k');
plot([analyzeTimes(end) analyzeTimes(end)],[-1 1],'k');
legend([h(1) j(1) k],{'Forward bottom start','Backward bottom start','Analysis time window'},'Location','NW')

subplot(423); hold on;
plot(diff(p(1).A),'color',[.2 .8 .2])
plot(diff(p(2).A),'color',[.8 .2 .2])
ylabel('angular hand velocity')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])

subplot(425); hold on;
plot(tq(1).A,'color',[.2 .8 .2])
plot(tq(2).A,'color',[.8 .2 .2])
ylabel('torque')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])

subplot(427); hold on;
plot(w(1).A,'color',[.2 .8 .2])
plot(w(2).A,'color',[.8 .2 .2])
ylabel('wrist angle')
set(gca,'XTick',ticks,'XTickLabels',tickLabels,'XLim',[0 T])
xlabel('time from movement onset (s)')


subplot(422); hold on;
plot(vp(1).A,'color',[.2 .8 .2])
plot(vp(2).A,'color',[.8 .2 .2])
ylabel('vertical hand position')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])

subplot(424); hold on;
plot(hp(1).A,'color',[.2 .8 .2])
plot(hp(2).A,'color',[.8 .2 .2])
ylabel('horizontal hand position')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])

subplot(426); hold on;
plot(vv(1).A,'color',[.2 .8 .2])
plot(vv(2).A,'color',[.8 .2 .2])
ylabel('vertical hand velocity')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])

subplot(428); hold on;
plot(hv(1).A,'color',[.2 .8 .2])
plot(hv(2).A,'color',[.8 .2 .2])
ylabel('horizontal hand velocity')
set(gca,'XTick',ticks,'XTickLabels',[],'XLim',[0 T])
set(gca,'XTick',ticks,'XTickLabels',tickLabels)
xlabel('time from movement onset (s)')

suptitle('Example kinematics from 8 recording sessions')

%% Analyze tangling
load M1_sampleData
load EMG_sampleData
samplingRate = 1/1000; % inherent to the data. Do not change


% parameters to play with
% note: these parameters and will be chosen automatically if not specified

numPCs = 2; % pick a number that will capture most of the variance
% note: if you call as follows, tangleAnalysis will choose how many PCs to use based on how much variance you want captured 
% Q_m1 = tangleAnalysis(D_m1, samplingRate,'softenNorm',5,'fracVar',.9); 

timeStep = 20; % smaller values will yield high resolution at the expense of computation time, default will sample at 20ms
% code will complain if you pick a value that it isn't happy with

withinConditionsOnly = false; % if true, will only analyze tanlging for times within the same condition

% run analysis
[Q_m1, out_m1] = tangleAnalysis(D_m1, samplingRate,'softenNorm',5 ,'timeStep',timeStep,'withinConditionsOnly',withinConditionsOnly); % soft normalize neural data
[Q_emg, out_emg] = tangleAnalysis(D_emg, samplingRate,'softenNorm',0,'timeStep',timeStep,'withinConditionsOnly',withinConditionsOnly); % fully normalize EMG

formattedScatter(Q_m1, Q_emg, {'Q_{M1}','Q_{EMG}'});

%% Visualize tangling in EMG
% Here is some code to aid visualization of tangling on the sample
% datasets. The first two dimensions of neural trajectory (X) are plotted
% verses one another. A time point with high tangling (t1) was chosen and
% is plotted as an arrow in blue. In black arrow is another time point
% (t2). The direction of the arrow indicates the direction of the
% derivative in two dimensions at that time point while the size of the
% arrow represents the magnitude of the derivative (e.g. the speed).
% Tangling between two pairs of points will be high if the locations of t1
% and t2 are similar but the derivatives are different in terms of
% maginutude or direction. Note that these plots are in two dimensions but
% tangling was caluclated on however many dimensions was specified above
% and that a constant was added to the denominator to prevent hyperbolic
% growth.
%
% Press (or hold) the down arrow to step through t2 and see how tangling
% changes as X and X_dot change with respect to one another. If
% withinConditionsOnly was set to false, you will see t2 from the same
% condition as well as t2 from the other condition (just keep steping
% through t2). The visualization will pause briefly when tangling becomes
% particularly high. Press CTR+C to exit early if you like, otherwise it
% will end when it has stepped through all t2s. Note that I wrote this code
% specifically for the demo so some things have been hard coded (e.g. what
% is considered "particularly high" tangling) and I can't guarentee it will
% look super pretty for all datasets. Hopefully, it is useful for gaining
% intuition for the metric and exactly what it measures. (read: apologies
% for the jankiness).
%

tangle_visualize( out_emg ) 


%% Visualize tangling in M1

tangle_visualize( out_m1 ) 

