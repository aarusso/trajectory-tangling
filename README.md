# trajectory-tangling
MATLAB code, Demo, and sample data to analyze data for trajectory tangling as in Russo et al Neuron 2018 (https://www.sciencedirect.com/science/article/pii/S0896627318300072?via%3Dihub)

# Datasets
**M1_sampleData.mat**: trial-avergaged firing rates from 116 neurons

**EMG_sampleData.mat**: trial-avergaged activity from 29 muscles

**kinematic_sampleData.mat**: trial-avergaged kinematic data corresponding to a variety of behavioral parameters such as hand position and velocity 

The data included in this demo were recorded in the Churchland lab at Columbia University. Single-unit neural data were recorded sequentially from a macaque monkey during the performance of a hand-pedaling task.  Monkeys grasped a hand crank and cycled through a virtual environment for a number of prescribed cycles as indicated by a visual cue. The data included here correspond to two seven-cycle conditions. EMG data were recorded with intramuscular electrodes and are rectified and filtered to obtain an envelop of the muscle response. Neural and EMG were filtered with a Gaussian kernel with a standard deviation of 25ms and trial averaged. For more information about data processing and the task, see Russo et al Neuron 2018. 

Data should be formatted as for jPCA. Each data structure (e.g. D_m1, located in M1_sampleData.mat) contains C elements corresponding to the number of conditions (here, 2). The t x n matrix 'A' contains the trial-averaged firing rates as a function of time for each neuron. All other fields are optional. The 'times' field indicates the time course of each sample in A and the 'analyzeTimes' field indicates which of these times should be analyzed for tangling.

# Code

**tangleAnalysis.m**

**tangling_demo.m**

**tangle_visualize.m**

**Supporing functions**: makes pretty scatter plots
Here is some code to aid visualization of tangling on the sample
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
