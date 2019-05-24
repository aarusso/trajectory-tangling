# trajectory-tangling
MATLAB code, demo, and sample data to analyze data for trajectory tangling as in [Russo et al Neuron 2018](https://www.sciencedirect.com/science/article/pii/S0896627318300072?via%3Dihub)

**Author**: Abigail Russo, May 2019

# Datasets
The data included in this demo were recorded in the Churchland lab at Columbia University. Single-unit neural data were recorded sequentially from a macaque monkey during the performance of a hand-pedaling task.  Monkeys grasped a hand crank and cycled through a virtual environment for a number of prescribed cycles as indicated by a visual cue. The data included here correspond to two seven-cycle conditions. EMG data were recorded with intramuscular electrodes and are rectified and filtered to obtain an envelop of the muscle response. Neural and EMG were filtered with a Gaussian kernel with a standard deviation of 25ms and trial averaged. For more information about data processing and the task, see [Russo et al Neuron 2018](https://www.sciencedirect.com/science/article/pii/S0896627318300072?via%3Dihub).

Data should be formatted as for [jPCA](https://www.dropbox.com/sh/2q3m5fqfscwf95j/AAC3WV90hHdBgz0Np4RAKJpYa?dl=0&preview=NOTES.pdf). Each data structure (e.g. D_m1, located in **M1_sampleData.mat**) contains C elements corresponding to the number of conditions (here, 2). The t x n matrix 'A' contains the trial-averaged firing rates as a function of time for each neuron. All other fields are optional. The 'times' field indicates the time course of each sample in A and the 'analyzeTimes' field indicates which of these times should be analyzed for tangling. For more information on formatting, see the **tangleAnalysis** code comments.

* **M1_sampleData.mat**: trial-avergaged firing rates from 116 neurons
* **EMG_sampleData.mat**: trial-avergaged activity from 29 muscles
* **kinematic_sampleData.mat**: trial-avergaged kinematic data corresponding to a variety of behavioral parameters such as hand position and velocity 

# Code

* **tangleAnalysis.m**: compute trajectory tangling on data formatted as described above
* **tangling_demo.m**: visualize task kinematics, compute tangling on sample datasets, play around with parameters, visualize tangling in two dimensions
* **tangle_visualize.m**: aids visualization of tangling on the sample datasets. 
* **formattedScatter.m**: makes pretty scatter plots

*Supporing functions*: **AxisMMC.m** (makes pretty axes, used by **formattedScatter**), **arrowMMC.m** (makes pretty arrows, used by **tangle_visualize**)

