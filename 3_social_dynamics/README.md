## Monte Carlo simulations
1. First run Matlab script "segregation.m" to obtain the minimum of rho at equilibrium. We will use this equilibrium as stop criteria for Monte Carlo simulations.
2. The Matlab scripts that solves the population density model are stored in directory "cluster_computation". The names of the sub-directories indicate the noise amplitude used. The Matlab scripts that solves the population density model all have the same name "segregation_with_noise_ETD.m". There are two key parameters: 1) Nevent, after how many events (segregation) the simulation is stopped, 2) seed, seed used to generate random numbers. For small noise amplitude, the waiting time is long, so to reduce computational time, we use a smaller Nevent and run the same script with different seed simulteneously.
3. Calculate the mean waiting time for segregation.
## Prefactor of Eyring-Kramers Law
1. Run Matlab script "GAD_SEG_JB_const_mob.m" to calculate the saddle shape using Gentlest Ascent Dynamics and Exponential Differencing method. The saddle shape is stored in file "s_128.mat".
2. Run Matlab script "plot_EK.m" to calculate the prefactor and make initial plot.
3. Run Python script "seg-MC.py" for the final plot.
## Plot segregation process
1. Run "segregation_with_noise_ETD.m". This is modified to get one event and save one trajectory to "SEG_process.txt".
2. Run "seg_process.py" for plot.

