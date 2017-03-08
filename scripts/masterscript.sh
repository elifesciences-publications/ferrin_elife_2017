# Author: Arvind R. Subramaniam

# this script uses pre-run simulation data; you need to uncomment
# two highlighted lines below to run fresh simulations (takes several hours)

# compile whole-cell simulation
g++ -g -w wholecell_simulation.cpp -o wholecell_simulation -lboost_regex
date; echo "Compiled whole cell simulation successfully."

# compile reporter simulation
g++ -g -w reporter_simulation.cpp -o reporter_simulation -lboost_regex
date; echo "Compiled reporter simulation successfully."

# create directories for storing input mrna and stall elongation rate files
mkdir --parents ../annotations/simulations
cd ../annotations/simulations
for runnumber in `seq 8; seq 11 16`;
do
  mkdir --parents run$runnumber;
done;
cd ../../scripts

# create directories for writing results from 16 simulation runs
mkdir --parents ../rawdata/simulations
cd ../rawdata/simulations
for runnumber in `seq 9; seq 11 16`;
do
  mkdir --parents run$runnumber;
done;
cd ../../scripts

# create directories for writing processed data and figures
mkdir --parents ../processeddata/platereader
mkdir --parents ../processeddata/simulations
mkdir --parents ../figures
date; echo "Made directories for storing generated data and figures."

# plot figure 1 fluorescence, polysome qpcr, mrna qpcr data
jupyter nbconvert --log-level=ERROR --to notebook --execute plot_fig1_fluorescence.ipynb
date; echo "Plotted Fig 1 fluorescence data."
jupyter nbconvert --log-level=ERROR --to notebook --execute plot_fig1_polysome_qpcr.ipynb
date; echo "Plotted Fig 1 polysome qPCR data."
jupyter nbconvert --log-level=ERROR --to notebook --execute plot_fig1_mrna_qpcr.ipynb
date; echo "Plotted Fig 1 mRNA qPCR data."

# analyze plate reader data to use as input for simulations and plotting
jupyter nbconvert --log-level=ERROR --to notebook --execute plate_reader_analysis.ipynb
date; echo "Analyzed all plate reader data."

# run whole-cell simulation and parameter sweep of reporter simulation.
# this submits jobs to a SLURM cluster. .
# Comment out indicated lines in the jupyter notebook if running locally.
# Uncomment the simulation run below if you want to generate fresh
# simulation data
#################################################################################################################
# jupyter nbconvert  --log-level=ERROR --to notebook --execute run_simulations_whole_cell_parameter_sweep.ipynb #
#################################################################################################################
# note that the commands below will not run till the simulations from the 
# previous step are complete

# fit simulation results to experiment data.
# if you ran fresh simulations, you need to uncomment a few lines in 
# this notebook to read the fresh data
jupyter nbconvert --log-level=ERROR --to notebook --execute fit_simulation_to_experiment.ipynb
date; echo "Fitted simulation results to experimental data."

# create new simulation inputs based on experiment fits of previous simulations
jupyter nbconvert --log-level=ERROR --to notebook --execute simulation_inputs_based_on_experiment_fits.ipynb
date; echo "Created input files for simulations based on experimental fits."

# run reporter simulation based on experimental fits.
# this submits jobs to a SLURM cluster. 
# Comment out indicated lines in the jupyter notebook if running locally.
# Uncomment the simulation run below if you want to generate fresh
# simulation data
##############################################################################################################
# jupyter nbconvert --log-level=ERROR --to notebook --execute run_simulations_based_on_experiment_fits.ipynb #
##############################################################################################################

# note that the commands below will not run till the simulations from the 
# previous step are complete 

# plot simulation results along with experiment comparison.
# if you ran fresh simulations, you need to uncomment a few lines in 
# this notebook to read the fresh data
jupyter nbconvert  --log-level=ERROR --to notebook --execute plot_simulation_results_figs_3_to_7.ipynb
date; echo "Plotting simulation results and comparison to experiments."

# remove nbconvert copies of jupyter notebooks
rm *.nbconvert.ipynb
