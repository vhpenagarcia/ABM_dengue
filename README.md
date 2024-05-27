# ABM_dengue
Agent-based model to assess the relative contribution of different urban environments in dengue transmission

There are three folders containing codes and data. Those are organizing as:

•	populations: 
Contains the information and codes to create the synthetic populations. The files ages_cities.csv, hrs_times.csv and houses_containers.csv are the data related to the age probabilities, spending time range for each structure type and frequencies of presence for different container sizes, respectively. The files pop_parameters.jl includes the codes to call the necessary information and its preparation for the creation of populations. File functions_pops.jl includes functions necessary to execute codes included in the file populations.jl, which creates the populations.

•	model: 
Contains the agents-based model, which core code can be found in the file parallel_abm_model.jl. The files human_pop_20k.csv and vector_pop_20k.csv include information of human and vector populations, respectively. The file weather2.csv includes the daily weather information. Finally, the functions needed to execute the models are found in source_info_movimm.jl.
In addition, a folder named vector_control is included comprising two code files as example for preventive (parallel_abm_preventive.jl) and reactive (parallel_abm_reactive.jl) scenario simulations. The shared examples includes a code for preventive control when it is applied on small containers in households by eliminating 100% of containers, and a code for reactive control when it is applied on day 50 by eliminating 75% of small containers on both household and non-household environments.

•	results: 
Includes information produced by the model. The model without any control can be found in the file infections_100.csv and infections_100c.csv (the latter includes information extracted from children). For control strategies many simulations were executed and hence many files produced. The file summary.jl is a code to call all the simulations (two files with 200 replicates each for a single simulation), create summaries and putting them together in a single summary file. A file was created for preventive strategies (summary_preventive_iqr.csv) and another for reactive strategy results (summary_reactive_iqr.csv). To execute summary.jl is necessary to have a control database file, those are included in the folder as Control_preventive_ABM.csv and Control_reactive_ABM.csv for preventive and reactive summaries, respectively.

# ABM_urban
This folder contains a modification of the original model to make it spatially explicit and answer the questions related to influence of the distribution of non-household environments on the burden of dengue. Additionally, this question was evaluated along with the effect of movement of people and vectors.
To evaluate different distance regimes, three population files were created each of them represented one distance regime out of three: closest, at least 500 meters, and at least 1000 meters. The code to create these files from the original population file and the file of the mosquitoes can be found in the code spatial_assignment.jl. On the same file can be found the code to create the file with the structures and the spatial coordinates. Three different spatial conformations were evaluated: 
* Scattered: when non-household environments are randomly distributed in space. In the file was coded as random 
* Centered: When most of non-household environments are concentrated in a centered cluster
* Clustered: When most of the non-household environments are concentrated in three different clusters, coded as scat in the file

The level of movement for humans and vectors are designated in the code. Originally evaluated in three levels:

For humans: 100% (with no modification from the original model, the file abm_urban_rand.jl is the code for this type of simulations), 50%, and 20% of human movement. The latter levels are specified in the same code. The file abm_urban_rand50.jl is coded to simulate 50% of movement. The value can be seen/modified in lines 68 and 124.

For mosquitoes: 100%, 50%, or 10% (default), which can be specified in the code in line 53.

The distance in the code files in the folder are simulating a distance of at least 1000 meters. To change the file to include other distance regime (like the closest or at least 500 meters), the proper file should be specified in the code (under @everywhere statement, defining the population object).
