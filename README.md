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
