Pneumo-ABM
==========

Pneumo-ABM simulates an agent-based model of the dynamics of multiple pneumococcal serotypes and, optionally, a single strain of *Haemophilus influenzae*. Underlying the transmission model is a host life-history model that simulates birth, death, and the formation and destruction of households. The user can choose to ignore this demographic structure or to make transmission probabilities depend on host age and/or household membership. In addition, the model can simulate the effects of a pneumococcal vaccine targeting any subset of serotypes. It also includes a very simple, optional optimizer to fit the transmission rate to obtain a given prevalence in children.


Relevant citations
------------------
Cobey, S. and M. Lipsitch. 2012. Niche and neutral effects of acquired immunity permit coexistence of pneumococcal serotypes. *Science* 335(6074):1376-1380. doi: 10.1126/science.1215947

Cobey, S. and M. Lipsitch. 2013. Pathogen coexistence through hidden regimes of apparent competition. *American Naturalist* 181(1):12-24. doi: 10.1086/668598


Good things to know
-------------------
This software is not optimized for what it does. It was developed with different assumptions in mind. If you plan to use this code extensively, refactor. Certain features--like complex spatial structure in the host population--have been deprecated and mostly stripped from the code here. Some relics have been flagged in the Parameters.h file or have been left under the hood. You might wonder why they are there, if you start poking around. A few are left as inspiration. 


Requirements
------------
1. GNU C++ compiler (developed with gcc version 4.2.1)
2. Boost libraries (developed with boost 1.42; http://www.boost.org/)
3. (optional) Matlab (for viewing and processing results; http://www.mathworks.com/products/matlab/)


Compilation and execution
--------------------------
The simulation runs in C++. It is not multithreaded.

Sample compilation and execution on a Mac:

    g++ -O3 -o pneumo_abm.out -I /Applications/boost_1_42_0/ Host.cpp Simulation.cpp main.cpp Rdraws.cpp SimPars.cpp
    ./pneumo_abm.out

Nothing will happen if you do just this. Read the next section.

Outputs of individual simulations are summarized in Matlab by `ABM_processor.m` and further by `treatment_processor.m`.


Inputs
------
The simulation takes arguments from three sources: values inputted to the screen, values hard-coded in `Parameters.h`, and values in the input files.

###Screen inputs
Immediately after executing, the program expects (i) a numeric treatment identifier, (ii) a value of serotype-specific immunity (the "treatment", parameter sigma in Cobey & Lipsitch, 2012, *Science*), and (iii) a numeric simulation (replicate) identifier, which also serves as the seed for the PRNG. Practically, these inputs are best supplied via a bash script. One is included (`script.sh`) that runs multiple simulations (replicates) for multiple values of serotype-specific immunity (treatments) on a LSF-based cluster. If you use the script, you need to include two files:

1. `simulations.txt`: Lists the numeric identifiers of the simulations to be run for each treatment.
2. `treatments.txt`: Lists the values of serotype-specific immunity for which to run simulations in `simulations.txt`.


###Hard-coded inputs
Please see `Parameters.h`.

###File inputs
Sample input files have been included in the release. There is minimal validation of input files, so check yours carefully.

* `Associated_treatments.txt`<br>The values of serotype-specific immunity for which estimates of beta, the transmission rate, are given in `Betas_used.txt`. The program will compare its serotype-specific immunity (the treatment value) to the ones listed here and select the index whose value is closest. If you want to run a simulation with serotype-specific immunity of 0.3 and a particular transmission rate, add 0.3 to this file and the transmission rate in `Betas_used.txt` (at the same index). If you are fitting a new transmission rate, put one or more educated guesses in these files. 

* `Betas_used.txt`<br>Transmission rates corresponding to the values of serotype-specific immunity listed in `Associated_treatments.txt`. If `MATCH_PREVALENCE` is defined in `Parameters.h`, one of these will be the starting value of beta (the one whose associated treatment most closely matches the defined value of serotype-specific immunity). If `MATCH_PREVALENCE` is not defined, the simulation will be run with one of these values.

* `BIRTH_AGE_PMF.txt`<br>Probability mass function for age at giving birth. Used regardless of birth number. Indexed by age in years.

* `FLEDGE_PMF.txt`<br>Probability mass function for age at leaving household of origin and starting new household. Can be preempted by partnering. Indexed by age in years.

* `HFLU_PROBS.txt`<br>The probability that each serotype will be cleared in the presence of *H. influenzae*. Indexed by serotype.

* `IMMIG_RATES.txt`<br>Immigration rates for each serotype and *H. influenzae* (*w* in the Cobey & Lipsitch, 2012, *Science*). Indexed by serotype/strain.

* `INIT_AGE_PMF.txt`<br>The initial age distribution of the host population at the start of the simulation. Indexed by age in years.

* `INIT_INFECTEDS.txt`<br>The initial fraction of the host population colonized with each serotype (and *H. influenzae*, if simulated). Indexed by serotype/strain.

* `LSPAN_PMF.txt`<br>Probability mass function of host age at death. Indexed by age in years.

* `M_DURATION_INFECTION.txt`<br>The mean intrinsic durations of carriage (in days) of each serotype and *H. influenzae*. Indexed by serotype/strain.

* `NEIGHBORHOODS.txt`<br>Once upon a time, this file was used to specify relative contact rates between neighborhoods. It's still read in, but the code for >1 neighborhood has been mostly deprecated.

* `PAIR_PMF.txt`<br>Probability mass function of initiating partnership at each age. Indexed by age in years.

* `PARITY_PMF.txt`<br>The lifetime probability of having 0, 1, 2,… offspring. Indexed by number of offspring (starting at 0) for an *average* (50% male) individual.

* `WAIFW.txt`<br>Contact matrix ("who acquires infection from whom"), from column *j* to row *i*. Ignored unless age-assortative transmission is turned on (`\\#DEFINE NO_AGE_ASSORT` in `Parameters.h`). Will be normalized. Values of zero are reset to just above zero (in `SimPars.cpp`) to ensure that two hosts in the same household always have a nonzero probability of transmission, in case household transmission is also turned on (`\\#DEFINE NO_HHOLD`).

###Modifying model assumptions 
To simulate random mixing, define `NO_AGE_ASSORT` and `NO_HHOLDS` in `Parameters.h`. To add vaccination, define `SIM_PCV` in the same file. To not simulate *H. influenzae*, ensure that the last "serotype" in the input files has an immigration rate of 0 and no initial infecteds. 


Outputs
-------
Most simulation outputs are saved to files.

###To screen 
An abbreviated summary of model assumptions, key inputs, progress, and total simulation time are printed to screen. If fitting the transmission rate, this information is shown for every simulation. 

###To files 
All output files are prefixed by the treatment identifier *X* and the replicate/simulation number *Y*, `tr_X_sim_Y_...`).

* `age_dist_neighborhood_0`<br>The current number of hosts of each age (in y, column) for each of the times in `dem_times` (row).

* `BETA`<br>The transmission rates used for each serotype/strain.

* `coinfection_dist_hflu-pneumo`<br>The fraction of kids <5 y carrying *H. influenzae* also colonized with each serotype (column) for every time in `epid_times` (row).

* `coinfection_dist_pneumo`<br>The number of kids <5 y colonized with 1,2,… serotypes of pneumococcus (column) for every time in `epid_times` (row).

* `dem_times`<br>The times (in days) at which demographic data are output.

* `epid_times`<br>The times (in days) at which epidemiological data are output.

* `hh_dist`<br>The distribution of households containing 1,2,… members (column) for every time in `dem_times` (row).

* `infecteds_i_neighborhood_0`<br>The number of people in each age (in y, column) carrying serotype *i* for every time in `epid_times` (row).

* `infections_i_neighborhood_0`<br>The number of current colonizations with serotype *i* in hosts of each age (in y, column) for every time in `epid_times` (row).

* `theta`<br>Lists, at the end of the simulation, the age (in days) of each child <=5 y old and the cumulative number of cleared pneumococcus colonizations per child.

* `totCarriage`<br>Total number of hosts carrying pneumococcus at every time in `epid_times`.

* `XI`<br>The matrix of serotype-specific (cross-)immunity between each serotype *i* and *j*. As only homologous serotype-specific immunity is currently included in the model, only the diagonals might be nonzero.   

###Processing results
The Matlab files produce a few self-explanatory plots summarizing features of the simulations. Run `treatment_processor.m`, which calls `ABM_processor.m`.


License
--------
Copyright © 2012, Sarah Cobey

The software and its documentation are in the public domain and are furnished "as is." The author makes no warranty, express or implied, as to the usefulness of the software and documentation for any purpose. The author assumes no responsibility for the use of the software and documentation or to provide technical support to users. 


Support
-------
There is no support. I will nonetheless try to reply to reasonable requests for clarification. Email scobey@hsph.harvard.edu.
