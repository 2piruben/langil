# Analysis of mRNA distributions in a proliferating cell population

Analysis done in *R. Perez-Carrasco, C. Beentjes, R. Grima, Effects of cell cycle variability on lineage and population measurements of mRNA abundance. BioRxiv (2019)* to study the distribution of mRNA in a cell analyzing effects of cell cycle duration variability in a transcription model that takes into account bipartition, DNA replication, and a structured cell cycle.

## Files

* [cyclo_2states.cpp](cyclo_2states.cpp) source code to run trajectories (lineage) of an individual cell using the `langil` class
* [cyclo_2states_pop.cpp](cyclo_2states_pop.cpp) source code to run trajectories of a proliferating population of cells using the `langil` class and stores the final snapshot.
* [runcyclo.py](runcyclo.py) python code defining functions to create input files and run `cyclo_2states` and `cyclo_2states_pop`
* [manuscriptfigures.py](manuscriptfigures.py) python code used to create the figures of the manuscript 
* [bursttrans.py](bursttrans.py) definition of the analytical expressions for the factorial moments of the distribution of mRNAs


