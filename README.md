# Software of PONE-D-19-31181 DOI: 10.1371/journal.pone.0226094
Software running the simulations made for PONE-D-19-31181 manuscript


The simulations are divided in 4 folders:

  ANLS:             the simulation of the model in the ANLS case during normal conditions
  
  NALS:             the simulation of the model in the NALS case during normal conditions
  
  ANLS_exrsice:     the simulation of the model in the ANLS case during exercise conditions
  
  NALS_exrsice:     the simulation of the model in the NALS case during exercise conditions


Each folder contains a script file in order to run the simulation called script.sh.
The main program is called SUB_main.f and the Makefile contains the appropriate flags for compiling and the connection of the subroutines.
The rest files are subroutines called in SUB_main.f or other subroutines.

The file Com_sol.f in ANLS_exercise and NALS_exercise runs independently after the run of the simulation with script.sh. This file needs the solution under normal condition (Asol.dat, which is obtained from the simulation of ANLS and NALS folders) and the solution under exercise conditions (Asol.dat, which is obtained from the simulation of ANLS_exercise and NALS_exercise folders) in order to compare them and generate Figures 9 and 11.
