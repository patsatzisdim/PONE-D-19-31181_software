# PCOMPBIOL-D-19-00752R1_software
Software running the simulations made for PCOMPBIOL-D-19-00752R1 manuscript


The simulations are divided in 4 folders:
  ANLS:             the simulation of the model in the ANLS case during normal conditions
  NALS:             the simulation of the model in the NALS case during normal conditions
  ANLS_exrsice:     the simulation of the model in the ANLS case during exercise conditions
  NALS_exrsice:     the simulation of the model in the NALS case during exercise conditions

Each folder contains a script.sh file in order to run the simulation
The main program is called SUB_main.f and the Makefile contains the appropriate flags for compiling and the connection of the subroutines.
The rest files are subroutines called in SUB_main.f or other subroutines.
