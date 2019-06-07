# Setting up your IDE #

## Compiler and linker flags ##
Please check the provided Makefile (located in the main directory) to setup your IDE with the required compiler flags.

## Configfile paths ##
Each IDE uses their own path for creating an output file. Please make sure that all the configuration files in the git resourse bin directory are coppied to your IDE output file path.

## Create a directory for your results ##
You need to create 1 directory and 3 subdirectories next to your output file:

outputfile.out

/Results/		xyz trajectories of the main simulation are stored here

/Results/Relaxation/	xyz of the initial relaxation of the system will be saved here.

/Results/Reports/	These are the information the code uses to run the simulation. This file can be later used to generate a configuraion file (if the original is lost).

/Results/Resumes/	The VC can use these files to resume any interupted simulation.

