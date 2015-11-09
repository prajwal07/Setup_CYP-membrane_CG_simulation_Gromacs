# Setup_CYP-membrane_CG_simulation_Gromacs
Bash script to set-up and Run the CYP-membrane Coarse-grain MD simualtion using Martini FF and Gromacs Software

Please provide the default name you want to use as Prefix for all files 
Command ./run_cg_simulation.sh input_PDB_filename prefix_for_files 
Example: ./run_cg_simulation.sh CYP1A1.pdb cg_CYP_ori1 

Please set the path for GROMACS executables 
"Example export PATH="/hits/fast/mcm/app/gromacs/gromacs5.0.6/compiled/bin":$PATH 

Requirements:
1) Gromacs sfotware installed at PATH.
2) DSSP tool for get seconadary structure elements data file from PDB file.
3) Full length all-atom PDB structure of CYP protein.
4) VMD should be installed and accessible via command "vmd" 

Sample data files for test run:
Sample files are present in "test" folder

Requried input files/scripts:
Input files/scripts are present in "input_files" folder

Executable script is "run_cg_simulation.sh"
