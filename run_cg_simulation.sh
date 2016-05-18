#!/bin/bash

#### Please provide the default name you want to use as Prefix for all files ####
#### Command ./run_cg_simulation.sh input_PDB_filename prefix_for_files #############
##### Example: ./run_cg_simulation.sh CYP1A1.pdb cg_CYP_ori1 #####

#### Please set the path for GROMACS executables ##########
### "Example export PATH="/hits/fast/mcm/app/gromacs/gromacs5.0.6/compiled/bin":$PATH ####

export PATH="/hits/fast/mcm/app/gromacs/gromacs5.0.6/compiled/bin":$PATH

input_files_PATH=/hits/fast/mcm/nandekpl/software/Setup_CYP-membrane_CG_simulation_Gromacs/input_files;
echo $input_files_PATH;
echo "Copy necessary files from directory $input_files_PATH";

cp $input_files_PATH/martini* .;

export in_pdb=$1;
export fname=$2;

echo "Dear User your files will generate with prefix name $1 "
echo "Running MARTINIZE.PY script on all-atom protein-protein model with elastic network ON"
echo "Input file is CYP.pdb and Outputs are cg_CYP.top and cg_CYP.pdb"

./martinize.py -f $in_pdb -o cg_CYP.top -x cg_CYP.pdb -dssp /sw/mcm/app/dssp/2002/dssp -p backbone -ff martini22 -elastic -ef 500 -el 0.5 -eu 0.9 -ea 0 -ep 0;

cat << EOF > superimpose.vmd
mol new $input_files_PATH/CYP51_cg.pdb
mol new cg_CYP.pdb

set s0 [ atomselect 1 "all" ]
# set s1 [ atomselect 0 "name BB SC1 and resname PRO and resid 29" ]
# set s2 [ atomselect 1 "name BB SC1 and resname PRO and resid 2" ]

set s1 [ atomselect 0 "name BB and resid 2 to 21 and not resname W WF POP" ]
### CHANGE THE RESIDUE NUMBERS FOR TM HELIX ###

set s2 [ atomselect 1 "name BB and resid 11 to 30" ]

set mat [ measure fit \$s2 \$s1 ]
\$s0 move \$mat

##set r [ measure rmsd \$s2 \$s1 ]
##puts "RMSD is \$r"

set s3 [atomselect 0 "resname POP"]

\$s0 writepdb CYP_cg_si.pdb
\$s3 writepdb cg_lipids.pdb

quit

EOF

/sw/mcm/app/vmd/vmd-1.9.1-amd64/vmd -e superimpose.vmd -dispdev none;


cat CYP_cg_si.pdb cg_lipids.pdb > ${fname}.pdb;

sed -i 's/END/TER/g' ${fname}.pdb;

echo "Check the orientation membrane bound model PDB file Example: pymol ${fname}.pdb ";
echo "";

echo "Solvate the membrane-bound CYP system in a water box with the same x and y dimension as the membrane";
echo "Define BOX SIZE";
editconf -f ${fname}.pdb -box 14.3912 13.9582 20 -o ${fname}.box.pdb -center 7.1210  6.9067 7

echo "Add water box of Non-polarizable water model";

gmx solvate -cp ${fname}.box.pdb -cs $input_files_PATH/waterMix-1bar-303K.gro -o ${fname}.w.pdb -radius 0.21
    
echo "Check generated PDB files in pymol Example: pymol ${fname}.w.pdb ";
echo "";

echo "REMOVE WATER molecules";

cat << EOF > remove_water.vmd
mol new ${fname}.w.pdb

### Change the protein residue numbers as per protein ###
set rem_wat [atomselect 0 "(resname W WF) and within 5 of (resname POP or (not resname W WF POP and resid 1 to 513))"]
\$rem_wat num

set resid [\$rem_wat get resid]

set system [atomselect top "not ((resname W WF) and (within 5 of (resname POP or (not resname W WF POP and resid 1 to 513))))"]
\$system writepdb ${fname}.w_removed.pdb

quit

EOF

/sw/mcm/app/vmd/vmd-1.9.1-amd64/vmd -e remove_water.vmd -dispdev none;


echo "Check generated PDB file in pymol Example: pymol ${fname}.w_removed.pdb";
echo "";

echo "Make index file";
make_ndx -f ${fname}.w_removed.pdb -o ${fname}.w_removed.no.ion.ndx << EOF
q
EOF

echo "Count number of water atoms";
grep "W " ${fname}.w_removed.pdb | wc -l > temp;
read w < temp;
echo "Number of W water residues: $w";
rm temp;
grep "WF" ${fname}.w_removed.pdb | wc -l > temp;
read wf < temp;
echo "Number of WF water residues: $wf";
rm -r temp;

/usr/bin/python $input_files_PATH/remove_elastic_network.py Protein.itp Protein_EN_bead_removed.itp;

### Check for "cg_CYP1A1_CPR_cl.no.ion.top" file in folder ###
cat << EOF > ${fname}.no.ion.top
#include "martini_v2.2.itp"     
#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_ions.itp"    
#define RUBBER_BANDS  
#include "Protein_EN_bead_removed.itp"

[ system ]
; name
Martini system from CYP.pdb

[ molecules ]
; name        number
Protein          1
POPC            594
W               $w
WF              $wf

EOF


echo "Create folder for minimization";
mkdir ${fname}_min;
echo "Copy minimization input file";
cp $input_files_PATH/NPW_min.mdp ${fname}_min/;

echo "Check system with GROMPP"
grompp -f ./${fname}_min/NPW_min.mdp -c ${fname}.w_removed.pdb -p ${fname}.no.ion.top -n ${fname}.w_removed.no.ion.ndx -o ${fname}.w_min_0.tpr -maxwarn 1;

echo "Create copy of recent toplogy file";
cp ${fname}.no.ion.top ${fname}.top;

### Check the Total charge of system ###

echo "";
echo "CHECK Total charge on system";
echo "Add Counter IONS in PDB file to Neutralize it";
genion -s ${fname}.w_min_0.tpr -n ${fname}.w_removed.no.ion.ndx -p ${fname}.top -o ${fname}.ion.pdb -neutral -pname NA+ -nname CL- << EOF
14
EOF

echo "Modify the PDB file to change the atom name from NA to NA+"
sed -i 's/NA  NA+/NA+ NA+/g' "${fname}.ion.pdb";
sed -i 's/CL  CL-/CL- CL-/g' "${fname}.ion.pdb";

make_ndx -f ${fname}.ion.pdb -o ${fname}.ion.ndx << EOF
chain A
q
EOF

echo "Go to minimization directory";
cd ${fname}_min;
i
echo "Run Gromacs commands";
make_ndx -f ../${fname}.ion.pdb -o ${fname}_min.ndx << EOF
q
EOF

grompp -f NPW_min.mdp -c ../${fname}.ion.pdb -p ../${fname}.top -n ${fname}_min.ndx -o ${fname}_min.tpr -maxwarn 1;
mdrun -deffnm "${fname}_min";

echo "First minimization finished";
echo "Second minimization started";
grompp -f NPW_min.mdp -c ${fname}_min.gro -p ../${fname}.top -n ${fname}_min.ndx -o ${fname}_min2.tpr;
mdrun -deffnm "${fname}_min2";
echo "Second minimization finished";


echo "Create folder for Equilibration";
cd ..;
mkdir ${fname}_eq;
mkdir ${fname}_md;

cd ${fname}_md;
cp $input_files_PATH/connect_pdb_atoms.pl .;
cp $input_files_PATH/NPW_md.mdp .;
cp $input_files_PATH/calculate_parameters.vmd .;

cat << EOF > md_${fname}_submit.sh
#!/bin/bash
#$ -S /bin/bash
#$ -pe mvapich16 96
#$ -q bridge.q,bridge_fat.q
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00

source /etc/profile.d/modules.sh
module load mbm
module delete gcc
module load sge
module load gromacs/5.0.5/fftw/3.3.3/gcc/4.9.2/openmpi/1.8.5/bridge

# Run command
fname="${fname}";

mpirun mdrun_mpi -deffnm "${fname}_md" -maxh 24;

exit

EOF

chmod 777 md_${fname}_submit.sh;

cd ../${fname}_eq;
cp $input_files_PATH/NPW_eq.mdp .;


echo "Run Gromacs commands";
editconf -f ../${fname}_min/${fname}_min2.gro -o ${fname}_min2.pdb;

make_ndx -f ${fname}_min2.pdb -o ${fname}_eq.ndx << EOF
q
EOF

grompp -f NPW_eq.mdp -c ${fname}_min2.pdb -p ../${fname}.top -n ${fname}_eq.ndx -o ${fname}_eq.tpr -maxwarn 1 > grompp.log;


cat << EOF > eq_${fname}_submit.sh
#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=4:00:00
#$ -pe mvapich12 12
#$ -q sandy.q

. /etc/profile.d/modules.sh

module load mbm
module delete gcc
module load gromacs/5.0.5/fftw/3.3.3/gcc/4.9.2/openmpi/1.8.5/magny

# Run command
fname="${fname}"
# Equlibration at 310K
mpirun mdrun_mpi -deffnm "${fname}_eq" -maxh 24;

### GO TO MD RUN DIRECTORY ####

cd ../${fname}_md/;
echo "Run Gromacs commands";
editconf -f ../${fname}_eq/${fname}_eq.gro -o ${fname}_eq.pdb;

make_ndx -f ${fname}_eq.pdb -o ${fname}_md.ndx << EOF
q
\EOF

grompp -f NPW_md.mdp -c ${fname}_eq.pdb -p ../${fname}.top -n ${fname}_md.ndx -o ${fname}_md.tpr -t ../${fname}_eq/${fname}_eq.cpt -maxwarn 1 > grompp.log;

echo "GOT TO DIRECTORY"
pwd;
echo "EDIT the md_${fname}_submit.sh File";
echo "Submit the JOB for MD Production Run";
echo "";

EOF

chmod 777 eq_${fname}_submit.sh;

echo "";
echo "GO TO Directory";
pwd;
echo "";
echo "Edit the eq_${fname}_submit.sh File and Submit to cluster";
echo "";


















