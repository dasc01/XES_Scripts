#!/bin/bash

TmpDir=`pwd`
OutDir=$1

cd $OutDir

. ./TMP_INPUT.in${PBS_JOBID}

NODE_POOL="${TmpDir}/Nodes"
MyDir=`echo $0 | sed -e 's/SCF.sh//g'`

#Arg1: directory of the output of the SCF  calculation
#Arg2: Indicator whether run should be parallel (1= parallel; 0= not parallel)

#get hostfile information
XAS_HOSTFILE=`ls | grep "Node-${PBS_JOBID}" | sed -n '1p'`
PPN=`sed 's/,/ /g' "$XAS_HOSTFILE" | wc -w`
if [ -z $CRAY_SITE_LIST_DIR ]; then
  XAS_HOSTS=`cat "${XAS_HOSTFILE}"`
  TMP_PARA_PREFIX="$TMP_PARA_PREFIX -host ${XAS_HOSTS} -np $PPN"
else
  mppnppn=`qstat -f $PBS_JOBID | grep Resource_List.mppnppn | awk '{print $3}'`
  TMP_PARA_PREFIX="$TMP_PARA_PREFIX -n $PPN -N $mppnppn"
fi

function F_SCF {
local NAME=$TMP_MOLNAME
local PREFIX="$NAME.scf"
local INPUT="$OutDir/$PREFIX.in"
local OUTPUT="$OutDir/$PREFIX.out"

TMP_PW_POSTFIX=`echo TMP_PW_POSTFIX | sed -e "s/.*npool.*//"`
# Check if parallel run
if [[ $2 -eq 1 ]] ; then	#parallel

local SCF_HOSTFILE=`ls|grep "Node-${PBS_JOBID}" |sed -n '1p' `
local PPN=`sed 's/,/ /g' "$SCF_HOSTFILE" | wc -w`


local SCF_HOSTS=`cat "${SCF_HOSTFILE}"` 
local PW_COMMAND="$TMP_PARA_PREFIX $TMP_BIN_DIR/pw.x $TMP_PARA_POSTFIX $TMP_PW_POSTFIX $TMP_SCF_POSTFIX"

else  #not parallel

local PW_COMMAND="$TMP_PARA_PREFIX $TMP_BIN_DIR/pw.x $TMP_PARA_POSTFIX $TMP_PW_POSTFIX $TMP_SCF_POSTFIX"

fi 

if [[ $3 -eq 1 ]]; then
  startingpot="startingpot='file'"
else
  startingpot=""
fi

#Check if XCH or FCH
local tot_charge=$TMP_tot_charge
if [[ -z $TMP_tot_charge ]]; then
  tot_charge=0.0
fi
if test "$TMP_CHAPPROX" == "FCH" ; then
    tot_charge=`echo "$tot_charge+1.0" | bc -l`
fi
echo tot_charge = $tot_charge
if [ -z $TMP_electron_maxstep ]
then
  TMP_electron_maxstep=30
fi

#Get the number of atomic species
echo "$TMP_ATOMIC_POSITIONS"
local NTYP=`echo "$TMP_ATOMIC_POSITIONS"| awk '{if(NR>1){print $1}}'| sort | uniq |wc -l`

local LATTICE="ibrav=$TMP_IBRAV"
if [[ -n $TMP_A ]] ; then LATTICE="${LATTICE} , a=$TMP_A"; fi
if [[ -n $TMP_B ]] ; then LATTICE="${LATTICE} , b=$TMP_B"; fi
if [[ -n $TMP_C ]] ; then LATTICE="${LATTICE} , c=$TMP_C"; fi
if [[ -n $TMP_COSAB ]] ; then LATTICE="${LATTICE} , cosab=$TMP_COSAB"; fi
if [[ -n $TMP_COSAC ]] ; then LATTICE="${LATTICE} , cosac=$TMP_COSAC"; fi
if [[ -n $TMP_COSBC ]] ; then LATTICE="${LATTICE} , cosbc=$TMP_COSBC"; fi
if [[ -n $TMP_CELLDM1 ]] ; then LATTICE="${LATTICE} , celldm(1)=$TMP_CELLDM1"; fi
if [[ -n $TMP_CELLDM2 ]] ; then LATTICE="${LATTICE} , celldm(2)=$TMP_CELLDM2"; fi
if [[ -n $TMP_CELLDM3 ]] ; then LATTICE="${LATTICE} , celldm(3)=$TMP_CELLDM3"; fi
if [[ -n $TMP_CELLDM4 ]] ; then LATTICE="${LATTICE} , celldm(4)=$TMP_CELLDM4"; fi
if [[ -n $TMP_CELLDM5 ]] ; then LATTICE="${LATTICE} , celldm(5)=$TMP_CELLDM5"; fi
if [[ -n $TMP_CELLDM6 ]] ; then LATTICE="${LATTICE} , celldm(6)=$TMP_CELLDM6"; fi

cat > $INPUT <<EOFSCF
&control
    calculation='scf'
    prefix='$NAME'
    pseudo_dir='$TMP_PSEUDO_DIR'
    outdir='$OutDir'
    wf_collect=$TMP_wf_collect
    disk_io='$TMP_disk_io'
    restart_mode='$TMP_RM'
    $TMP_FORCESTRESS
/ 
&system
    $LATTICE
    nat=$TMP_NAT, ntyp=$NTYP, tot_charge=$tot_charge
    nbnd=$TMP_NBND, occupations='smearing', degauss=$TMP_DEGAUSS
    ecutwfc=$TMP_ECUT_WFC, ecutrho=$TMP_ECUT_RHO, assume_isolated=$TMP_assume_isolated, input_dft=$TMP_input_dft, screening_parameter=$TMP_screening_parameter
    $TMP_SPIN 
    $TMP_LDAU 
/
&electrons
    electron_maxstep=$TMP_electron_maxstep
    conv_thr=$TMP_ELEC_CONV_THR
    mixing_beta=$TMP_ELEC_MIXING_BETA
    $startingpot
/ 
$TMP_ATOMIC_SPECIES
$TMP_ATOMIC_POSITIONS
$TMP_K_POINTS
$TMP_CELL_PARAMETERS
EOFSCF



$PW_COMMAND < $INPUT > $OUTPUT

}

F_SCF $1 $2

count=0
while ! grep -q 'convergence has been achieved' $1/$TMP_MOLNAME.scf.out
do
  count=$((count + 1))
  echo restart SCF $count
  cp $1/$TMP_MOLNAME.scf.out $1/$TMP_MOLNAME.scf.out-$count
  # restart
  F_SCF $1 $2 1
  if [[ $count -gt 2 ]]; then break; fi
done

if [[ $2 -eq 1 ]] ; then	# If SCF parallel
${MyDir}ReleaseNode.sh "$NODE_POOL" "$1"

fi 

exit
