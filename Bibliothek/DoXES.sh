#!/bin/bash

#Arg1: directory of the output of the XAS calculation
#Arg2: index of the excited atom
#Arg3: location of NodePool
TmpDir=`pwd`
OutDir=$1
XATOM=$2
XATOM_SYMBOL=$3
NODE_POOL=$4

cd $OutDir

. ./TMP_INPUT.in${PBS_JOBID}


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


function release {

if [ -z $BibDir ]
then
  ./ReleaseNode.sh "$NODE_POOL" "$OutDir"
else
  ${BibDir}/ReleaseNode.sh "$NODE_POOL" "$OutDir"
fi

}

#Get the number of atomic species
NTYP=`echo "$TMP_ATOMIC_POSITIONS"| awk '{print $1}'| sort | uniq -c |wc -l`
NTYP="$(( ${NTYP} - 1 ))"
#echo $NTYP

#Functions

function error_handler {
  exitstatus=$?
  echo =================================================================
  echo error_handler
  echo problem with an mpirun job
  echo exit status $exitstatus
  if [ $exitstatus -ne 0 ]
  then
    echo exiting
    echo =================================================================
    release
    exit 1
  fi
  echo =================================================================
}

function xes {
# ----------------------------------------------------------------------
local PREFIX=$TMP_MOLNAME
local INPUT=$PREFIX.xes.in
local OUTPUT=$PREFIX.${XATOM_SYMBOL}${XATOM}.xes.out

local SUFFIX=`echo ${TMP_PSEUDO_POT_ES_POST_FIX} | sed 's/\.UPF.*$//'`

echo TMP_EA_HOMO = ${TMP_EA_HOMO}
echo TMP_EA_CORE = ${TMP_EA_CORE}
ECORE=`head -1 ${TMP_MOLNAME}.${XATOM_SYMBOL}${XATOM}.cpot | awk '{ print $2 }'`;
echo ECORE = ${ECORE}

Delta=0;
setdelta=1;
if [ ! -z ${TMP_EA_HOMO} ] ; then 
    Delta=`echo "$Delta + $TMP_EA_HOMO" | bc -l`;
else
#    echo "Warning TMP_EA_HOMO is empty!" ; 
    setdelta=0;
fi
if [ ! -z ${TMP_EA_CORE} ] ; then
    Delta=`echo "$Delta - $TMP_EA_CORE" | bc -l`;
else
#    echo "Warning TMP_EA_CORE is empty!" ;  
    setdelta=0;
fi
if [ ! -z ${ECORE} ] ; then
    Delta=`echo " -1.0*$ECORE - $Delta " | bc -l`;
else
#    echo "Warning ECORE is empty!" ;     
    setdelta=0;
fi
if [ $setdelta -gt 0 ] ; then
    Delta=`echo "13.60569*$Delta" | bc -l`;
    echo DeltaE = $Delta ;
else
    Delta=0;
    echo "Warning Delta is zero";
fi
cat > $INPUT <<EOFXES
&inputpp
  prefix='${TMP_MOLNAME}',
  outdir='$OutDir',
  Nener=5000
  Emin=-40
  Emax=10
  broadening=${TMP_SIGMA}
  DeltaE=${Delta}
  chapprox='GS'
  readcorerep=.true.,   
  spectype='XES'
/
COREREP
  $NTYP $TMP_NAT
  1
  $NTYP  $XATOM  '$TMP_PSEUDO_DIR/${XATOM_SYMBOL}.${SUFFIX}.pos'
EOFXES

#Calcualting of the Delta
local XAS="$TMP_PARA_PREFIX $TMP_BIN_DIR/pw2xas.x  $TMP_SCF_POSTFIX"
$XAS < $INPUT > $OUTPUT || error_handler

}

function core_eig {
local INPUT=${TMP_MOLNAME}.pp.in
local FILEOUT=${TMP_MOLNAME}.${XATOM_SYMBOL}${XATOM}.cpot
local OUTPUT=${TMP_MOLNAME}.${XATOM_SYMBOL}${XATOM}.pp.out

NLIN=`echo "$XATOM + 1" | bc` ;
#echo "grep -A${TMP_NAT} 'positions (alat units)' ${TMP_MOLNAME}.scf.out | tail -n$NLIN | tail -1";
xtmp=`grep -A${TMP_NAT} 'positions (alat units)' ${TMP_MOLNAME}.scf.out | head -n$NLIN | tail -1 | awk '{print $7,":",$8,":",$9}'`;
IFS=':' read -a x0 <<< "$xtmp";
#x0=(${xtmp//:/ });
#echo ${x0[0]},${x0[1]},${x0[2]};#[0],$x0[1],$x0[2];

cat > $INPUT <<EOFPP
&INPUTPP
prefix='${TMP_MOLNAME}',
outdir='./',
filplot='vnuc.tmp',
plot_num=1,
/
EOFPP

local PPRUN="$TMP_PARA_PREFIX $TMP_BIN_DIR/pp.x  $TMP_SCF_POSTFIX"

$PPRUN < $INPUT > $OUTPUT || error_handler

local INPUT=${TMP_MOLNAME}.pplot.in

cat > $INPUT <<EOFPP
&INPUTPP
/

&PLOT
nfile=1,
filepp(1)='vnuc.tmp',
weight(1)=1,
iflag=1,
output_format=0,
fileout='${FILEOUT}',
e1(1)=1
e1(2)=0
e1(3)=0
x0(1)=${x0[0]}
x0(2)=${x0[1]}
x0(3)=${x0[2]}
nx=101
/
EOFPP

local PPRUN="$TMP_PARA_PREFIX $TMP_BIN_DIR/pp.x"

$PPRUN < $INPUT > $OUTPUT || error_handler

}


# Executable section
#Check what actually needs to be done and define JOB
PPOUT=${TMP_MOLNAME}.${XATOM_SYMBOL}${XATOM}.pp.out
if ! ( [ -f $PPOUT ] && grep -q 'Output format' $PPOUT ) ; then
echo "$XATOM_SYMBOL$XATOM : Starting PP..."
core_eig
#cpotev=`echo "$cpot*13.60569" | bc -l`;
#echo $cpotev ;
fi

XESOUT=$TMP_MOLNAME.${XATOM_SYMBOL}${XATOM}.xes.out

if ! ( [ -f $XESOUT ] && grep -q 'done spec' $XESOUT ) ; then
echo "$XATOM_SYMBOL$XATOM : Starting XES..."
xes $TMP_XAS_ARG
if [ -f $TMP_MOLNAME.pw.xes ] ; then
    mv $TMP_MOLNAME.pw.xes $TMP_MOLNAME.${XATOM_SYMBOL}${XATOM}.xes.dat
fi
echo "done"

fi

release


exit
