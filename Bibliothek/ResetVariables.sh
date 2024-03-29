#!/bin/bash

#Script makes a copy of the inputblock.in file 
#Arg1 (optional): defines the directory, where the TMP_INNPUT.in will be saved

function F_ResetVariables {

if [[ $# -eq 0 ]] ; then

local CallingDir=`pwd`

else

local CallingDir=$1

fi 

. "$CallingDir/Input_Block.in"

cat > "${CallingDir}/TMP_INPUT.in${PBS_JOBID}"<<EOF
TMP_PREFIX="$PREFIX"
TMP_BIN_DIR="$BIN_DIR"
TMP_PSEUDO_DIR="$PSEUDO_DIR"
TMP_TMP_DIR="$TMP_DIR"
TMP_PARA_PREFIX="$PARA_PREFIX"
TMP_PARA_POSTFIX="$PARA_POSTFIX"
TMP_PW_POSTFIX="$PW_POSTFIX"
TMP_SCF_POSTFIX="$SCF_POSTFIX"
TMP_NSCF_POSTFIX="$NSCF_POSTFIX"
TMP_MOLNAME="$MOLNAME"
TMP_BIN_LIST="$BIN_LIST"	
TMP_PSEUDO_LIST="$PSEUDO_LIST"
TMP_PSEUDO_POT_ES_POST_FIX="$PSEUDO_POT_ES_POST_FIX"
TMP_DT=$DT
TMP_NSTEP=$NSTEP
TMP_CHAPPROX="$CHAPPROX"
TMP_N="$N"
TMP_XCH_STEP=$XCH_STEP
TMP_A=$A
TMP_B=$B
TMP_C=$C
TMP_COSAB=$COSAB
TMP_COSAC=$COSAC
TMP_COSBC=$COSBC
TMP_CELLDM1=$CELLDM1
TMP_CELLDM2=$CELLDM2
TMP_CELLDM3=$CELLDM3
TMP_CELLDM4=$CELLDM4
TMP_CELLDM5=$CELLDM5
TMP_CELLDM6=$CELLDM6
TMP_NAT=$NAT
TMP_NTYP=$NTYP
TMP_IBRAV=$IBRAV
TMP_ECUT_WFC=$ECUT_WFC
TMP_ECUT_RHO=$ECUT_RHO
TMP_NELEC=$NELEC
TMP_tot_charge=$tot_charge
TMP_SMEARING='$SMEARING'
TMP_DEGAUSS='$DEGAUSS'
TMP_DEGAUSS_STEP='$DEGAUSS_STEP'
TMP_DEGAUSS_FAC='$DEGAUSS_FAC'
TMP_DIAG='$DIAG'
TMP_DIAG_NSCF='$DIAG_NSCF'
TMP_ELEC_CONV_THR='$ELEC_CONV_THR'
TMP_ELEC_MIXING_BETA=$ELEC_MIXING_BETA
TMP_mixing_mode='$mixing_mode'
TMP_electron_maxstep=$electron_maxstep
TMP_diago_thr_init=$diago_thr_init
TMP_NBND=$NBND
TMP_NBND_FAC=$NBND_FAC
TMP_ION_TEMPERATURE='$ION_TEMPERATURE'
TMP_POT_EXTRAPOLATION='$POT_EXTRPOLATION'
TMP_WFC_EXTRAPOLATION='$WFC_EXTRPOLATION'
TMP_TEMPW='$TEMPW'
TMP_DELTA_T=$DELTA_T
TMP_NRAISE=$NRAISE
TMP_TOLP='$TOLP'
TMP_RM='from_scratch'
TMP_wf_collect=.TRUE.
TMP_disk_io='high'
TMP_BASIS_TRACE_TOL="$BASIS_TRACE_TOL"
TMP_DELTA=0
TMP_XAS_ARG="$XAS_ARG"
TMP_ATOMIC_SPECIES="`echo "$ATOMIC_SPECIES" | sed '/^ *$/d'`"
TMP_ATOMIC_POSITIONS="`echo "$ATOMIC_POSITIONS" | sed '/^ *$/d'`"
TMP_K_POINTS="`echo "$K_POINTS" | sed '/^ *$/d'`"
TMP_CELL_PARAMETERS="`echo "$CELL_PARAMETERS" | sed '/^ *$/d'`"
TMP_SPIN="`echo "$SPIN" | sed '/^ *$/d'`"
TMP_MAG_SPEC="`echo "$MAG_SPEC" | sed '/^ *$/d'`"
TMP_MAG_VAL="`echo "$MAG_VAL" | sed '/^ *$/d'`"
TMP_TMAG="$TMAG"
TMP_LDAU="`echo "$LDAU" | sed '/^ *$/d'`"
TMP_U_SPEC="`echo "$U_SPEC" | sed '/^ *$/d'`"
TMP_U_VAL="`echo "$U_VAL" | sed '/^ *$/d'`"
TMP_FORCESTRESS="$FORCESTRESS"
TMP_startingpot=$startingpot
VAR2=$VAR1 && [[ -n "$VAR1" ]] && VAR2="'$VAR1'"
TMP_assume_isolated=$assume_isolated && [[ -n "$assume_isolated" ]] && TMP_assume_isolated="'$assume_isolated'"
TMP_input_dft=$input_dft && [[ -n "$input_dft" ]] && TMP_input_dft="'$input_dft'"
TMP_screening_parameter=0.2 && [[ -n "$screening_parameter" ]] && TMP_screening_parameter="$screening_parameter"
TMP_esomo=$set_esomo
TMP_SIGMA=$SIGMA
EOF

}

F_ResetVariables $1

exit

