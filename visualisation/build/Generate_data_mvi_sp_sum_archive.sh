#!/bin/bash -l

########################################################################################################################
# Generate_data_mvi_sp_sum_archive.sh
#
# ARCHIVED / OBSOLETE — kept for reproducing data generated before the four formal
# fermion topo charge estimators (q_A, q_B, q_Bp, q_C) were derived.
#
# The spectral sum (sp_sum) computed here is:
#   sp_sum(x) = sum_n [ -chi_B^n(x) + 0.5*sign(mu_n)*sqrt(mu_n^2-m_f^2)*rho_n(x) ]
# This formula is NOT one of the four formal estimators defined in topo_charge.pdf
# and was superseded when those estimators were added to Compute_DWF_G5R5.cc.
# sp_sum has also been removed from Compute_DWF_G5R5.cc.
#
# Sections:
#   §A  Spectral sum movie (T-summed, configs animated)
#   §B  IP & corr: spectral sum vs gluonic TCD (FieldDensityFindRegion --compute_PCF)
########################################################################################################################


########################################################################################################################
####################  Input parameters  ###############################################################################
########################################################################################################################

dfiles=()

CDIR=$(pwd)
PDIR=/ccs/home/syamamoto/tmp/src/Grid_cleanedup_for_pullrequest/systems/Frontier/HMC
LDIR=/lustre/orion/phy157/proj-shared/phy157_dwf/syamamoto
HMC=32cube-rho0.124-tau4
HMC_DIR=$PDIR/$HMC
vol=32.32.32.32
REGEN=0   # set to 1 to regenerate; default 0 (safe no-op)

CONF_S=700
CONF_F=709

DATA_DIR=${HMC_DIR}/eigen


########################################################################################################################
### §A  Spectral sum movie (T-summed, configs animated)
# Reads sp_sum_tau_${tau}.${conf} files (written by the old Compute_DWF_G5R5.cc),
# produces a T-summed MPEG and compressed .dat over all configs.
########################################################################################################################

for tau in 0 4; do

    ext=T_summed

    fname=sp_sum_tau_${tau}
    dfile=${HMC}/${fname}_${ext}.dat
    dpath=${HMC_DIR}/${fname}_${ext}.dat
    mpeg=${HMC_DIR}/${fname}_${ext}.avi
    dfiles+=( $dfile )

    if [[ $REGEN == 1 ]] ; then
	F=""
	for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do F+=$DATA_DIR/${conf}/${fname}.${conf},;done
	Fs=${F%?}
	${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --sum T \
	       --mpeg $mpeg --isosurface -0.5 --save_data_to $dpath
    fi

    # evec_tensor is written separately by Compute_DWF_G5R5.cc
    dfile=${HMC}/eigen/evec_tensor_tau_$tau
    dfiles+=( $dfile )

done


########################################################################################################################
### §B  IP & corr: spectral sum vs gluonic TCD
# For each (TD_tau, tau, conf): computes Pearson corr coeff and inner product
# between sp_sum_tau_${tau}.${conf} and Top_dnsty_${TD_tau}_*.${conf} via
# FieldDensityFindRegion --compute_PCF.
# Output: data/corr_ip.dat
# Columns: TD_tau  tau  conf  value
# Rows: (TD_tau=0,4,16) x (tau=0,4) x (Corr, IP) x nconf
########################################################################################################################

dfile=${HMC}/data/corr_ip.dat
dpath=${HMC_DIR}/data/corr_ip.dat
dfiles+=( $dfile )

if [[ $REGEN == 1 ]] ; then

    >$dpath
    for TD_tau in 0 4 16; do
	for tau in 0 4; do
	    >foo
	    for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do
		Fs=${HMC_DIR}/dnsty/Top_dnsty_${TD_tau}_ckpoint_EODWF_lat_smr.${conf},${HMC_DIR}/eigen/${conf}/sp_sum_tau_${tau}.${conf}
		${CDIR}/FieldDensityFindRegion --files ${Fs} --grid $vol --compute_PCF | tee -a foo | cat
	    done
	    paste <(seq -f "%03g" $CONF_S 1 $CONF_F) <(grep "Corr Coeff"   foo) | awk -v tauTC=$TD_tau -v tauQ=$tau '{print tauTC, tauQ, $1, $NF}' >> $dpath
	    paste <(seq -f "%03g" $CONF_S 1 $CONF_F) <(grep "Inner Produc" foo) | awk -v tauTC=$TD_tau -v tauQ=$tau '{print tauTC, tauQ, $1, $NF}' >> $dpath
	done
    done
fi


########################################################################################################################
####################  Archive to Lustre  ##############################################################################
########################################################################################################################

tar -cvf ${LDIR}/sp_sum_archive.tar -C ${PDIR} ${dfiles[@]}


### Refs
# https://zenn.dev/shuh/articles/tar-command-use
