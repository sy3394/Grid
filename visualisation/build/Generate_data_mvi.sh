#!/bin/bash -l

dfiles=()

##### INPUT ####################
CDIR=$(pwd)
PDIR=/ccs/home/syamamoto/tmp/src/Grid_cleanedup_for_pullrequest/systems/Frontier/HMC
LDIR=/lustre/orion/phy157/proj-shared/phy157_dwf/syamamoto
HMC=32cube-rho0.124-tau4
HMC_DIR=$PDIR/$HMC
vol=32.32.32.32
REGEN=1
###############################


########################################################################################################################
####################      Analyze Dnsty & Evec Data for Each Config      ###############################################
########################################################################################################################
# Chiral Matrix
# E density & Topo charge
# evec dnsty
# chiral dnsty                                     <- not necessary?
# spectral reconstruction of TC dnsty
# Measure Similarity of Sp reconst & topo charge
# NOTE:
#   For now) data file format: (slow -> fast) : config, X, Y, Z where T summed over


########    INPUT: DDIR dependent   #########
CONF_S=700
CONF_F=709
nconv=15  #relevant for Evec
#############################################

### Chiral Matrix

for tau in 0 4; do
    dfile=${HMC}/eigen/chiral_matrix_real_tau_$tau
    >${PDIR}/$dfile

    for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do
	awk -v Nk=10 '{if(NR<Nk+1){for(i=1;i<Nk+1;i++) printf("%f ",$i)};printf("\n");}' ${HMC_DIR}/eigen/$conf/chiral_matrix_real_tau_${tau}_$conf >> ${PDIR}/$dfile
    done

dfiles+=( $dfile )
done

### E density & Topo charge

DATA_DIR=${HMC_DIR}/dnsty

for tau in 0 4; do
    for D_TYPE in E_dnsty Top_dnsty; do

	ext=T_summed
	
	fname=${D_TYPE}_${tau}_ckpoint_EODWF_lat_smr
	dpath=${HMC_DIR}/${fname}_${ext}.dat
	dfile=${HMC}/${fname}_${ext}.dat
	mpeg=${HMC_DIR}/${fname}_${ext}.avi
	dfiles+=( $dfile )
	
	if [[ $REGEN == 0 ]] ; then
	    F=""
	    for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do F+=$DATA_DIR/${fname}.${conf},;done
	    Fs=${F%?}
	    seq -f "%03g" $((CONF_S)) 1 $CONF_F > foo_ind
	    
	    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --omit_intcpts -1 \
		   --mpeg $mpeg --isosurface -0.75 --save_data_to $dpath --index_file foo_ind
	fi
    
    done
done

### Evec
# Assume: nconv=15

DATA_DIR=${HMC_DIR}/eigen

for tau in 0 4; do
    for i_evec in `seq 0 1 $nconv`; do
	
	ext=TLs_summed
	
	fname=evec_density_${i_evec}_tau_$tau
	dpath=${HMC_DIR}/${fname}_${ext}.dat
	dfile=${HMC}/${fname}_${ext}.dat
	mpeg=${HMC_DIR}/${fname}_${ext}.avi
	dfiles+=( $dfile )
	
	if [[ $REGEN == 0 ]] ; then
	    
	    F=""
	    for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do F+=$DATA_DIR/${conf}/${fname}.${conf},;done
	    Fs=${F%?}
	    seq -f "%03g" $((CONF_S)) 1 $CONF_F > foo_ind

	    sf=`awk -v t=$tau -v n=$i_evec -v N=$nconv 'BEGIN{print -1*( 0.3 + 0.4*t/4 + 0.9*n/N)}'`
	    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --dynm_dir 5 --omit_dirs 0.4 --omit_intcpts -1.-1 \
		   --mpeg $mpeg --isosurface $sf --save_data_to $dpath --index_file foo_ind
	fi
	
    done
done

### Chiral

for tau in 0 4; do
    for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do

	ext=ZT_summed
	
	fname=chiral_density_tau_${tau}_
	dfile=${HMC}/${fname}${ext}_${ext}_${conf}.dat
	dpath=${HMC_DIR}/${fname}_${ext}_${conf}.dat
	mpeg=${HMC_DIR}/${fname}_${ext}_${conf}.avi
	dfiles+=( $dfile )

	if [[ $REGEN == 0 ]] ; then
	    F=""
	    for f in `ls $DATA_DIR/${conf}/${fname}_* | grep tau_${tau} | sort -n -t _ -k3`; do F+=$f,; done
	    Fs=${F%?}
	    seq -f "%03g" $((CONF_S)) 1 $CONF_F > foo_ind
	    
	    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --dynm_dir 5 --omit_dirs 3.4 --omit_intcpts -1.-1 \
		   --mpeg $mpeg --isosurface -0.7 --save_data_to $dpath --index_file foo_ind
	fi

    done
done

### Spectral sum

for tau in 0 4; do

    ext=T_summed
    
    fname=sp_sum_tau_${tau}
    dfile=${HMC}/${fname}_${ext}.dat
    dpath=${HMC_DIR}/${fname}_${ext}.dat
    mpeg=${HMC_DIR}/${fname}_${ext}.avi
    dfiles+=( $dfile )
    
    if [[ $REGEN == 0 ]] ; then
	F=""
	for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do F+=$DATA_DIR/${conf}/${fname}.${conf},;done
	Fs=${F%?}
	${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --omit_intcpts -1 --mpeg $mpeg --isosurface -0.5 --save_data_to $dpath
    fi
    
    dfile=${HMC}/eigen/evec_tensor_tau_$tau
    dfiles+=( $dfile )
    
done


### Measure Similarity of Sp reconst & topo charge

dfile=${HMC}/data/corr_ip.dat
dpath=${HMC_DIR}/data/corr_ip.dat
dfiles+=( $dfile )

if [[ $REGEN == 0 ]] ; then
    
    >$dpath
    for TD_tau in 0 4 16; do
	for tau in 0 4; do
	    >foo
	    for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do
		Fs=${HMC_DIR}/dnsty/Top_dnsty_${TD_tau}_ckpoint_EODWF_lat_smr.${conf},${HMC_DIR}/eigen/${conf}/sp_sum_tau_${tau}.${conf}
		${CDIR}/FieldDensityFindRegion --files ${Fs} --grid $vol --compute_PCF | tee -a foo | cat
	    done
	    paste <(seq -f "%03g" $CONF_S 1 $CONF_F) <(grep "Corr Coeff" foo) | awk -v tauTC=$TD_tau -v tauQ=$tau '{print tauTC, tauQ, $1,$NF}' >> $dpath
	    paste <(seq -f "%03g" $CONF_S 1 $CONF_F) <(grep "Inner Produc" foo) | awk -v tauTC=$TD_tau -v tauQ=$tau '{print tauTC, tauQ, $1,$NF}' >> $dpath
	done
    done
fi

####################################################################################################################
######################    Analyze Evec & Dnsty Data along Trajectory     ###########################################
####################################################################################################################


###########   INPUT   ###########################
conf=702
NCUT=4
#################################################

DATA_DIR=${HMC_DIR}/eigen/${conf}

###################################################################
### H_DWF: Generate time lapse of evec density over a trajectory
###################################################################

for dof in smr lat; do
    for tau in 0 4; do

	#########   Separetely for Each Vector      #############
	
	for n in `seq 0 1 $NCUT`; do
	    	    
	    fname=evec_density_sorted_${n}_tau_${tau}_${dof}

	    ### (Ls, X, Y) => Z,T summed
	    
	    ext=${conf}_ZT_summed

	    mpeg=${HMC_DIR}/${fname}_${ext}.avi
            dpath=${HMC_DIR}/${fname}_${ext}.dat
            dfile=${HMC}/${fname}_${ext}.dat
            mfile=${HMC}/${fname}_${ext}.avi
            dfiles+=( $dfile ) #$mfile )
            echo $tau $dof $n

	    if [[ $REGEN == 1 ]] ; then
		F=""
		for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
		Fs=${F%?}
		#eigen/702/evec_density_sorted_0_tau_0_smr_702_3.345833
		iso="-0.01" #`awk -v t=$tau -v n=$i_evec -v N=$nconv 'BEGIN{print -1*( 0.3 + 0.4*t/4 + 0.9*n/N)}'`
		sep=${conf}_
		${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --dynm_dir 5 --omit_dirs 3.4 --omit_intcpts -1.-1 \
		       --mpeg $mpeg --isosurface $iso -use_fname_as_frame_counter $sep
            fi

	    ### (Ls, X, Y) => Z,T looped <- takes too much time

	    ext=${conf}_ZT_update

            mpeg=${HMC_DIR}/${fname}_${ext}.avi
            dpath=${HMC_DIR}/${fname}_${ext}.dat
            dfile=${HMC}/${fname}_${ext}.dat
            mfile=${HMC}/${fname}_${ext}.avi
            #dfiles+=( $dfile )

            if [[ $REGEN == 0 ]] ; then
                F=""
                for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
                Fs=${F%?}

                iso="-0.01"
                sep=${conf}_
                ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --dynm_dir 5 --omit_dirs 3.4 --xlate_omit_dirs 0.1 --omit_intcpts 0.0 \
                       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
            fi

	    ### (Ls, X, Y) => Z looped, T=23

	    ext=${conf}_Z_update_T23

            mpeg=${HMC_DIR}/${fname}_${ext}.avi
            dpath=${HMC_DIR}/${fname}_${ext}.dat
            dfile=${HMC}/${fname}_${ext}.dat
            mfile=${HMC}/${fname}_${ext}.avi
            #dfiles+=( $dfile )

            if [[ $REGEN == 0 ]] ; then
                F=""
                for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
                Fs=${F%?}

                iso="-0.01"
                sep=${conf}_
                ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --dynm_dir 5 --omit_dirs 3.4 --xlate_omit_dirs 0 --omit_intcpts 0.23 \
                       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
            fi
	    
	    ### (X, Y, Z) + T looped: for comparison with H_W modes

	    ext=T_update
	    
	    mpeg=${HMC_DIR}/${fname}_${ext}.avi
            dpath=${HMC_DIR}/${fname}_${ext}.dat
            dfile=${HMC}/${fname}_${ext}.dat
            mfile=${HMC}/${fname}_${ext}.avi
            #dfiles+=( $dfile ) #$mfile )
	    
	    if [[ $REGEN == 1 ]] ; then
		F=""
                for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
                Fs=${F%?}

		iso="-0.05" #`awk -v t=$tau -v n=$i_evec -v N=$nconv 'BEGIN{print -1*( 0.3 + 0.4*t/4 + 0.9*n/N)}'`
                sep=${conf}_
                ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --dynm_dir 5 --omit_dirs 0.4 --xlate_omit_dirs 1 --omit_intcpts -1.0 \
		       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
	    fi

	done


	#######   All Evecs upto NCUT Summed Over     ########################
	
	#### Sum over Ls and converged evecs for a given tau_W, dof, and t_MD
	    
	if [[ $REGEN == 0 ]] ; then
	    for t in `ls $DATA_DIR/evec_density_sorted_0_tau_${tau}_${dof}_${conf}_*| awk -F _ '{print $NF}' | sort -n`; do
		F=""
		for f in `ls $DATA_DIR/evec_density_sorted_*_tau_${tau}_${dof}_${conf}_${t} | sort -n` ; do F+=$f,; done
		Fs=${F%?}
		
		save_fname=$DATA_DIR/summed_evec_density_sorted_tau_${tau}_${dof}_${conf}_${t}
		${CDIR}/FieldDensityEigen --grid $vol --files2 $Fs --Ls 48 --sum_all_files $save_fname
	    done
	fi

		
	### (X, Y, Z) + T looped: for comparison with fermion force
	
	fname=summed_evec_density_sorted_tau_${tau}_${dof}
	ext=${conf}_T_update
	
	mpeg=${HMC_DIR}/${fname}_${ext}.avi
	dpath=${HMC_DIR}/${fname}_${ext}.dat
	dfile=${HMC}/${fname}_${ext}.dat
	mfile=${HMC}/${fname}_${ext}.avi
	
	if [[ $REGEN == 0 ]] ; then
	    F=""
            for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
            Fs=${F%?}
	    
	    iso="-0.05"
	    sep=${conf}_
	    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --xlate_omit_dirs 0 --omit_intcpts 0 \
		   --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
	fi
    done
done


#######################################################################################################
####  H_W: Evecs from Specflow  #######################################################################
#######################################################################################################

############  INPUT ##################
CONF=702
tau=0
######################################

DATA_DIR=${HMC_DIR}/eigen_Wilson/${CONF}/${CONF}

###### Sum all evec density at each M_5 along the trajectory    #############

if [[ $REGEN == 0 ]] ; then
    for d in `ls -d $DATA_DIR/U_smr_*| sort -n`; do
	for m5 in `ls $d/evec_*_0| awk -F _ '{print $(NF-1)}'`; do
	    F=""
	    for f in `ls $d/evec_${m5}_* | sort -n` ; do F+=$f,; done
	    Fs=${F%?}

	    save_fname=$d/evec_sum_${m5}_${CONF}
	    ${CDIR}/FieldDensityEigen --files1 $Fs --grid $vol --sum_all_files $save_fname
	    
	done
    done
fi

###### Summed evec density along M_5 for the initial config    ############

### (X, Y, Z) + T looped:

t_MD=0.000000
ext=T_update

fname=summed_evec_density_Wilson_tau_${tau}_t_${t_MD}_${CONF}
mpeg=${HMC_DIR}/${fname}_${ext}.avi
dpath=${HMC_DIR}/${fname}_${ext}.dat
dfile=${HMC}/${fname}_${ext}.dat
mfile=${HMC}/${fname}_${ext}.avi

if [[ $REGEN == 0 ]] ; then

    F=""
    for f in `ls $DATA_DIR/U_smr_${t_MD}/evec_sum_* | sort -n` ; do F+=$f,; done
    Fs=${F%?}

    ls $DATA_DIR/U_smr_${t_MD}/evec_sum_* | awk -F _ '{print $NF}' | sort -n > foo_ind
    iso="-0.01"
    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --xlate_omit_dirs 0 --omit_intcpts 0 \
	   --mpeg $mpeg --isosurface $iso --index_file foo_ind
fi

### (X, Y, Z) + T summed: Fix input mass to m=-1.8

m=-1.800000
ext=T_update

fname=summed_evec_density_Wilson_tau_${tau}_m_${m}_${CONF}
mpeg=${HMC_DIR}/${fname}_${ext}.avi
dpath=${HMC_DIR}/${fname}_${ext}.dat
dfile=${HMC}/${fname}_${ext}.dat
mfile=${HMC}/${fname}_${ext}.avi

if [[ $REGEN == 0 ]] ; then
    F=""
    for f in `ls $DATA_DIR/U_smr_*/evec_sum_-1.800000 | sort -n` ; do F+=$f,; done
    Fs=${F%?}

    ls $DATA_DIR/U_smr_*/evec_sum_-1.800000 | awk -F _ '{print $(NF-2)}' | awk -F '/' '{print $1}' | sort -n > foo_ind
    iso="-0.01"
    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --xlate_omit_dirs 0 --omit_intcpts 0 \
	   --mpeg $mpeg --isosurface $iso --index_file foo_ind
fi


###### 0^th Evec density along M_5 for tau_W = 0   ################

#########   INNPUT   #################
CONFS=( 702 7026 70201 70202 70203 70204 70205 703 70301 718 719 )
regens=( 0   0     0     0     0     0     0    0    0    0   0  )
######################################

for((i_conf=0; i_conf<${#CONFS[@]}; i_conf++)); do
    CONF=${CONFS[i_conf]}
    DATA_DIR=${HMC_DIR}/eigen_Wilson/${CONF}/
    
    ### (X, Y, Z) + T looped: 
    
    t_MD=0.000000
    ext=T_update

    fname=evec_density_Wilson_0_tau_${tau}_t_${t_MD}_${CONF}
    mpeg=${HMC_DIR}/${fname}_${ext}.avi
    dpath=${HMC_DIR}/${fname}_${ext}.dat
    dfile=${HMC}/${fname}_${ext}.dat
    mfile=${HMC}/${fname}_${ext}.avi
    #dfiles+=( $dfile )
    
    if [[ $REGEN == ${regens[i_conf]} ]] ; then
    
	F=""
	for f in `ls $DATA_DIR/U_smr_${t_MD}/evec_*_0 | sort -n` ; do F+=$f,; done
	Fs=${F%?}
    
	ls $DATA_DIR/U_smr_${t_MD}/evec_*_0 | awk -F _ '{print $(NF-1)}' | sort -n > foo_ind
	iso="-0.01"
	${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --xlate_omit_dirs 0 --omit_intcpts 0 \
	       --mpeg $mpeg --isosurface $iso --index_file foo_ind
    fi


    ###### 0^th Evec density over a trajectory at m=-1.8   ###############

    ### (X, Y, Z) + T summed: Fix input mass to m=-1.8

    m=-1.800000
    ext=T_summed
    
    fname=evec_density_Wilson_0_tau_${tau}_m_${m}_${CONF}
    mpeg=${HMC_DIR}/${fname}_${ext}.avi
    dpath=${HMC_DIR}/${fname}_${ext}.dat
    dfile=${HMC}/${fname}_${ext}.dat
    mfile=${HMC}/${fname}_${ext}.avi
    #dfiles+=( $dfile ) 
    if [[ $REGEN == ${regens[i_conf]} ]] ; then
    F=""
    for f in `ls $DATA_DIR/U_smr_*/evec_-1.800000_0 | sort -n` ; do F+=$f,; done
    Fs=${F%?}

    ls $DATA_DIR/U_smr_*/evec_-1.800000_0 | awk -F _ '{print $(NF-2)}' | awk -F '/' '{print $1}' | sort -n > foo_ind
    iso="-0.01"
    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --omit_intcpts -1 \
	   --mpeg $mpeg --isosurface $iso --save_data_to $dpath --index_file foo_ind
    fi
    
    ### (X, Y, Z) + T looped: Fix input mass to m=-1.8
    
    tau=0
    m=-1.800000
    ext=T_update
    
    fname=evec_density_Wilson_0_tau_${tau}_m_${m}_${CONF}
    mpeg=${HMC_DIR}/${fname}_${ext}.avi
    dpath=${HMC_DIR}/${fname}_${ext}.dat
    dfile=${HMC}/${fname}_${ext}.dat
    mfile=${HMC}/${fname}_${ext}.avi
    #dfiles+=( $dfile )
    
    if [[ $REGEN == ${regens[i_conf]} ]] ; then
	F=""
	for f in `ls $DATA_DIR/U_smr_*/evec_-1.800000_0 | sort -n` ; do F+=$f,; done
	Fs=${F%?}
	
	ls $DATA_DIR/U_smr_*/evec_-1.800000_0 | awk -F _ '{print $(NF-2)}' | awk -F '/' '{print $1}' | sort -n > foo_ind
	iso="-0.01"
	${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --xlate_omit_dirs 0 --omit_intcpts 0 \
	       --mpeg $mpeg --isosurface $iso --save_data_to $dpath --index_file foo_ind 
    fi
    
    # (X, Y, Z) + T looped: can we track poles of specflow of each mode and visualize the corresp evecs?
    # SKIP FOR NOW
done


#######################################################################################################
###########   E & TC Dnsty along the trajectory
#######################################################################################################


#########   INNPUT   #################
CONFS=( 702 7026 70201 70202 70203 70204 70205 703 70301 718 719 )
regens=( 1   0     0     0     0     0     0    0    0    0   0  )
######################################


##### Topo density over a trajectory (tau_T = 4)  ##########

DATA_DIR=${HMC_DIR}/dnsty

for((i_conf=0; i_conf<${#CONFS[@]}; i_conf++)); do
    conf=${CONFS[i_conf]}
    
    for tau in 0 4; do
	for dof in smr lat; do
	    
	    ext=${conf}_T_update #summed

	    flow_kernel= #_Iwasaki
	    skip=1
	    if [ $conf == 70201 -o $conf == 70202 ]; then
		N_frames=1000
	    elif [ $conf == 719 -o $conf == 703 -o $conf == 702 ] ; then
		N_frames=2000
		skip=4
		ext=${conf}_T_update
	    else
		N_frames=500
	    fi

	    fname=Top_dnsty_${tau}${flow_kernel}_${dof}
    
	    mpeg=${HMC_DIR}/${fname}_${ext}.avi
	    dpath=${HMC_DIR}/${fname}_${ext}.dat
	    dfile=${HMC}/${fname}_${ext}.dat
	    mfile=${HMC}/${fname}_${ext}.avi
	    #dfiles+=( $dfile ) #$mfile )
	    echo $fname $tau $dof
	    
	    if [[ $REGEN == ${regens[i_conf]} ]] ; then
		if [[ "$dof" == "smr" ]] ; then iso0=-0.3; else iso0=-0.4; fi
		if [[ "$tau" == "0" ]] ; then iso=`awk -v a=$iso0 'BEGIN{print a+0.1}'`; else iso=`awk -v a=$iso0 'BEGIN{print a+0.2}'`; fi
		if [ "$dof" == "smr" -a $tau -eq 0 ] ; then iso=-0.4; fi
		if [[ "${ext}" == *"summed"* ]] ; then
		    if [[ "$tau" == "0" ]] ; then
			iso=-0.01;
		    else
			iso=-0.04;
		    fi
		fi
		
		F=""
		for f in `ls $DATA_DIR/${conf}/${fname}.*|sort -t"." -nk2 | awk -v skip=$skip '{if((NR-1)%skip==0)print $0}' | tail -$N_frames`; do F+=$f,;done #tail -375 | head -100
		Fs=${F%?}
		
		if [[ "${ext}" == *"summed"* ]] ; then
		    sep="${dof}."
		    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --sum_omit_dir --use_fname_as_frame_counter $sep \
			   --isosurface $iso --mpeg $mpeg --save_data_to $dpath #$(( (tau+1)*5 ))
		else
		    sep="${dof}."
		    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --xlate_omit_dirs 0 --isosurface $iso \
			   --mpeg $mpeg --use_fname_as_frame_counter $sep
		fi
	    fi
	done
    done
done


#####  Force densities over a trajectory (702 -> 703: tau_T = 4)  ############

#########   INNPUT   #################
CONFS=( 702 7026 70201 70202 70203 70204 70205 703 70301 718 719 )
regens=( 1   0     0     0     0     0     0    0    0    0   0  )
######################################

DATA_DIR=${HMC_DIR}/snapshots

for((i_conf=0; i_conf<${#CONFS[@]}; i_conf++)); do
    conf=${CONFS[i_conf]}
    
    for force in ExactOneFlavourRatioPseudoFermionAction TwoFlavourEvenOddRatioPseudoFermionActiondet_0.5_det_1 TwoFlavourEvenOddRatioPseudoFermionActiondet_0.25_det_0.5 TwoFlavourEvenOddRatioPseudoFermionActiondet_0.1_det_0.25 TwoFlavourEvenOddRatioPseudoFermionActiondet_0.05_det_0.1 TwoFlavourEvenOddRatioPseudoFermionActiondet_0.0047_det_0.05; do

	for dof in smr lat; do

	    ext=${conf}_T_update #summed
	    fname=F_${force}_${dof}
	

	    mpeg=${HMC_DIR}/${fname}_${ext}.avi
	    dpath=${HMC_DIR}/${fname}_${ext}.dat
	    dfile=${HMC}/${fname}_${ext}.dat
	    mfile=${HMC}/${fname}_${ext}.avi
	    dfiles+=( $dfile ) #$mfile )
	
	    if [[ $REGEN == ${regens[i_conf]} ]] ; then

		if [[ "$dof" == "smr" ]] ; then iso=-0.4; else iso=-0.45; fi
		if [[ $force != *"Two"* ]] ; then iso=-0.2; else iso=-0.3; fi
		if [ "$force" == "TwoFlavourEvenOddRatioPseudoFermionActiondet_0.0047_det_0.05" ] ; then iso=-0.03; fi
	    
		if [ $conf == 70201 -o $conf == 70202 ]; then
                    N_frames=1000
		elif [ $conf == 719 -o $conf == 703 ] ; then
                    N_frames=2000
		else
                    N_frames=70
		fi

	    
		F=""
		for f in `ls $DATA_DIR/${conf}/${fname}.*|awk -F . '{print $NF, $0}' | sort  -nk1| cut -f2- -d' ' | tail -$N_frames`; do F+=$f,;done
		Fs=${F%?}
		tail -$N_frames traj_times_top_${conf} > foo_ind
	    
		${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --xlate_omit_dirs 0 --isosurface $iso \
		       --mpeg $mpeg --save_data_to $dpath --index_file foo_ind
	    fi
	done
    done
done

for force in IwasakiGaugeAction JacobianAction; do

    if [[ "$force" == "IwasakiGaugeAction" ]] ; then iters="smr lat"; else iters=lat; fi
    for dof in $iters; do

	fname=F_${force}_${dof}
	ext=T_update #summed
		
	mpeg=${HMC_DIR}/${fname}_${ext}.avi
	dpath=${HMC_DIR}/${fname}_${ext}.dat
	dfile=${HMC}/${fname}_${ext}.dat
	mfile=${HMC}/${fname}_${ext}.avi
	if [[ "$ext" == *"update"* ]] ; then
	    dfiles+=( $mfile )
	else
	    dfiles+=( $dfile $mfile )
	fi

	if [[ "$dof" == "smr" ]] ; then iso=-0.45; else iso=-0.58; fi
	if [[ "$force" == "JacobianAction" ]] ; then iso=-0.61; fi
	
	if [[ $REGEN == 0 ]] ; then
            F=""
            for f in `ls $DATA_DIR/${conf}/${fname}.*|awk -F . '{print $NF, $0}' | sort  -nk1| cut -f2- -d' ' | tail -375 |head -100`; do F+=$f,;done
            Fs=${F%?}
	    tail -376 traj_times|head -100 > foo_ind #traj_times includes 4.00, at whcich force is not computed
	    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --dynm_dir 4 --omit_dirs 3 --xlate_omit_dirs 0 --isosurface $iso \
		   --mpeg $mpeg --save_data_to $dpath --index_file foo_ind
	fi
    done
done


###  Compute zero-modes filtering of gluonic TCD   ##############  moved to somehere else

DATA_DIR1=${HMC_DIR}/dnsty
DATA_DIR2=${HMC_DIR}/eigen

dfile=${HMC}/data/filter_TCD.dat
dpath=${HMC_DIR}/data/filter_TCD.dat
dfiles+=( $dfile )

>$dpath
if [[ $REGEN == 0 ]] ; then
    for TD_tau in 0 4; do
	fname2=evec_density_0_tau_${TD_tau}
	F1s=""
	for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do  F1s+=$DATA_DIR1/Top_dnsty_${TD_tau}_ckpoint_EODWF_lat_smr.${conf},;done
	F1s=${F1s%?}
	F2s=""
	for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do  F2s+=$DATA_DIR2/$conf/$fname2.${conf},;done
	F2s=${F2s%?}
	
	${CDIR}/FieldDensityEigen --files1 $F1s --files2 $F2s --grid $vol --Ls 48 --compareTCD_defs --cut 0.00001 | tee -a foo2 | cat
	grep "Filtered Sum" foo2 | awk -v tau=$TD_tau '{print tau, $3, $4, $8}' >> $dpath
    done
fi
wait

tar -cvf ${LDIR}/eigen_data.tar -C ${PDIR} ${dfiles[@]}


### Refs
# https://zenn.dev/shuh/articles/tar-command-use
# https://stackoverflow.com/questions/50338201/how-to-compress-and-tar-a-folder-in-linux
# https://askubuntu.com/questions/392885/how-can-i-view-the-contents-of-tar-gz-file-without-extracting-from-the-command-l
