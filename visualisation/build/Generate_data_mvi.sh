#!/bin/bash -l

########################################################################################################################
# Generate_data_mvi.sh
#
# Orchestration script run on Frontier (ORNL) to process eigenvector and topological
# charge density data produced by the DWF HMC simulation, generate VTK movies (.avi),
# write compressed data files (.dat), and archive everything to Lustre.
#
# INDEX
# =====
#  §1   Configuration & global input parameters
#  §2   Ensemble averages (configs CONF_S .. CONF_F)
#       §2.1  Chiral matrix
#       §2.2  Gluonic E density & topo charge (T-summed, configs as animation axis)
#       §2.3  H_DWF eigenvector density (TLs-summed, configs as animation axis)
#       §2.4  Chiral density (ZT-summed, configs as animation axis)
#       §2.5  Fermion topo charge density: compute q_A/B/C, IP & corr vs gluonic TCD,
#             and 3-panel T-animated movie of all 3 definitions
#  §3   Trajectory analysis — H_DWF evecs, config 702
#       §3.1  Per-evec density movies: ZT-summed, ZT-updated, T-updated (5D & 4D)
#       §3.2  Sum over all modes at each tau_MD snapshot (FieldDensityEigen)
#       §3.3  Summed-evec density movie (T-updated, 4D)
#  §4   Trajectory analysis — H_W (Wilson) evecs via specflow for config 702
#       §4.1  Sum evec density at each M_5 step (initial config t_MD=0)
#       §4.2  Summed evec density movies (T-updated) for t_MD=0 and m=-1.8
#       §4.3  0th evec density over trajectory at m=-1.8 (multi-config loop)
#  §5   Gluonic E & TC density along trajectory (multi-config loop)
#  §6   Fermion force density along trajectory (multi-config loop)
#  §7   Gauge action force density (Iwasaki, Jacobian)
#  §8   Gluonic TCD filtering via fermion zero modes (FieldDensityEigen --compareTCD_defs)
#  §9   Archive to Lustre
########################################################################################################################


########################################################################################################################
####################  §1  Configuration & global input parameters  ####################################################
########################################################################################################################

dfiles=()

CDIR=$(pwd)
PDIR=/ccs/home/syamamoto/tmp/src/Grid_cleanedup_for_pullrequest/systems/Frontier/HMC
LDIR=/lustre/orion/phy157/proj-shared/phy157_dwf/syamamoto
HMC=32cube-rho0.124-tau4
HMC_DIR=$PDIR/$HMC
vol=32.32.32.32
REGEN=1


########################################################################################################################
####################  §2  Ensemble averages (configs CONF_S .. CONF_F)  ##############################################
########################################################################################################################

########    INPUT: DDIR dependent   #########
CONF_S=700
CONF_F=709
nconv=15  # number of converged H_DWF eigenvectors
#############################################


########################################################################################################################
### §2.1  Chiral matrix
# Reads chiral_matrix_real_tau_* per config, concatenates into one file per tau.
########################################################################################################################

for tau in 0 4; do
    dfile=${HMC}/eigen/chiral_matrix_real_tau_$tau
    >${PDIR}/$dfile

    for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do
	awk -v Nk=10 '{if(NR<Nk+1){for(i=1;i<Nk+1;i++) printf("%f ",$i)};printf("\n");}' ${HMC_DIR}/eigen/$conf/chiral_matrix_real_tau_${tau}_$conf >> ${PDIR}/$dfile
    done

dfiles+=( $dfile )
done


########################################################################################################################
### §2.2  Gluonic E density & topo charge (T-summed, configs animated)
# Animates over configs; T is summed; saves compressed 4D (X,Y,Z) data and MPEG.
########################################################################################################################

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

	    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --sum T \
		   --mpeg $mpeg --isosurface -0.75 --save_data_to $dpath --index_file foo_ind
	fi

    done
done


########################################################################################################################
### §2.3  H_DWF eigenvector density (TLs-summed, configs animated)
# 5D input; sums over Ls and T to display 3D (X,Y,Z) density per mode.
########################################################################################################################

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
	    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --animate configs --sum Ls --sum T \
		   --mpeg $mpeg --isosurface $sf --save_data_to $dpath --index_file foo_ind
	fi

    done
done


########################################################################################################################
### §2.4  Chiral density (ZT-summed, per-config movies)
# 5D input; sums over Z and T; per-config loop (Ls modes as animation axis).
########################################################################################################################

for tau in 0 4; do
    for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do

	ext=ZT_summed

	# Note: fname ends with trailing '_'; do NOT insert extra '_' before ext
	fname=chiral_density_tau_${tau}_
	dfile=${HMC}/${fname}${ext}_${conf}.dat
	dpath=${HMC_DIR}/${fname}${ext}_${conf}.dat
	mpeg=${HMC_DIR}/${fname}${ext}_${conf}.avi
	dfiles+=( $dfile )

	if [[ $REGEN == 0 ]] ; then
	    F=""
	    for f in `ls $DATA_DIR/${conf}/${fname}_* | grep tau_${tau} | sort -n -t _ -k3`; do F+=$f,; done
	    Fs=${F%?}
	    seq -f "%03g" $((CONF_S)) 1 $CONF_F > foo_ind

	    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --animate configs --sum Z --sum T \
		   --mpeg $mpeg --isosurface -0.7 --save_data_to $dpath --index_file foo_ind
	fi

    done
done

########################################################################################################################
### §2.5  Fermion topo charge density: q_A/B/C per config, IP & corr vs gluonic TCD,
###       and 3-panel T-animated movie of all 3 definitions
#
# For each (conf, tau_W):
#   Step 1 + 2 (single FieldDensityEigen call):
#     --files2 evec density files (5D, Ncut modes)  --Ls 48  --evals eigenvalue file
#     --topo_out  writes topo_q_{A,B,Bp,C}_<tau>_<dof>.<conf>  (no extension)
#     --files1 gluonic TCD at TD_tau=0,4,16  --topo_compare  --conf_id <conf>
#       => prints "# q_X" label lines + "Topo PCF Corr/IP: TD_tau conf value" to stdout
#       => shell splits by label into 4 per-def output files (pure numeric, notebook-ready)
#   Step 3 (FieldDensityAnimateMultiFiles):
#     --animate T  with 3 input files => 3-panel side-by-side T-animated movie
#
# Output data files (one per def, format: TD_tau  conf  corr  — matching notebook schema):
#   data/corr_ip_q_A.dat   data/corr_ip_q_B.dat   data/corr_ip_q_C.dat
#   Rows: (TD_tau=0,4,16) x (tau=0,4) x (Corr, IP) per conf  => 12*nconf rows total
#
# NOTE: EVALS_FILE is required for m_gap/mu_n weighting in all estimators q_A/B/C.
#       Without it, m_gap falls back to mass_f (fermion mass parameter).
#
# --comp_file (optional): qlat stochastic DWF TCD reference fields (SCIDAC).
#   Located in visualisation/data/topo_field_<idx>.scidac (converted by pickle_to_scidac.py).
#   If present, each is compared against q_A/B/C; results go to data/comp_ref_<def>.dat.
#   Format: comp_idx  conf  Q_evec  Q_ref  Corr  IP  rms_diff
########################################################################################################################

##########   INPUT   ######################################
CONFS=(    700  701  702  703  704  705  706  707  708  709 )
REGENS_25=(  0    0    1    0    0    0    0    0    0    0 )  # 1=run, 0=skip
# 702: has comp_file data (topo_field_*.scidac) and correctly ordered evec files.
# Add more confs here as data becomes available.
###########################################################

DATA_DIR_dnsty=${HMC_DIR}/dnsty
COMP_DIR=$(cd "$(dirname "$0")/.." && pwd)/data   # qlat reference SCIDAC files: visualisation/data/topo_field_<idx>.scidac

# Output IP/corr files — pure numeric, no string columns
# WEIGHTS: comma-separated list of weight tokens passed to --weights.
#   Each token is "sign", "mgap", or "label=value" (e.g. mext=0.011).
#   Default: "sign,mgap"  (both tracks always).
#   Example: WEIGHTS=sign,mgap,mext=0.011 ./Generate_data_mvi.sh
WEIGHTS="${WEIGHTS:-sign,mgap}"

# Parse WEIGHTS -> _weight_labels array (label part before any '=')
_weight_labels=()
IFS=',' read -ra _wtoks <<< "$WEIGHTS"
for tok in "${_wtoks[@]}"; do
    _weight_labels+=("${tok%%=*}")
done

# Build _q_defs_all from {A,B,C} x _weight_labels
_q_defs_all=""
for _lbl in "${_weight_labels[@]}"; do
    for _def in A B C; do
        _q_defs_all+="q_${_def}_${_lbl} "
    done
done
_q_defs_all="${_q_defs_all% }"   # trim trailing space

for q_def in $_q_defs_all; do
    >${HMC_DIR}/data/corr_ip_${q_def}.dat
    dfiles+=( ${HMC}/data/corr_ip_${q_def}.dat )
done

# Output comp-ref files — appended each run; to reset: rm ${HMC_DIR}/data/comp_ref_*.dat
# Format: comp_idx  tau  conf  Q_evec  Q_ref  Corr  IP  rms_diff
for q_def in $_q_defs_all; do
    touch ${HMC_DIR}/data/comp_ref_${q_def}.dat
done

# Output stochastic-vs-gluonic file — appended each run
# Format: comp_idx  TD_tau  tau  conf  Q_ref  Q_gluon  Corr  IP
touch ${HMC_DIR}/data/corr_ip_stoch.dat
dfiles+=( ${HMC}/data/corr_ip_stoch.dat )

for i_conf in "${!CONFS[@]}"; do
    conf=${CONFS[$i_conf]}
    REGEN=${REGENS_25[$i_conf]}
    if [[ $REGEN == 1 ]] ; then

        DATA_DIR_eigen=${HMC_DIR}/eigen/${conf}
        DATA_DIR_topo=${HMC_DIR}/eigen/${conf}   # Top_dnsty_q_X_<tau>_smr.<conf> lives here

        for tau in 0 4; do

            ### Gather eigenvector density files (n = 0 .. nconv-1)
            F2=""
            for n in $(seq 0 1 $nconv); do
                f=${DATA_DIR_eigen}/evec_density_${n}_tau_${tau}.${conf}
                [[ -f $f ]] && F2+=$f,
            done
            F2s=${F2%?}
            [[ -z "$F2s" ]] && { echo "No evec files for conf=$conf tau=$tau, skipping"; continue; }

            ### Gluonic TCD files at TD_tau = 0, 4, 16  (smr hardcoded: actual Frontier naming)
            F1=""
            for TD_tau in 0 4 16; do
                f=${DATA_DIR_dnsty}/Top_dnsty_${TD_tau}_ckpoint_EODWF_lat_smr.${conf}
                [[ -f $f ]] && F1+=$f, || echo "Warning: gluonic TCD not found: $f"
            done
            F1s=${F1%?}

            ### Eigenvalue file (enables m_gap/mu_n weighting; written by Compute_DWF_G5R5)
            EVALS_FILE=${DATA_DIR_eigen}/eigenvalues_tau_${tau}.${conf}
            eval_opt=""
            [[ -f "$EVALS_FILE" ]] && eval_opt="--evals $EVALS_FILE"

            ### Weight tracks: controlled by WEIGHTS env var (parsed above into _weight_labels).
            ### Pass the full token list to --weights so FieldDensityEigen knows which
            ### tracks to activate (sign, mgap, and/or any literal-value tracks).
            weights_opt="--weights $WEIGHTS"

            ### qlat reference SCIDAC files (optional; requires comp_file data + ordered evecs)
            comp_opt=""
            comp_files=""
            for idx in 0 1; do
                f=${COMP_DIR}/topo_field_${idx}.scidac
                [[ -f $f ]] && comp_files+=$f,
            done
            comp_files=${comp_files%?}   # strip trailing comma
            [[ -n "$comp_files" ]] && comp_opt="--comp_file $comp_files"

            ### Steps 1+2: compute topo fields + IP/corr vs gluonic TCD in one call
            ### --topo_out template: C++ substitutes {def} with A_mgap, B_mgap, C_mgap,
            ###                      A_sign, B_sign, C_sign  (6 output files per conf/tau)
            ### --comp_file (if present): compares each qlat field vs all 6 defs
            ### Output lines starting with "Topo PCF" and "CompRef" parsed below
            scratch=${HMC_DIR}/tmp_topo_pcf_${conf}_${tau}
            # MPI ranks write directly to fd 1, bypassing any shell pipe, so
            # piping through tr is ineffective.  Instead: redirect stdout to a
            # raw scratch file, then strip null bytes in a second pass.
            # Null bytes appear mid-line when Grid's LIME reader flushes a binary
            # read buffer to the same fd as GridLogMessage on Frontier.
            scratch_raw=${scratch}.raw
            ${CDIR}/FieldDensityEigen \
                --grid $vol \
                --files2 $F2s --Ls 48 $eval_opt $weights_opt \
                --topo_out ${DATA_DIR_topo}/Top_dnsty_q_{def}_${tau}_smr.${conf} \
                --files1 $F1s --topo_compare --conf_id $conf \
                $comp_opt \
                > $scratch_raw 2>&1
            tr -d '\000' < $scratch_raw > $scratch
            cat $scratch          # echo cleaned output to terminal/batch log
            rm -f $scratch_raw

            # Split PCF output into per-def files; active set matches WEIGHTS token list
            # Each block is preceded by "# q_X_wmode"; lines are "Topo PCF Corr/IP: TD_i conf val"
            # Reformat to: TD_tau  tau  conf  value  (matching notebook MultiIndex schema)
            for q_def in $_q_defs_all; do
                awk -v q=$q_def -v tau=$tau -v tds="0 4 16" '
                    /^# / { active=($2==q) }
                    active && /^Topo PCF Corr:/ { split(tds,td," "); print td[$3+1], tau, $4, $5 }
                    active && /^Topo PCF IP:/   { split(tds,td," "); print td[$3+1], tau, $4, $5 }
                ' $scratch >> ${HMC_DIR}/data/corr_ip_${q_def}.dat
            done

            # Parse CompRef lines (present only when --comp_file fired)
            # Line format: "CompRef: def comp_idx conf Q_evec Q_ref Corr IP rms_diff"
            # Output: comp_idx  tau  conf  Q_evec  Q_ref  Corr  IP  rms_diff
            for q_def in $_q_defs_all; do
                awk -v q=$q_def -v tau=$tau '
                    /^CompRef:/ && $2==q { print $3, tau, $4, $5, $6, $7, $8, $9 }
                ' $scratch >> ${HMC_DIR}/data/comp_ref_${q_def}.dat
            done

            # Parse TopoCompRef lines: stochastic qlat field vs gluonic TCD
            # Line format: "TopoCompRef: comp_idx TD_tau_idx conf Q_ref Q_gluon Corr IP"
            # Output: comp_idx  TD_tau  tau  conf  Q_ref  Q_gluon  Corr  IP
            awk -v tau=$tau -v tds="0 4 16" '
                /^TopoCompRef:/ { split(tds,td," "); print $2, td[$3+1], tau, $4, $5, $6, $7, $8 }
            ' $scratch >> ${HMC_DIR}/data/corr_ip_stoch.dat

            rm -f $scratch

            # Register topo density files for archiving
            pfx=${DATA_DIR_topo}/Top_dnsty_q
            for _lbl in "${_weight_labels[@]}"; do
                for _def in A B C; do
                    dfiles+=( ${HMC}/eigen/${conf}/Top_dnsty_q_${_def}_${_lbl}_${tau}_smr.${conf} )
                done
            done

            ### Step 3: T-animated movie — one column per weight track (A/B/C rows) + qlat reference
            Fs_all=""
            for _lbl in "${_weight_labels[@]}"; do
                for _def in A B C; do
                    Fs_all+="${pfx}_${_def}_${_lbl}_${tau}_smr.${conf},"
                done
            done
            for idx in 0 1; do
                f=${COMP_DIR}/topo_field_${idx}.scidac
                [[ -f $f ]] && Fs_all+="${f},"
            done
            Fs_all="${Fs_all%,}"   # strip trailing comma
            mpeg_all=${HMC_DIR}/Top_dnsty_all_defs_${conf}_tau${tau}.avi
            _first_lbl="${_weight_labels[0]}"
            if [[ -f ${pfx}_A_${_first_lbl}_${tau}_smr.${conf} ]]; then
                ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs_all --grid $vol --animate T \
                       --mpeg $mpeg_all --isosurface -0.01
            fi

        done   # tau
    fi     # REGEN
done       # i_conf


########################################################################################################################
####################  §3  Trajectory analysis — H_DWF evecs, config 702  #############################################
########################################################################################################################


###########   INPUT   ###########################
CONFS=(  702  719  71902 )
REGENS=(   0    1      1 )   # set 0 to skip a config
NCUT=4
#################################################

for i_conf in "${!CONFS[@]}"; do
    conf=${CONFS[$i_conf]}
    REGEN=${REGENS[$i_conf]}

DATA_DIR=${HMC_DIR}/eigen/${conf}


########################################################################################################################
### §3.1  Per-evec density movies (ZT-summed, ZT-updated, T-updated)
# Each mode n=0..NCUT is shown in four presentations:
#   (a) ZT-summed: 5D, sum Z+T → show (Ls,X,Y); tau_MD animated
#   (b) ZT-updated: 5D, cycle Z and T → show (Ls,X,Y) at each (Z,T); tau_MD animated
#   (c) Z-updated, T=23: 5D, cycle Z, fix T=23 → show (Ls,X,Y); tau_MD animated
#   (d) T-updated: 5D, sum Ls, cycle T → show (X,Y,Z) at each T; tau_MD animated
########################################################################################################################

for dof in smr lat; do
    for tau in 0 4; do

	#########   Separately for each eigenvector      #############

	for n in `seq 0 1 $NCUT`; do

	    fname=evec_density_sorted_${n}_tau_${tau}_${dof}

	    ### (a) (Ls, X, Y): Z,T summed

	    ext=${conf}_ZT_summed

	    mpeg=${HMC_DIR}/${fname}_${ext}.avi
            dpath=${HMC_DIR}/${fname}_${ext}.dat
            dfile=${HMC}/${fname}_${ext}.dat
            mfile=${HMC}/${fname}_${ext}.avi
            dfiles+=( $dfile ) #$mfile )
            echo $tau $dof $n

	    if [[ $REGEN == 0 ]] ; then
		F=""
		for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
		Fs=${F%?}
		#eigen/702/evec_density_sorted_0_tau_0_smr_702_3.345833
		iso="-0.01" #`awk -v t=$tau -v n=$i_evec -v N=$nconv 'BEGIN{print -1*( 0.3 + 0.4*t/4 + 0.9*n/N)}'`
		sep=${conf}_
		${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --animate configs --sum Z --sum T \
		       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
            fi

	    ### (b) (Ls, X, Y): Z,T looped (takes too much time — disabled)

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
                ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --animate configs --cycle Z=0 --cycle T=0 \
                       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
            fi

	    ### (c) (Ls, X, Y): Z looped, T=23 fixed

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
                ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --animate configs --cycle Z=0 --fix T=23 \
                       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
            fi

	    ### (d) (X, Y, Z): T looped; for comparison with H_W modes

	    ext=T_update

	    mpeg=${HMC_DIR}/${fname}_${ext}.avi
            dpath=${HMC_DIR}/${fname}_${ext}.dat
            dfile=${HMC}/${fname}_${ext}.dat
            mfile=${HMC}/${fname}_${ext}.avi
            #dfiles+=( $dfile ) #$mfile )

	    if [[ $REGEN == 0 ]] ; then
		F=""
                for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
                Fs=${F%?}

		iso="-0.05" #`awk -v t=$tau -v n=$i_evec -v N=$nconv 'BEGIN{print -1*( 0.3 + 0.4*t/4 + 0.9*n/N)}'`
                sep=${conf}_
                ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --Ls 48 --animate configs --sum Ls --cycle T=0 \
		       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
	    fi

	done


	########################################################################################################################
	### §3.2  Sum over all modes at each tau_MD snapshot (FieldDensityEigen --sum_all_files)
	# Sums evec density across all converged modes at a given tau_MD and dof, writing one file per snapshot.
	########################################################################################################################

	if [[ $REGEN == 0 ]] ; then
	    for t in `ls $DATA_DIR/evec_density_sorted_0_tau_${tau}_${dof}_${conf}_*| awk -F _ '{print $NF}' | sort -n`; do
		F=""
		for f in `ls $DATA_DIR/evec_density_sorted_*_tau_${tau}_${dof}_${conf}_${t} | sort -n` ; do F+=$f,; done
		Fs=${F%?}

		save_fname=$DATA_DIR/summed_evec_density_sorted_tau_${tau}_${dof}_${conf}_${t}
		${CDIR}/FieldDensityEigen --grid $vol --files2 $Fs --Ls 48 --sum_all_files $save_fname
	    done
	fi


	########################################################################################################################
	### §3.3  Summed-evec density movie (T-updated, 4D)
	# Animates summed evec density over trajectory; T cycles; for comparison with fermion force.
	########################################################################################################################

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
	    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --cycle T=0 \
		   --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
	fi
    done
done

done  # end CONFS loop


########################################################################################################################
####################  §4  Trajectory analysis — H_W (Wilson) evecs via specflow  #####################################
########################################################################################################################

############  INPUT ##################
CONF=702
tau=0
######################################

DATA_DIR=${HMC_DIR}/eigen_Wilson/${CONF}/${CONF}


########################################################################################################################
### §4.1  Sum evec density at each M_5 (initial config, all t_MD directories)
# For each t_MD step and each M_5 mass, sums all eigenvector density files via FieldDensityEigen.
########################################################################################################################

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


########################################################################################################################
### §4.2  Summed evec density movies (T-updated, initial t_MD and fixed m=-1.8)
# (a) Animates over M_5 at fixed t_MD=0; T cycles.
# (b) Animates over t_MD steps at fixed m=-1.8; T cycles.
########################################################################################################################

### (a) T-updated: animate over M_5 values at t_MD=0

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
    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --cycle T=0 \
	   --mpeg $mpeg --isosurface $iso --index_file foo_ind
fi

### (b) T-updated: animate over t_MD at fixed m=-1.8

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
    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --cycle T=0 \
	   --mpeg $mpeg --isosurface $iso --index_file foo_ind
fi


########################################################################################################################
### §4.3  0th evec density over trajectory at m=-1.8 (multi-config loop)
# For each config in CONFS: animate over M_5 (T-updated) and over t_MD (T-summed and T-updated) at m=-1.8.
########################################################################################################################

#########   INPUT   #################
CONFS=( 702 7026 70201 70202 70203 70204 70205 703 70301 718 719 )
regens=( 0   0     0     0     0     0     0    0    0    0   0  )
######################################

for((i_conf=0; i_conf<${#CONFS[@]}; i_conf++)); do
    CONF=${CONFS[i_conf]}
    DATA_DIR=${HMC_DIR}/eigen_Wilson/${CONF}/

    ### T-updated: animate over M_5 at t_MD=0

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
	${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --cycle T=0 \
	       --mpeg $mpeg --isosurface $iso --index_file foo_ind
    fi


    ### T-summed: animate over t_MD at m=-1.8

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
    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --sum T \
	   --mpeg $mpeg --isosurface $iso --save_data_to $dpath --index_file foo_ind
    fi

    ### T-updated: animate over t_MD at m=-1.8

    tau=0   # only tau=0 data available
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
	${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --cycle T=0 \
	       --mpeg $mpeg --isosurface $iso --save_data_to $dpath --index_file foo_ind
    fi

    # (X, Y, Z) + T looped: track specflow poles, visualize corresp. evecs — SKIP FOR NOW
done


########################################################################################################################
####################  §5  Gluonic E & TC density along trajectory (multi-config loop)  ###############################
########################################################################################################################

#########   INPUT   #################
CONFS=( 702 7026 70201 70202 70203 70204 70205 703 70301 718 719 )
regens=( 0   0     0     0     0     0     0    0    0    0   0  )
######################################

# Gluonic topological charge density files: ${HMC_DIR}/dnsty/${conf}/${fname}.${tau_MD}
# tau_MD is the molecular-dynamics time used as the animation frame index.

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
		    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --sum T --use_fname_as_frame_counter $sep \
			   --isosurface $iso --mpeg $mpeg --save_data_to $dpath #$(( (tau+1)*5 ))
		else
		    sep="${dof}."
		    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --cycle T=0 --isosurface $iso \
			   --mpeg $mpeg --use_fname_as_frame_counter $sep
		fi
	    fi
	done
    done
done


########################################################################################################################
####################  §6  Fermion force density along trajectory (multi-config loop)  ################################
########################################################################################################################

#########   INPUT   #################
CONFS=( 702 7026 70201 70202 70203 70204 70205 703 70301 718 719 )
regens=( 0   0     0     0     0     0     0    0    0    0   0  )
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

		${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --cycle T=0 --isosurface $iso \
		       --mpeg $mpeg --save_data_to $dpath --index_file foo_ind
	    fi
	done
    done
done


########################################################################################################################
####################  §7  Gauge action force density (Iwasaki, Jacobian)  ############################################
########################################################################################################################

for force in IwasakiGaugeAction JacobianAction; do

    if [[ "$force" == "IwasakiGaugeAction" ]] ; then iters="smr lat"; else iters=lat; fi
    for dof in $iters; do

	fname=F_${force}_${dof}
	ext=T_update #summed

	mpeg=${HMC_DIR}/${fname}_${ext}.avi
	dpath=${HMC_DIR}/${fname}_${ext}.dat
	dfile=${HMC}/${fname}_${ext}.dat
	mfile=${HMC}/${fname}_${ext}.avi
	dfiles+=( $dfile )   # archive .dat only; .avi files are too large

	if [[ "$dof" == "smr" ]] ; then iso=-0.45; else iso=-0.58; fi
	if [[ "$force" == "JacobianAction" ]] ; then iso=-0.61; fi

	if [[ $REGEN == 0 ]] ; then
            F=""
            for f in `ls $DATA_DIR/${conf}/${fname}.*|awk -F . '{print $NF, $0}' | sort  -nk1| cut -f2- -d' ' | tail -375 |head -100`; do F+=$f,;done
            Fs=${F%?}
	    tail -376 traj_times|head -100 > foo_ind #traj_times includes 4.00, at which force is not computed
	    ${CDIR}/FieldDensityAnimateMultiFiles --files $Fs --grid $vol --animate configs --cycle T=0 --isosurface $iso \
		   --mpeg $mpeg --save_data_to $dpath --index_file foo_ind
	fi
    done
done


########################################################################################################################
####################  §8  Gluonic TCD filtering via fermion zero modes  ###############################################
########################################################################################################################
# Compares gluonic topo charge density with the lowest H_DWF eigenvector densities
# using FieldDensityEigen --compareTCD_defs.  Output: filter_TCD.dat.
########################################################################################################################

DATA_DIR_dnsty=${HMC_DIR}/dnsty
DATA_DIR_eigen=${HMC_DIR}/eigen

dfile=${HMC}/data/filter_TCD.dat
dpath=${HMC_DIR}/data/filter_TCD.dat
dfiles+=( $dfile )

>$dpath
if [[ $REGEN == 0 ]] ; then
    for TD_tau in 0 4; do
	fname2=evec_density_0_tau_${TD_tau}
	F1s=""
	for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do  F1s+=$DATA_DIR_dnsty/Top_dnsty_${TD_tau}_ckpoint_EODWF_lat_smr.${conf},;done
	F1s=${F1s%?}
	F2s=""
	for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do  F2s+=$DATA_DIR_eigen/$conf/$fname2.${conf},;done
	F2s=${F2s%?}

	${CDIR}/FieldDensityEigen --files1 $F1s --files2 $F2s --grid $vol --Ls 48 --compareTCD_defs --cut 0.00001 | tee -a foo2 | cat
	grep "Filtered Sum" foo2 | awk -v tau=$TD_tau '{print tau, $3, $4, $8}' >> $dpath
    done
fi
wait


########################################################################################################################
####################  §9  Archive to Lustre  ##########################################################################
########################################################################################################################

tar -cvf ${LDIR}/eigen_data.tar -C ${PDIR} ${dfiles[@]}


### Refs
# https://zenn.dev/shuh/articles/tar-command-use
# https://stackoverflow.com/questions/50338201/how-to-compress-and-tar-a-folder-in-linux
# https://askubuntu.com/questions/392885/how-can-i-view-the-contents-of-tar-gz-file-without-extracting-from-the-command-l
