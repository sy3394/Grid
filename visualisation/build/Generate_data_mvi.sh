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

# Convenience wrappers — strip stdout null bytes only, pass everything else
# through.  c-lime / Grid SCIDAC I/O emits NUL bytes via direct write()
# syscalls that bypass C-level redirects (a known issue on Cray/HPE — see
# git history of "Fix binary garbage" / "Eliminate binary garbage" commits).
# `tr -d '\0'` drops the NULs and keeps all printable text intact, so log_G
# remains greppable while Topo PCF / FermFerm / AlphaSweep diagnostics still
# reach the operator.
# §2.5 Step 2 (SCIDAC-write call) uses ${CDIR}/FieldDensityEigen directly
# with its own > $topo_write_log redirect, so it bypasses this filter; its
# binary lands in topo_write_log, never in log_G.
FDE()  { ${CDIR}/FieldDensityEigen            "$@" 2> >(tr -d '\0' >&2) | tr -d '\0'; }
FDAM() { ${CDIR}/FieldDensityAnimateMultiFiles "$@" 2> >(tr -d '\0' >&2) | tr -d '\0'; }
PDIR=/ccs/home/syamamoto/tmp/src/Grid_cleanedup_for_pullrequest/systems/Frontier/HMC
LDIR=/lustre/orion/phy157/proj-shared/phy157_dwf/syamamoto
HMC=32cube-rho0.124-tau4
HMC_DIR=$PDIR/$HMC
vol=32.32.32.32
REGEN=1
REGEN_MOVIE=0   # set to 1 to also generate the T-animated multi-panel movie in §2.5 Step 3
PANEL_PX=1536   # per-panel pixel size for FDAM movies (total = PANEL_PX * ncols × PANEL_PX * nrows)


########################################################################################################################
####################  §2  Ensemble averages (configs CONF_S .. CONF_F)  ##############################################
########################################################################################################################

########    INPUT: DDIR dependent   #########
CONF_S=700
CONF_F=709
nconv=15  # number of converged H_DWF eigenvectors

# ── Two distinct fermion masses are needed by FieldDensityEigen ────────────────
# MASS_EVEC: the m_f used by Compute_DWF_G5R5.cc when generating the eigenvectors.
#   Used to recover the kinetic eigenvalue lambda^(0) = sqrt((lambda^H)^2 - m_f^2).
#   Must MATCH the value used to generate the eigenvectors, else q_naive's second
#   term has a wrong kinetic eigenvalue.  These eigenvectors were generated at m=0.
MASS_EVEC=0
# BC_MASS: the m_f used as Banks–Casher Lorentzian regulator inside Sigma_low(x):
#     Sigma_low(x) = sum_n [ m_f / ((lambda^(0)_n)^2 + m_f^2) ] rho_n(x)
#   This is a *probe* mass for the §4.1 scalar-foil diagnostic, independent of
#   the evec generation mass.  Setting BC_MASS=0 collapses Sigma_low to zero
#   (the Banks–Casher chiral limit Sigma=pi*rho(0) is delta-function-like and
#   not pointwise-defined).  Pick a small physical-ish value, e.g. 0.01.
BC_MASS=0.01
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
	#dfiles+=( $dfile )   # too large — skip tarball

	if [[ $REGEN == 0 ]] ; then
	    F=""
	    for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do F+=$DATA_DIR/${fname}.${conf},;done
	    Fs=${F%?}
	    seq -f "%03g" $((CONF_S)) 1 $CONF_F > foo_ind

	    FDAM --files $Fs --grid $vol --animate configs --sum T \
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
	#dfiles+=( $dfile )   # too large — skip tarball

	if [[ $REGEN == 0 ]] ; then

	    F=""
	    for conf in `seq -f "%03g" $((CONF_S)) 1 $CONF_F`; do F+=$DATA_DIR/${conf}/${fname}.${conf},;done
	    Fs=${F%?}
	    seq -f "%03g" $((CONF_S)) 1 $CONF_F > foo_ind

	    sf=`awk -v t=$tau -v n=$i_evec -v N=$nconv 'BEGIN{print -1*( 0.3 + 0.4*t/4 + 0.9*n/N)}'`
	    FDAM --files $Fs --grid $vol --Ls 48 --animate configs --sum Ls --sum T \
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
	#dfiles+=( $dfile )   # too large — skip tarball

	if [[ $REGEN == 0 ]] ; then
	    F=""
	    for f in `ls $DATA_DIR/${conf}/${fname}_* | grep tau_${tau} | sort -n -t _ -k3`; do F+=$f,; done
	    Fs=${F%?}
	    seq -f "%03g" $((CONF_S)) 1 $CONF_F > foo_ind

	    FDAM --files $Fs --grid $vol --Ls 48 --animate configs --sum Z --sum T \
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
CONFS=(    700  701  702  703  704  705  706  707  708  709  795 )
REGENS_25=(  0    0    1    0    0    0    0    0    0    0    1 )  # 1=run, 0=skip
# NTOPO_25[i] = number of leading near-zero modes to exclude from the direct
# bulk remainder B_bulk = sum_{n in bulk} chi_n^B (passed as --n_topo to
# FieldDensityEigen).  Set = |Q| of each config (or 2|Q| if the solver lists
# each zero mode as a degenerate +-pair); VERIFY against the per-mode
# "TopoContrib evec= c mu_n=" printout that the first NTOPO modes are the
# |mu_n|~m_gap ones.  0 => B_bulk is the full unweighted chiral sum.
NTOPO_25=(   1    1    1    1    1    1    1    1    1    1    1 )  # |Q| per conf
# 702: Q=-1; has comp_file data (topo_field_*.scidac) and correctly ordered evec files.
# 795: Q=+1; the sign-flip falsification config (sec:autocorr / sec:sp_sum_bonus
#      predict q_B^mgap uniformly POSITIVE and PCF(q_B^mgap,Sigma_low) ~ +0.63).
#      Evecs generated with HMC/Compute_DWF_G5R5.cc; needs its own gluonic
#      Top_dnsty_*.795 (Wilson flow) and, optionally, its own Luchang reference
#      topo_field_*_795.scidac for the comp comparison (see conf-aware lookup below).
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

# All append-mode .dat files written by FieldDensityEigen are truncated here at
# script start (>file).  This way reruns do not accumulate duplicate rows; the
# user does not have to manually `rm` anything between runs.
#
# Truncate-and-register pattern: ">file" empties the file (creating if absent),
# then dfiles+=(file) registers it for the final tarball.

# corr_ip_q_*.dat — per-(weight,definition) PCF/IP vs gluonic TCD at each TD_tau.
# Format: comp_idx  conf  Corr  (one section per TD_tau, separated by labels)
# B_bulk and q_B_unit are the direct bulk remainder and the unit-weight bulk
# estimator; FieldDensityEigen writes their corr_ip/comp_ref like the rest, so
# truncate and register them here too (else they accumulate and miss the tarball).
for q_def in $_q_defs_all q_naive Sigma_low B_bulk q_B_unit; do
    >${HMC_DIR}/data/corr_ip_${q_def}.dat
    dfiles+=( ${HMC}/data/corr_ip_${q_def}.dat )
done

# comp_ref_*.dat — per-(weight,definition) comparison vs qlat reference field(s).
# Format: comp_idx  tau_wf  conf  Q_evec  Q_ref  Corr  IP  rms_diff
for q_def in $_q_defs_all q_naive Sigma_low B_bulk q_B_unit; do
    >${HMC_DIR}/data/comp_ref_${q_def}.dat
    dfiles+=( ${HMC}/data/comp_ref_${q_def}.dat )
done

# Band-pass (--band_pass) corr_ip/comp_ref files have dynamic names
# (q_{A,B,C}_bp<lc>); clear any stale ones so the per-run append starts fresh.
# They are globbed back into the tarball just before tar (see end of script).
if [[ -n "$BAND_PASS" ]]; then
    rm -f ${HMC_DIR}/data/corr_ip_q_*_bp*.dat ${HMC_DIR}/data/comp_ref_q_*_bp*.dat
fi

# corr_ip_stoch.dat — qlat stochastic field vs gluonic TCD at each TD_tau.
# Format: comp_idx  TD_tau  tau_wf  conf  Q_ref  Q_gluon  Corr  IP
>${HMC_DIR}/data/corr_ip_stoch.dat
dfiles+=( ${HMC}/data/corr_ip_stoch.dat )

# alpha_sweep.dat — PCF/IP of q_B^mgap + alpha*B vs gluonic q at each TD_tau.
# Format: alpha  tau_wf  TD_tau  conf  Corr  IP
>${HMC_DIR}/data/alpha_sweep.dat
dfiles+=( ${HMC}/data/alpha_sweep.dat )

# fermferm.dat — pairwise PCF/IP between fermion-q definitions.
# Format: def1  def2  tau_wf  conf  Corr  IP
>${HMC_DIR}/data/fermferm.dat
dfiles+=( ${HMC}/data/fermferm.dat )

# bk_sweep.dat — cumulative-spectral-sum PCF (sec:bk_sweep): how concentrated
# in the lowest modes is the chirality structure of B = q_naive - q_B^mgap?
# Emitted only when --bk_sweep is passed to FieldDensityEigen.
# Format: k  tau_wf  TD_tau  conf  Corr  IP
>${HMC_DIR}/data/bk_sweep.dat
dfiles+=( ${HMC}/data/bk_sweep.dat )

# smear_sweep.dat — Gaussian density-smearing sweep (sec:smear_sweep,
# Luchang's diagnostic): apply 4D Gaussian smear of width sigma to BOTH
# fermionic q and gluonic q (and q_L if --comp_file present), measure PCF
# vs sigma.  Emitted only when --smear_sweep is passed.
# Format: kind  name  sigma  tau_wf  TD_tau  conf  Corr  IP
#   kind in {fermion, comp};  name is q_B_mgap, q_naive, comp_0, comp_1, ...
>${HMC_DIR}/data/smear_sweep.dat
dfiles+=( ${HMC}/data/smear_sweep.dat )

# Optional flags for the new diagnostics — controlled by env vars.
#   BK_SWEEP (default 1)              -> pass --bk_sweep to FieldDensityEigen;
#                                        set BK_SWEEP=0 to disable.
#   SMEAR_SWEEP (default "0.5,1,1.5,2,3,5")
#                                     -> pass --smear_sweep "<value>" (comma-sep);
#                                        set SMEAR_SWEEP="" to disable.
#   SIGNED_MGAP=1         -> pass --signed_mgap (LEGACY signed m_gap/mu_n weight,
#                           sign-blind to Q; default is the |mu_n| weight). Use
#                           only for cross-checks.
#   BAND_PASS="0.05,0.1,0.2" -> pass --band_pass "<lc,...>": add Gaussian even-window
#                           weight tracks w(mu)=exp(-mu^2/2 lc^2), one per lambda_c.
#                           Integer-exact members of the bulk-estimator family; add
#                           q_{A,B,C}_bp<lc> rows to the corr_ip PCF output (Step 1)
#                           and SCIDAC density fields (Step 2).  Empty = off.
#   PER_MODE_OUT=1        -> pass --per_mode_out: write per-mode mu_n, int chi_n^B,
#                           int rho_n, int rho_n^2, IPR_n to the Step-2 topo_out
#                           prefix "...permode.dat" (mobility-edge / R5 tests).
BK_SWEEP="${BK_SWEEP:-1}"
SMEAR_SWEEP="${SMEAR_SWEEP:-0.5,1,1.5,2,3,5}"
SIGNED_MGAP="${SIGNED_MGAP:-0}"
BAND_PASS="${BAND_PASS:-0.005,0.01,0.02,0.03,0.05,0.07,0.1,0.15,0.2,0.3}"  # 10-pt scan into+above the localized band; set "" to disable
PER_MODE_OUT="${PER_MODE_OUT:-1}"        # default ON; set PER_MODE_OUT=0 to disable
bk_opt=""
smr_opt=""
smgap_opt=""
bp_opt=""
pm_opt=""
[[ "$BK_SWEEP" == "1" ]] && bk_opt="--bk_sweep"
[[ -n "$SMEAR_SWEEP" ]] && smr_opt="--smear_sweep $SMEAR_SWEEP"
[[ "$SIGNED_MGAP" == "1" ]] && smgap_opt="--signed_mgap"
[[ -n "$BAND_PASS" ]] && bp_opt="--band_pass $BAND_PASS"
[[ "$PER_MODE_OUT" == "1" ]] && pm_opt="--per_mode_out"

for i_conf in "${!CONFS[@]}"; do
    conf=${CONFS[$i_conf]}
    regen=${REGENS_25[$i_conf]}
    n_topo=${NTOPO_25[$i_conf]:-1}      # |Q| for this conf (B_bulk zero-mode exclusion)
    if [[ $regen == 1 ]] ; then

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

            ### Gluonic TCD files at TD_tau = 0, 4, 16  (smr hardcoded: actual Frontier naming)
            ### Build the TD_tau label list in LOCKSTEP with the files that actually
            ### exist and pass it as --td_taus below.  --td_taus is POSITIONAL, so a
            ### missing file (e.g. no Top_dnsty_0 for conf 795) would otherwise
            ### mislabel the survivors (Top_dnsty_4 tagged "0", Top_dnsty_16 tagged
            ### "4").  Keeping F1 and TDL in step guarantees every PCF row carries
            ### its true tau_WG.
            F1=""; TDL=""
            for TD_tau in 0 4 16; do
                f=${DATA_DIR_dnsty}/Top_dnsty_${TD_tau}_ckpoint_EODWF_lat_smr.${conf}
                if [[ -f $f ]]; then F1+=$f,; TDL+=$TD_tau,;
                else echo "Warning: gluonic TCD not found: $f"; fi
            done
            F1s=${F1%?}
            TDLs=${TDL%?}

            ### qlat reference SCIDAC files (optional; requires comp_file data + ordered evecs)
            ### Conf-aware lookup: prefer topo_field_<idx>_<conf>.scidac so each
            ### config compares against ITS OWN Luchang reference.  The legacy
            ### unsuffixed topo_field_<idx>.scidac is conf 702's reference, so it
            ### is used only as a fallback for conf 702 — never for other confs
            ### (otherwise e.g. 795 would be compared against 702's field, which
            ### is physically meaningless).
            comp_opt=""
            comp_files=""
            for idx in 0 1; do
                f=${COMP_DIR}/topo_field_${idx}_${conf}.scidac
                [[ ! -f $f && "$conf" == "702" ]] && f=${COMP_DIR}/topo_field_${idx}.scidac
                [[ -f $f ]] && comp_files+=$f,
            done
            comp_files=${comp_files%?}   # strip trailing comma
            [[ -n "$comp_files" ]] && comp_opt="--comp_file $comp_files"

            ### Step 1b (tau=0 only, no evec files): pure qlat-vs-gluonic comparison.
            ### When no evec density files exist for tau=0 (F2s empty), provide the
            ### τ_FW=0 baseline: Luchang reference vs gluonic TCD without any eigenvector
            ### smearing.  This is the diagnostic from the analysis document — comparing
            ### the qlat stochastic field against Wilson-flowed gluonic TCD at tau_WF=0
            ### to separate resolution mismatch from operator-definition disagreement.
            ### Writes to corr_ip_stoch.dat only (no fermion-q fields available).
            ### (When F2s is non-empty, Step 1 handles this comparison via do_gluon_comp.)
            if [[ "$tau" == "0" && -z "$F2s" && -n "$comp_opt" && -n "$F1s" ]]; then
                echo "Step 1b: tau=0 qlat-vs-gluonic (no evec files): conf=$conf"
                FDE \
                    --grid $vol \
                    --files1 $F1s --topo_compare --conf_id $conf \
                    --tau_wf $tau --td_taus ${TDLs:-0,4,16} \
                    --data_dir ${HMC_DIR}/data \
                    $comp_opt
            fi

            ### Skip eigenvector-dependent steps if no evec files found
            [[ -z "$F2s" ]] && { echo "No evec files for conf=$conf tau=$tau, skipping eigenvec steps"; continue; }

            ### Eigenvalue file (enables m_gap/mu_n weighting; written by Compute_DWF_G5R5)
            ### For tau=0: if eigenvalues_tau_0 is absent but eigenvalues_tau_4 exists,
            ### fall back to tau=4 evals as a proxy (eigenvalues shift <1% with tau_WF).
            EVALS_FILE=${DATA_DIR_eigen}/eigenvalues_tau_${tau}.${conf}
            eval_opt=""
            if [[ -f "$EVALS_FILE" ]]; then
                eval_opt="--evals $EVALS_FILE"
            elif [[ "$tau" == "0" && -f "${DATA_DIR_eigen}/eigenvalues_tau_4.${conf}" ]]; then
                echo "Note: conf=$conf eigenvalues_tau_0 absent — using tau=4 evals as proxy"
                eval_opt="--evals ${DATA_DIR_eigen}/eigenvalues_tau_4.${conf}"
            else
                ### No per-config eigenvalue file: q_B^mgap/q_B^sign weights fall back
                ### to mu_n=0 and the SIGN-weighted estimators are MISCOMPUTED on
                ### topological configs.  Warn loudly so it is never silent.
                ### (NB: this is a DIFFERENT failure mode from the conf-795 sign-flip
                ###  anomaly, where the eval file is present but the mu~0 zero mode's
                ###  sign is convention-determined; see the predictions appendix.)
                echo "WARNING: conf=$conf tau=$tau — eigenvalue file '$EVALS_FILE' MISSING."
                echo "         sgn(mu_n)/m_gap weights unavailable: q_B^{mgap,sign} will be"
                echo "         miscomputed (sign-blind).  Generate eigenvalues_tau_${tau}.${conf}"
                echo "         (Compute_DWF_G5R5) before trusting the Gamma_5-weighted fields."
            fi

            ### Weight tracks: controlled by WEIGHTS env var (parsed above into _weight_labels).
            ### Pass the full token list to --weights so FieldDensityEigen knows which
            ### tracks to activate (sign, mgap, and/or any literal-value tracks).
            weights_opt="--weights $WEIGHTS"

            ### Step 1: IP/corr stats — text only, no SCIDAC writes, no LIME binary.
            ### Writes corr_ip_q_*.dat, comp_ref_q_*.dat, corr_ip_stoch.dat directly.
            ### When eval_opt is empty (old evec files without embedded eigenvalues and
            ### no tau=4 fallback), the fermion-q comparisons are skipped with a warning,
            ### but the qlat-vs-gluonic comparison (TopoCompRef) still runs via do_gluon_comp.
            FDE \
                --grid $vol \
                --files2 $F2s --Ls 48 $eval_opt $weights_opt \
                --mass $MASS_EVEC --bc_mass $BC_MASS --n_topo $n_topo $smgap_opt \
                --files1 $F1s --topo_compare --conf_id $conf \
                --tau_wf $tau --td_taus ${TDLs:-0,4,16} \
                --data_dir ${HMC_DIR}/data \
                $comp_opt $bk_opt $smr_opt $bp_opt

            ### Step 2: write topo density SCIDAC fields — binary LIME output.
            ### Redirected to a dedicated per-conf/tau log to keep log_G clean.
            topo_write_log=${HMC_DIR}/tmp_topo_write_${conf}_${tau}.log
            ${CDIR}/FieldDensityEigen \
                --grid $vol \
                --files2 $F2s --Ls 48 $eval_opt $weights_opt \
                --mass $MASS_EVEC --bc_mass $BC_MASS --n_topo $n_topo $smgap_opt \
                --topo_out ${DATA_DIR_topo}/Top_dnsty_q_{def}_${tau}_smr.${conf} \
                $bp_opt $pm_opt \
                > $topo_write_log 2>&1

            # NOTE: the topo-density SCIDAC fields (q_{A,B,C}_{sign,mgap},
            # q_naive, q_Sigma, q_Bbulk) are LARGE binaries (~16 MB each) and are
            # transferred by direct rsync of eigen/<conf>/, NOT bundled into the
            # .dat tarball.  They are therefore deliberately NOT added to dfiles
            # (which is the tarball list, see below).  Cf. the skipped E_dnsty/
            # Top_dnsty/evec_density entries above.
            pfx=${DATA_DIR_topo}/Top_dnsty_q

            ### Step 3: T-animated movie — one column per weight track (A/B/C rows) + qlat reference
            Fs_all=""
            for _lbl in "${_weight_labels[@]}"; do
                for _def in A B C; do
                    Fs_all+="${pfx}_${_def}_${_lbl}_${tau}_smr.${conf},"
                done
            done
            # q_naive (legacy sp_sum, sign-correct), Sigma_low (Banks-Casher
            # locality), and Bbulk (direct bulk remainder): single fields, no
            # weight-track loop. Added once per tau.
            for _name in naive Sigma Bbulk; do
                f="${pfx}_${_name}_${tau}_smr.${conf}"
                [[ -f $f ]] && Fs_all+="${f},"
            done
            for idx in 0 1; do
                f=${COMP_DIR}/topo_field_${idx}.scidac
                [[ -f $f ]] && Fs_all+="${f},"
            done
            Fs_all="${Fs_all%,}"   # strip trailing comma
            mpeg_all=${HMC_DIR}/Top_dnsty_all_defs_${conf}_tau${tau}.avi
            _first_lbl="${_weight_labels[0]}"
            if [[ $REGEN_MOVIE == 1 && -f ${pfx}_A_${_first_lbl}_${tau}_smr.${conf} ]]; then
                FDAM --files $Fs_all --grid $vol --animate T \
                       --mpeg $mpeg_all --isosurface -0.01 \
                       --panel_size ${PANEL_PX}
            fi

        done   # tau
    fi     # regen
done       # i_conf


########################################################################################################################
####################  §3  Trajectory analysis — H_DWF evecs, config 702  #############################################
########################################################################################################################


###########   INPUT   ###########################
CONFS=(  702  719  71902 )
regens=(   0    0      0 )   # set 1 to run a config
NCUT=4
#################################################

for i_conf in "${!CONFS[@]}"; do
    conf=${CONFS[$i_conf]}
    regen=${regens[$i_conf]}

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
            #dfiles+=( $dfile )   # too large / not always generated — skip tarball
            echo $tau $dof $n

	    if [[ $regen == 1 ]] ; then
		F=""
		for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
		Fs=${F%?}
		#eigen/702/evec_density_sorted_0_tau_0_smr_702_3.345833
		iso="-0.01" #`awk -v t=$tau -v n=$i_evec -v N=$nconv 'BEGIN{print -1*( 0.3 + 0.4*t/4 + 0.9*n/N)}'`
		sep=${conf}_
		FDAM --files $Fs --grid $vol --Ls 48 --animate configs --sum Z --sum T \
		       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
            fi

	    ### (b) (Ls, X, Y): Z,T looped (takes too much time — disabled)

	    ext=${conf}_ZT_update

            mpeg=${HMC_DIR}/${fname}_${ext}.avi
            dpath=${HMC_DIR}/${fname}_${ext}.dat
            dfile=${HMC}/${fname}_${ext}.dat
            mfile=${HMC}/${fname}_${ext}.avi
            #dfiles+=( $dfile )

            if [[ $regen == 1 ]] ; then
                F=""
                for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
                Fs=${F%?}

                iso="-0.01"
                sep=${conf}_
                FDAM --files $Fs --grid $vol --Ls 48 --animate configs --cycle Z=0 --cycle T=0 \
                       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
            fi

	    ### (c) (Ls, X, Y): Z looped, T=23 fixed

	    ext=${conf}_Z_update_T23

            mpeg=${HMC_DIR}/${fname}_${ext}.avi
            dpath=${HMC_DIR}/${fname}_${ext}.dat
            dfile=${HMC}/${fname}_${ext}.dat
            mfile=${HMC}/${fname}_${ext}.avi
            #dfiles+=( $dfile )

            if [[ $regen == 1 ]] ; then
                F=""
                for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
                Fs=${F%?}

                iso="-0.01"
                sep=${conf}_
                FDAM --files $Fs --grid $vol --Ls 48 --animate configs --cycle Z=0 --fix T=23 \
                       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
            fi

	    ### (d) (X, Y, Z): T looped; for comparison with H_W modes

	    ext=T_update

	    mpeg=${HMC_DIR}/${fname}_${ext}.avi
            dpath=${HMC_DIR}/${fname}_${ext}.dat
            dfile=${HMC}/${fname}_${ext}.dat
            mfile=${HMC}/${fname}_${ext}.avi
            #dfiles+=( $dfile ) #$mfile )

	    if [[ $regen == 1 ]] ; then
		F=""
                for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
                Fs=${F%?}

		iso="-0.05" #`awk -v t=$tau -v n=$i_evec -v N=$nconv 'BEGIN{print -1*( 0.3 + 0.4*t/4 + 0.9*n/N)}'`
                sep=${conf}_
                FDAM --files $Fs --grid $vol --Ls 48 --animate configs --sum Ls --cycle T=0 \
		       --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
	    fi

	done


	########################################################################################################################
	### §3.2  Sum over all modes at each tau_MD snapshot (FieldDensityEigen --sum_all_files)
	# Sums evec density across all converged modes at a given tau_MD and dof, writing one file per snapshot.
	########################################################################################################################

	if [[ $regen == 1 ]] ; then
	    for t in `ls $DATA_DIR/evec_density_sorted_0_tau_${tau}_${dof}_${conf}_*| awk -F _ '{print $NF}' | sort -n`; do
		F=""
		for f in `ls $DATA_DIR/evec_density_sorted_*_tau_${tau}_${dof}_${conf}_${t} | sort -n` ; do F+=$f,; done
		Fs=${F%?}

		save_fname=$DATA_DIR/summed_evec_density_sorted_tau_${tau}_${dof}_${conf}_${t}
		FDE --grid $vol --files2 $Fs --Ls 48 --sum_all_files $save_fname
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

	if [[ $regen == 1 ]] ; then
	    F=""
            for f in `ls $DATA_DIR/${fname}_* | sort -n` ; do F+=$f,; done
            Fs=${F%?}

	    iso="-0.05"
	    sep=${conf}_
	    FDAM --files $Fs --grid $vol --animate configs --cycle T=0 \
		   --mpeg $mpeg --isosurface $iso --use_fname_as_frame_counter $sep
	fi
    done
done

done  # end CONFS loop §3


########################################################################################################################
####################  §4  Trajectory analysis — H_W (Wilson) evecs via specflow  #####################################
########################################################################################################################

############  INPUT ##################
CONF=702
tau=0
regen_4=0   # set 1 to run §4.1/§4.2 (independent of global REGEN)
######################################

DATA_DIR=${HMC_DIR}/eigen_Wilson/${CONF}/${CONF}


########################################################################################################################
### §4.1  Sum evec density at each M_5 (initial config, all t_MD directories)
# For each t_MD step and each M_5 mass, sums all eigenvector density files via FieldDensityEigen.
########################################################################################################################

if [[ $regen_4 == 1 ]] ; then
    for d in `ls -d $DATA_DIR/U_smr_*| sort -n`; do
	for m5 in `ls $d/evec_*_0| awk -F _ '{print $(NF-1)}'`; do
	    F=""
	    for f in `ls $d/evec_${m5}_* | sort -n` ; do F+=$f,; done
	    Fs=${F%?}

	    save_fname=$d/evec_sum_${m5}_${CONF}
	    FDE --files1 $Fs --grid $vol --sum_all_files $save_fname

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

if [[ $regen_4 == 1 ]] ; then

    F=""
    for f in `ls $DATA_DIR/U_smr_${t_MD}/evec_sum_* | sort -n` ; do F+=$f,; done
    Fs=${F%?}

    ls $DATA_DIR/U_smr_${t_MD}/evec_sum_* | awk -F _ '{print $NF}' | sort -n > foo_ind
    iso="-0.01"
    FDAM --files $Fs --grid $vol --animate configs --cycle T=0 \
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

if [[ $regen_4 == 1 ]] ; then
    F=""
    for f in `ls $DATA_DIR/U_smr_*/evec_sum_-1.800000 | sort -n` ; do F+=$f,; done
    Fs=${F%?}

    ls $DATA_DIR/U_smr_*/evec_sum_-1.800000 | awk -F _ '{print $(NF-2)}' | awk -F '/' '{print $1}' | sort -n > foo_ind
    iso="-0.01"
    FDAM --files $Fs --grid $vol --animate configs --cycle T=0 \
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

    if [[ ${regens[i_conf]} == 1 ]] ; then

	F=""
	for f in `ls $DATA_DIR/U_smr_${t_MD}/evec_*_0 | sort -n` ; do F+=$f,; done
	Fs=${F%?}

	ls $DATA_DIR/U_smr_${t_MD}/evec_*_0 | awk -F _ '{print $(NF-1)}' | sort -n > foo_ind
	iso="-0.01"
	FDAM --files $Fs --grid $vol --animate configs --cycle T=0 \
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
    if [[ ${regens[i_conf]} == 1 ]] ; then
    F=""
    for f in `ls $DATA_DIR/U_smr_*/evec_-1.800000_0 | sort -n` ; do F+=$f,; done
    Fs=${F%?}

    ls $DATA_DIR/U_smr_*/evec_-1.800000_0 | awk -F _ '{print $(NF-2)}' | awk -F '/' '{print $1}' | sort -n > foo_ind
    iso="-0.01"
    FDAM --files $Fs --grid $vol --animate configs --sum T \
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

    if [[ ${regens[i_conf]} == 1 ]] ; then
	F=""
	for f in `ls $DATA_DIR/U_smr_*/evec_-1.800000_0 | sort -n` ; do F+=$f,; done
	Fs=${F%?}

	ls $DATA_DIR/U_smr_*/evec_-1.800000_0 | awk -F _ '{print $(NF-2)}' | awk -F '/' '{print $1}' | sort -n > foo_ind
	iso="-0.01"
	FDAM --files $Fs --grid $vol --animate configs --cycle T=0 \
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

	    if [[ ${regens[i_conf]} == 1 ]] ; then
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
		    FDAM --files $Fs --grid $vol --animate configs --sum T --use_fname_as_frame_counter $sep \
			   --isosurface $iso --mpeg $mpeg --save_data_to $dpath #$(( (tau+1)*5 ))
		else
		    sep="${dof}."
		    FDAM --files $Fs --grid $vol --animate configs --cycle T=0 --isosurface $iso \
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
	    #dfiles+=( $dfile )   # too large / not always generated — skip tarball

	    if [[ ${regens[i_conf]} == 1 ]] ; then

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

		FDAM --files $Fs --grid $vol --animate configs --cycle T=0 --isosurface $iso \
		       --mpeg $mpeg --save_data_to $dpath --index_file foo_ind
	    fi
	done
    done
done


########################################################################################################################
####################  §7  Gauge action force density (Iwasaki, Jacobian)  ############################################
########################################################################################################################

#########   INPUT   #################
CONFS=( 702 7026 70201 70202 70203 70204 70205 703 70301 718 719 )
regens=( 0   0     0     0     0     0     0    0    0    0   0  )
######################################

DATA_DIR=${HMC_DIR}/snapshots

for((i_conf=0; i_conf<${#CONFS[@]}; i_conf++)); do
    conf=${CONFS[i_conf]}

    for force in IwasakiGaugeAction JacobianAction; do

        if [[ "$force" == "IwasakiGaugeAction" ]] ; then iters="smr lat"; else iters=lat; fi
        for dof in $iters; do

            fname=F_${force}_${dof}
            ext=T_update #summed

            mpeg=${HMC_DIR}/${fname}_${ext}.avi
            dpath=${HMC_DIR}/${fname}_${ext}.dat
            dfile=${HMC}/${fname}_${ext}.dat
            #dfiles+=( $dfile )   # too large / not always generated — skip tarball

            if [[ "$dof" == "smr" ]] ; then iso=-0.45; else iso=-0.58; fi
            if [[ "$force" == "JacobianAction" ]] ; then iso=-0.61; fi

            if [[ ${regens[i_conf]} == 1 ]] ; then
                F=""
                for f in `ls $DATA_DIR/${conf}/${fname}.*|awk -F . '{print $NF, $0}' | sort  -nk1| cut -f2- -d' ' | tail -375 |head -100`; do F+=$f,;done
                Fs=${F%?}
                tail -376 traj_times|head -100 > foo_ind #traj_times includes 4.00, at which force is not computed
                FDAM --files $Fs --grid $vol --animate configs --cycle T=0 --isosurface $iso \
                       --mpeg $mpeg --save_data_to $dpath --index_file foo_ind
            fi
        done
    done

done   # i_conf §7


########################################################################################################################
####################  §8  Gluonic TCD filtering via fermion zero modes  ###############################################
########################################################################################################################
# Compares gluonic topo charge density with the lowest H_DWF eigenvector densities
# using FieldDensityEigen --compareTCD_defs.  Output: filter_TCD.dat.
########################################################################################################################

#########   INPUT   #################
CONFS=( 702 7026 70201 70202 70203 70204 70205 703 70301 718 719 )
regens=( 0   0     0     0     0     0     0    0    0    0   0  )
######################################

DATA_DIR_dnsty=${HMC_DIR}/dnsty
DATA_DIR_eigen=${HMC_DIR}/eigen

dfile=${HMC}/data/filter_TCD.dat
dpath=${HMC_DIR}/data/filter_TCD.dat
dfiles+=( $dfile )

>$dpath

for((i_conf=0; i_conf<${#CONFS[@]}; i_conf++)); do
    conf=${CONFS[i_conf]}
    if [[ ${regens[i_conf]} == 1 ]] ; then
        for TD_tau in 0 4; do
            fname2=evec_density_0_tau_${TD_tau}
            f1=$DATA_DIR_dnsty/Top_dnsty_${TD_tau}_ckpoint_EODWF_lat_smr.${conf}
            f2=$DATA_DIR_eigen/${conf}/${fname2}.${conf}
            [[ ! -f $f1 ]] && { echo "Warning §8: gluonic TCD not found: $f1"; continue; }
            [[ ! -f $f2 ]] && { echo "Warning §8: evec file not found: $f2"; continue; }

            ${CDIR}/FieldDensityEigen --files1 $f1 --files2 $f2 --grid $vol --Ls 48 --compareTCD_defs --cut 0.00001 | tee -a foo2 | cat
            grep "Filtered Sum" foo2 | awk -v tau=$TD_tau -v c=$conf '{print tau, c, $3, $4, $8}' >> $dpath
        done
    fi
done   # i_conf §8

wait


########################################################################################################################
####################  §9  Archive to Lustre  ##########################################################################
########################################################################################################################

# Register opt-in diagnostic outputs with dynamic names into the tarball.
#   BAND_PASS: corr_ip_q_{A,B,C}_bp<lc>.dat / comp_ref_*  (PCF of the band-pass
#              family members vs gluonic/qlat reference).
#   PER_MODE_OUT: Top_dnsty_q_permode_<tau>_smr.<conf>.dat under eigen/<conf>/.
if [[ -n "$BAND_PASS" ]]; then
    for f in ${HMC_DIR}/data/corr_ip_q_*_bp*.dat ${HMC_DIR}/data/comp_ref_q_*_bp*.dat; do
        [[ -f $f ]] && dfiles+=( ${HMC}/data/$(basename "$f") )
    done
fi
if [[ "$PER_MODE_OUT" == "1" ]]; then
    for f in ${HMC_DIR}/eigen/*/Top_dnsty_q_permode_*.dat; do
        [[ -f $f ]] && dfiles+=( ${HMC}/eigen/$(basename "$(dirname "$f")")/$(basename "$f") )
    done
fi

tar -cvf ${LDIR}/eigen_data.tar -C ${PDIR} ${dfiles[@]}


### Refs
# https://zenn.dev/shuh/articles/tar-command-use
# https://stackoverflow.com/questions/50338201/how-to-compress-and-tar-a-folder-in-linux
# https://askubuntu.com/questions/392885/how-can-i-view-the-contents-of-tar-gz-file-without-extracting-from-the-command-l
