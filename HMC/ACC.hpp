// ACC.hpp — error-estimation routines for spatial autocovariance fields G(x,t).
//
// Naming conventions used throughout this file (and the autocovariance_*.cc
// drivers that include it):
//
//   * "binning"  refers to MD-chain partition into n_bin segments.
//   * "blocking" / "sparsening" refer to SPATIAL-LATTICE coarsening.
//   * "ACC"      = autocorrelation coefficient ρ(t) = G(t)/G(0).
//   * "MFCOV"    = the unnormalized autocovariance G(t) itself.
//   * G_cent / G_conn distinguish the centered vs connected forms of G;
//     the choice is made in the driver, not here. Routines below are
//     agnostic and apply to either form.
//
// Four error estimators are exposed; for the full discussion see §4.2 of
// Master_Field_Type_Autocorrelation/main.tex.
//
//   (1) MF_approx          — Master-Field, single chain. Direct integration of
//                            the spatial covariance density Cov[G(x,t),G(y,t)]
//                            over |x−y| ≤ R. Saturation at R_sat is empirically
//                            observable; l_B is data compression only.
//
//   (2) MS_approx          — Local Madras-Sokal, single chain. Per-site MS
//                            variance from a four-point sum truncated at
//                            inner cutoff W (Lüscher 2005 Eq. E.11; W ≥ 100).
//                            Combining per-site estimates into a spatial-
//                            average error requires assumed site independence,
//                            which is data-circular. Sparse spatial input
//                            only — block-avg MS would need the unknown
//                            intra-block cross-site covariance.
//
//   (3) Block-first per-block binning ("Block-First") — multi-bin estimator
//                            that combines per-cell errors under inter-cell
//                            independence. Same circularity as (2) at the
//                            combination step. Not implemented as a separate
//                            routine here; legacy concept retained for context.
//
//   (4) binning_avg_cov / binning_avg_rho — bin-first per-bin spatial avg.
//                            PREFERRED. Spatial average inside each bin
//                            before inter-bin variance — no spatial-
//                            independence assumption anywhere. Block-avg
//                            spatial input is the safer default within each
//                            bin (sparse loses statistics at large l_B).
//
// The two binning sub-variants of method (4) are:
//
//   binning_avg_cov  — pool covariances first: G_pool(t) = ⟨G_b(t)⟩_b.
//                      Then ρ̂ = G_pool(t)/G_pool(0); variance via delta-
//                      method error propagation using inter-bin
//                      Cov(G_b(t), G_b(0)).
//
//   binning_avg_rho  — per-bin ratio first: ρ_b = G_b(t)/G_b(0).
//                      Then ρ̂ = ⟨ρ_b⟩_b; variance from inter-bin sample
//                      variance of ρ_b. PREFERRED — variance computed
//                      directly from a sample, no linearisation. The two
//                      sub-variants agree at leading order in the delta-
//                      method linearisation but are NOT mathematically
//                      equivalent at finite n_bin; disagreement at the
//                      working n_bin is itself a diagnostic on the
//                      linearisation.
//
// A separate, prior design choice operative throughout: we use the
// autocovariance G_x(t) as the basic statistical variable, not the per-site
// ACC ρ_x = G_x(t)/G_x(0) — at single sites, G_x(0) is a small, fluctuating
// denominator that ruins the signal-to-noise of ρ_x. Spatial averages are
// always taken on G first ("ratio of averages" at the spatial level).

namespace Grid{
  struct WFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(WFParameters,
	    int, tau,
	    std::string, data_name,
	    std::string, path);


    template <class ReaderClass >
    WFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "WilsonFlow", *this);
    }

  };

  struct ConfParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConfParameters,
           std::string, conf_prefix,
	   int, StartConfiguration,
	   int, EndConfiguration);

    template <class ReaderClass >
    ConfParameters(Reader<ReaderClass>& Reader){
      read(Reader, "Configurations", *this);
    }

  };

  struct ACFParameters: Serializable {
    GRID_SERIALIZABLE_CLASS_MEMBERS(ACFParameters,
	   int, MDtime_div_fac,
	   int, MScut,	          // Four-point inner cutoff W (= Λ in Lüscher 2005 Eq. E.11) for MS_approx; W ≥ 100
	   std::vector<int>, space_block_sizes,
	   int, R,                // Summation radius of Master field tecnnique
	   int, isFullTimeAvg);

    template <class ReaderClass >
    ACFParameters(Reader<ReaderClass>& Reader){
      read(Reader, "Autocorrelations", *this);
    }

  };
}

template <class T> void readFile(T& out, std::string const fname){
  Grid::emptyUserRecord record;
  Grid::ScidacReader RD;
  RD.open(fname);
  RD.readScidacFieldRecord(out,record);
  RD.close();
}

// MS_approx — Madras-Sokal variance formula at fixed lag t, applied per-site.
// Implements the four-point variance estimator with inner cutoff W (= Λ in
// Lüscher 2005 Eq. E.11; W ≥ 100 throughout). For the volume-summed scalar
// or for G_x(t) at any single fixed site this is the standard, theoretically
// clean estimator. Combining per-site results into an error on the spatial
// average requires an assumption of site independence not ensured by the
// data — validating it via Bienaymé scaling σ_ρ ∝ l_B^{d/2} is circular.
// Sparse spatial input is the only acceptable choice (block-avg would need
// the unknown intra-block cross-site covariance). Retained as a cross-check.
template <class L, typename A> void MS_approx(Grid::GridBase *Coarse, std::vector<L,A>  const& in, int T, int W, int block_size, int tau,
					      std::string const e_name, std::string ACC_Type ){
  using namespace Grid;
  // applicable only if blocking the src fields
  // Otherwise, for block averaged G_B(t), its error is hard to estimate as we need to consider correlation of G_x(t) at diff. sites wihting a block

  RealD sumG[in.size()], ACC;
  for (int t=0; t<T; t++) {
    sumG[t] = TensorRemove(sum(in[t])).real();
    std::cout << GridLogMessage << ACC_Type + " MFACC " + e_name + " (Madras-Sokal Approx): " << tau << " " << block_size << " " << 1 << " " << t << " "
              << sumG[t]/sumG[0] << std::endl;
    std::cout << GridLogMessage << ACC_Type + " MFCOV "  + e_name + " (Madras-Sokal Approx): " << tau << " " << block_size << " " << 1 << " " << t << " "
	      << sumG[t] << std::endl;
  }

  L var(Coarse), tmp(Coarse), Gt(Coarse), G0(Coarse);
  for (int t=0; 2*t<T-W; t++){ // \sigma_{\rho}(t) uses \rho(t') up to 2t+W
    ACC = sumG[t]/sumG[0];
    var = (1/sumG[t])*in[t] - (1/sumG[0])*in[0];
    var = (2.0/RealD(T))*var*var;
    for (int k=1; k<W+t; k++) {
      tmp = (1.0/sumG[t]/sumG[t]+2.0/sumG[0]/sumG[0])*in[k]*in[k] +
	in[std::abs(k-t)]*( (1.0/sumG[t]/sumG[t])*in[k+t] - (2.0/sumG[0]/sumG[t])*in[k]) -(2.0/sumG[0]/sumG[t])*in[k]*in[k+t];
      var = var + (2.0/RealD(T))*tmp;
    }
    std::cout << GridLogMessage << ACC_Type + " Variance " + e_name + " (Madras-Sokal Approx): " << tau << " " << block_size << " " << 1 << " " << t << " "
	      << TensorRemove(sum(var)).real()*ACC*ACC << std::endl;
  }
}

// binning_avg_cov  (formerly "binning") — average G across MD-time bins,
// then take the ratio.  Variance of the ratio via delta-method error
// propagation using inter-bin Cov(G_b(t), G_b(0)).
//
// Two distinct point estimators of ρ(t).  binning_avg_cov and binning_avg_rho
// agree at leading order in the delta-method linearisation but are NOT
// mathematically equivalent at finite n_bin.  Used as a cross-check on the
// linearisation; disagreement at the working n_bin is itself a diagnostic.
// See the file header for the recommended choice (binning_avg_rho preferred).
template <class L, typename A> void binning_avg_cov(Grid::GridBase *Coarse, std::vector<L,A>  const& in, int T, int n_bin, int block_size, int tau,
						    std::string const e_name, std::string ACC_Type){
  using namespace Grid;
  // error of G_b(t) is first computed via binning for each block b
  // G_b(t) with diff. b is considered as indep. measurement of G(t) => sigma_{G_(t)}^2 is reduced by V_b for the avg over blocks
  // Then, proceed to compute error of G(t)/G(0)
  // Computing ACC_b(t) = G_b(t)/G_b(0),

  L avg(Coarse), avg0(Coarse), var(Coarse), var0(Coarse), cov(Coarse), tmp(Coarse);
  for (int t=0; t<T; t++){

    avg = Zero(); var = Zero(); cov = Zero();

    // Find avg field over binns
    for (int b=0; b<n_bin; b++)
      avg = avg + (1/RealD(n_bin))*in[b*T + t];
    if ( t== 0 ) avg0 = avg;
    RealD ACC = TensorRemove(sum(avg)).real()/TensorRemove(sum(avg0)).real();
    std::cout << GridLogMessage << ACC_Type + " MFACC " + e_name + " (binning_avg_cov): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << ACC << std::endl;

    // Find variance of Gt/G0
    for	(int b=0; b<n_bin; b++){
      tmp = in[b*T + t] - avg;
      var = var + (1/RealD((n_bin-1)*n_bin))*(tmp*tmp);
      cov = cov + (1/RealD((n_bin-1)*n_bin))*tmp*(in[b*T] - avg0);
    }
    RealD G0 = real(sum(avg0)), Gt = real(sum(avg)); // the factor of 1/N_B is cancel out by 1/V^2 from Cov[G0,Gt] in the expression of var of Gt/G0
    if ( t == 0 ) var0 = var;
    std::cout << GridLogMessage << ACC_Type + " Variance " + e_name + " (binning_avg_cov): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << (real(sum(var))/Gt/Gt + real(sum(var0))/G0/G0 - 2.0*real(sum(cov))/G0/Gt)*ACC*ACC
	      << " " << (real(sum(var))/Gt/Gt + real(sum(var0))/G0/G0)*ACC*ACC << " " << - 2.0*real(sum(cov))/G0/Gt*ACC*ACC << std::endl;

  }
}

// binning_avg_rho  (formerly "binning2") — per-bin ratio ρ_b = G_b(t)/G_b(0),
// then average across bins.  Variance from inter-bin sample variance of ρ_b.
//
// PREFERRED.  More transparent (variance from direct sample, no delta-method
// linearization), more robust if per-bin G_b are non-Gaussian, and treats
// each MD-time bin as an independent replicate of the ACC estimate — which
// is the standard interpretation justified by the law of large numbers when
// bin width > τ_int.  Requires n_bin enough for sample variance to converge
// (≳ 20–50 in practice).
template <class L, typename A> void binning_avg_rho(Grid::GridBase *Coarse, std::vector<L,A>  const& in, int T, int n_bin, int block_size, int tau,
						    std::string const e_name, std::string ACC_Type){
  using namespace Grid;
  // if bin_size > 30, the mean over binns approx dist. like Gaussian; the validity of this  estimation rests on Central Limit Theorem
  // we avg G_x(t) over blocks first and then take the ratio G(t)/G(0) before avg over bins
  // block avg might give a better sig as it incorpolates more stats; both should be unbiased

  RealD avg, avg0, var, tmp;
  std::vector<RealD> in_sum(in.size());

  for (int i=0; i<in.size(); i++) {
    in_sum[i] = TensorRemove(sum(in[i])).real();
    std::cout << GridLogMessage << ACC_Type + " MFCOV " + e_name + " (binning_avg_rho): " << tau << " " << block_size << " " << n_bin << " " << i << " "
	      << in_sum[i] << std::endl;
  }

  // Find avg field over binns
  for (int t=0; t<T; t++){
    avg = 0.0; var = 0.0;
    for (int b=0; b<n_bin; b++)
      avg += (1/RealD(n_bin))*in_sum[b*T + t]/in_sum[b*T];
    std::cout << GridLogMessage << ACC_Type + " MFACC " + e_name + " (binning_avg_rho): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << avg << std::endl;

    // Find variance of Gt/G0
    for (int b=0; b<n_bin; b++){
      tmp = in_sum[b*T + t]/in_sum[b*T] - avg;
      var = var + (1/RealD((n_bin-1)*n_bin))*(tmp*tmp);
    }
    std::cout << GridLogMessage << ACC_Type + " Variance " + e_name + " (binning_avg_rho): " << tau << " " << block_size << " " << n_bin << " " << t << " "
              << var << " " << (1-avg*avg)*(1-avg*avg)/RealD(Coarse->gSites()-3)*RealD(n_bin) << std::endl;
  }
}

template <class L, typename A> void MF_approx(Grid::GridBase *Coarse, std::vector<L,A>  const& in, int T, int R, int block_size, int tau,
                                              std::string const e_name, std::string ACC_Type ){
  using namespace Grid;

  // Assume: block averaging

  int R_b = R/block_size;
  RealD sumG[in.size()], ACC;
  for (int t=0; t<T; t++) {
    sumG[t] = TensorRemove(sum(in[t])).real();
    std::cout << GridLogMessage << ACC_Type + " MFACC " + e_name + " (Master-Field Approx): " << tau << " " << block_size << " " << R << " " << t << " "
              << sumG[t]/sumG[0] << std::endl;
  }

  RealD var, cov;
  L tmp(Coarse), Gt(Coarse), G0(Coarse), one(Coarse); one = typename L::scalar_type(1.0,0.0);
  G0 = in[0];
  G0 = G0 - (sumG[0]/RealD(Coarse->gSites())) * one;
  //std::cout << GridLogMessage << "in MF" << T <<  " " << block_size << " "<< R_b << " " << R << std::endl;
  for(int t=0; t<T; t++) {
    Gt = in[t] - (sumG[t]/RealD(Coarse->gSites())) * one;
    tmp = Gt;
    var = TensorRemove(sum(Gt*Gt)).real(); // C_tt
    cov = TensorRemove(sum(Gt*G0)).real(); // C_t0

    for (int r=1; r<=R_b; r++){
      // Shift Gt by y s.t. |y| <= r
      for(int x=-r; x<=r; x++)
	for(int y=-r; y<=r; y++)
	  for(int z=-r; z<=r; z++)
	    for(int s=-r; s<=r; s++) {
	      int r2_c = x*x+y*y+z*z+s*s;
	      if ( (r-1)*(r-1) < r2_c && r2_c <= r*r) {

		int d[4] = {x,y,z,s};
		Gt = tmp;
		for(int mu=0; mu<Nd; mu++)
		  Gt = Cshift(Gt, mu, d[mu]);
		var += TensorRemove(sum(Gt*tmp)).real();
		cov += TensorRemove(sum(Gt*G0)).real();
	      }
	    }

      // Note: No division by the volume of the coarse lattice
      std::cout << GridLogMessage << ACC_Type + " Variance " + e_name + " (Master-Field Approx): " << tau << " " << block_size << " " << r << " " << t << " "
		<< sumG[t] << " " << var << " " << cov << std::endl;
    }
  }
}
