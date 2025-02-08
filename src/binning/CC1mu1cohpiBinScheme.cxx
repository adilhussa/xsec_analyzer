// XSecAnalyzer includes
#include "XSecAnalyzer/Binning/CC1mu1cohpiBinScheme.hh"

CC1mu1cohpiBinScheme::CC1mu1cohpiBinScheme() : BinSchemeBase( "CC1mu1cohpiBinScheme" ) {}

void CC1mu1cohpiBinScheme::DefineBlocks() {

  /////// Set some standard variables before managing the blocks

  // TTree name for the post-processed ntuples
  ntuple_ttree_name_ = "stv_tree";

  // Run numbers to use when plotting migration matrices
  runs_to_use_ = { 3 };

  // Prefix for the output bin and slice configuration text files
  out_config_prefix_ = "CC1mu1cohpi_";

  // Selection to use with this binning scheme
  selection_name_ = "CC1mu1cohpi";

  // TDirectory file name to use when producing the univmake output histograms
  out_tdir_name_ = "CC1mu1cohpi_coneangle_1D_bin";

  /////// Define the blocks of bins in both true and reco space

  // First block: cos_theta_mu in 1D
  std::vector< double > coneangle_1D_edges = {0.   , 0.02 , 0.041, 0.061, 0.082, 0.102, 0.122, 0.143, 0.163,
      0.184, 0.204, 0.224, 0.245, 0.265, 0.286, 0.306, 0.327, 0.347,
      0.367, 0.388, 0.408, 0.429, 0.449, 0.469, 0.49 , 0.51 , 0.531,
      0.551, 0.571, 0.592, 0.612, 0.633, 0.653, 0.673, 0.694, 0.714,
      0.735, 0.755, 0.776, 0.796, 0.816, 0.837, 0.857, 0.878, 0.898,
      0.918, 0.939, 0.959, 0.98 , 1.};

  Block1D* b1t = new Block1D( "CC1mu1cohpi_true_coneangle",
    "cos(cone angle)", "cos(cone angle)", coneangle_1D_edges,
    "CC1mu1cohpi_Selected", kSignalTrueBin );

  Block1D* b1r = new Block1D( "CC1mu1cohpi_reco_coneangle",
    "cos(cone angle)", "cos(cone angle)", coneangle_1D_edges,
    "!CC1mu1cohpi_MC_Signal", kOrdinaryRecoBin );

  vect_block.emplace_back( b1t, b1r );

  // Second block: p_mu in 1D
  std::vector< double > opnangle_1D_edges = {0.   , 0.01 , 0.02 , 0.03 , 0.04 , 0.051, 0.061, 0.071, 0.081,
      0.091, 0.101, 0.111, 0.121, 0.131, 0.141, 0.152, 0.162, 0.172,
      0.182, 0.192, 0.202, 0.212, 0.222, 0.232, 0.242, 0.253, 0.263,
      0.273, 0.283, 0.293, 0.303, 0.313, 0.323, 0.333, 0.343, 0.354,
      0.364, 0.374, 0.384, 0.394, 0.404, 0.414, 0.424, 0.434, 0.444,
      0.455, 0.465, 0.475, 0.485, 0.495, 0.505, 0.515, 0.525, 0.535,
      0.545, 0.556, 0.566, 0.576, 0.586, 0.596, 0.606, 0.616, 0.626,
      0.636, 0.646, 0.657, 0.667, 0.677, 0.687, 0.697, 0.707, 0.717,
      0.727, 0.737, 0.747, 0.758, 0.768, 0.778, 0.788, 0.798, 0.808,
      0.818, 0.828, 0.838, 0.848, 0.859, 0.869, 0.879, 0.889, 0.899,
      0.909, 0.919, 0.929, 0.939, 0.949, 0.96 , 0.97 , 0.98 , 0.99 ,
      1.};

  Block1D* b2t = new Block1D( "CC1mu1cohpi_true_coneangle",
    "Opening angle", "cos(cone angle)", opnangle_1D_edges,
    "CC1mu1cohpi_MC_Signal", kSignalTrueBin );

  Block1D* b2r = new Block1D( "CC1mu1cohpi_reco_coneangle",
    "cos(cone angle)", "cos(cone angle)", opnangle_1D_edges,
    "CC1mu1cohpi_Selected", kOrdinaryRecoBin );

  vect_block.emplace_back( b2t, b2r );
}
