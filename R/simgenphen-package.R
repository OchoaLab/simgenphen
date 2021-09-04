#' Simulate Genotypes and Phenotypes
#'
#' The `simgenphen` package is a wrapper around several more specialized packages: `bnpsd` for admixture, `simfam` for family structure, and `simtrait` for complex traits.
#' This package assumes reasonable defaults and generally saves on typing for some recurrent cases.
#' Key non-trivial functionality is drawing genotypes for admixed families whose loci are all polymorphic.
#'
#' @examples
#' ### one liner version
#' # admixture + family genotypes and phenotype
#' # run with default values, except smaller
#' data <- sim_gen_phen( n_ind = 50, G = 3, m_loci = 100, m_causal = 5 )
#'
#' # main objects of interest:
#' # genotype matrix
#' data$X
#' # trait vector
#' data$trait
#'
#' ### version with separate steps
#' # More efficient in certain loops where population structure is constant
#' # but want newly-drawn genotypes and phenotypes.
#' # Also returns additional parameters
#'
#' # first construct population parameters
#' # (relatedness but no loci)
#' k_subpops <- 3
#' data_pop <- sim_pop( n_ind = 50, G = 3, k_subpops = k_subpops )
#' 
#' # now draw genotypes!
#' m_loci <- 100
#' data_geno <- sim_geno(
#'     data_pop$admix_proportions_1,
#'     data_pop$inbr_subpops,
#'     data_pop$fam,
#'     data_pop$ids,
#'     m_loci
#' )
#'
#' # genotype matrix
#' data_geno$X
#' 
#' # now draw trait!
#' m_causal <- 5
#' herit <- 0.5
#' data_trait <- sim_trait_env(
#'     data_geno$X,
#'     data_geno$p_anc,
#'     m_causal,
#'     herit,
#'     env = 'gcat',
#'     k_subpops = k_subpops
#' )
#'
#' # the simulated trait
#' data_trait$trait
#' # indexes of randomly-selected causal loci
#' data_trait$causal_indexes
#' # coefficients of causal loci
#' data_trait$causal_coeffs
#' 
#' @docType package
#' @name simgenphen-package
#' @aliases simgenphen
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
