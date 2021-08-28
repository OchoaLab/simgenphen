#' Simulate genotypes and phenotypes
#'
#' Simulate entire parameters and random variables, particularly random genotypes and phenotypes, with desired dimensions and other base parameters.
#'
#' @param n_ind The number of individuals.
#' @param m_loci The number of loci.
#' @param k_subpops The number of intermediate subpopulations for admixture model.
#' @param fst The FST of the admixed individuals (for the founders if `G > 1`).
#' @param bias_coeff The bias coefficient of the admixed individuals (for the founders if `G > 1`).
#' @param G The number of generations for random family.
#' The `G = 1` case is no family structure (just admixture).
#' For `G > 1` first admixed founders are simulated, then the family structure is simulated from them.
#' @param m_causal The number of causal loci for the trait, selected randomly from among the simulated loci.
#' @param herit The trait heritability.
#' @param env A string describing environment model.
#' Only `NA` (no environment) or "gcta" are accepted.
#' @param env_var The variance of the environment effect.
#' Ignored if `env` is `NA`.
#' @param fes If `FALSE` (default) constructs trait from Random Coefficients (RC) model.
#' If `TRUE`, the Fixed Effect Sizes (FES) model is used instead.
#' @param n_chr Number of chromosomes to simulate.
#' Chromosome assignments are not biologically meaningful, as all loci are drawn independently (no LD).
#' @param verbose If `TRUE` reports progress, otherwise it is silent.
#'
#' @return A list containing the following elements:
#' - `X`: The simulated genotype matrix.
#' - `bim`: The variant info table.
#' - `kinship`: The true kinship of the joint admixed family model (final generation only).
#' - `admix_proportions`: The true admixture proportions of the joint admixed family model (final generation only).
#' - `trait`: The simulated trait vector.
#' - `causal_indexes`: The vector of randomly selected causal loci indexes.
#' - `causal_coeffs`: The vector of simulated regression coefficients for causal loci.
#'
#' @examples
#' # run with default values, except smaller
#' data <- sim_gen_phen( n_ind = 50, G = 3, m_loci = 100, m_causal = 5 )
#'
#' # main objects of interest:
#' # genotype matrix
#' data$X
#' # trait vector
#' data$trait
#'
#' @seealso
#' This function is a wrapper around [sim_pop()], [sim_geno()], [sim_bim()], and [sim_trait_env()], see those for more details.
#'
#' @export
sim_gen_phen <- function(
                         # params of population structure
                         n_ind = 1000, # number of individuals (columns)
                         m_loci = 100000, # number of loci (rows)
                         k_subpops = 3, # number of intermediate subpopulations
                         fst = 0.3, # desired final Fst
                         bias_coeff = 0.5, # bias coeff of standard Fst estimator
                         # the family code gets triggered if G > 1
                         G = 1, # number of generations (1 is unrelated founders)
                         # params of complex trait
                         m_causal = 100,
                         herit = 0.8,
                         env = NA,
                         env_var = 1,
                         fes = FALSE,
                         n_chr = 22,
                         verbose = TRUE # to print messages
                         ) {
    # this function is for generating random genotypes from an admixed population with 3 source subpopulations
    # meant to abstractly resemble Hispanics and African-Americans

    # things that can't be changed:
    # - Fst of intermediate subpopulations is a ramp proportional to 1:k
    # - Admixture proportions are from 1D geography model

    ##################
    ### 0: SIM POP ###
    ##################

    obj <- sim_pop(
        n_ind = n_ind,
        k_subpops = k_subpops,
        bias_coeff = bias_coeff,
        fst = fst,
        G = G,
        verbose = verbose
    )
    admix_proportions_1 <- obj$admix_proportions_1
    admix_proportions <- obj$admix_proportions
    inbr_subpops <- obj$inbr_subpops
    fam <- obj$fam
    ids <- obj$ids
    kinship <- obj$kinship
    
    ####################
    ### 1: GENOTYPES ###
    ####################

    obj <- sim_geno(
        admix_proportions_1 = admix_proportions_1,
        inbr_subpops = inbr_subpops,
        fam = fam,
        ids = ids,
        m_loci = m_loci,
        verbose = verbose
    )
    X <- obj$X
    p_anc <- obj$p_anc

    # BOLT requires more than one chr, this aims for human chrs (vaguely)
    bim <- sim_bim(
        m_loci = m_loci,
        n_chr = n_chr
    )
    
    ################
    ### 2: TRAIT ###
    ################

    obj <- sim_trait_env(
        X = X,
        p_anc = p_anc,
        m_causal = m_causal,
        herit = herit,
        env = env,
        env_var = env_var,
        k_subpops = k_subpops,
        fes = fes,
        verbose = verbose
    )
    trait <- obj$trait
    causal_indexes <- obj$causal_indexes
    causal_coeffs <- obj$causal_coeffs

    #################
    ### 3: RETURN ###
    #################

    # return a few items of interest
    # NOTES:
    # - inbr_subpops, p_anc also not used outside
    # - could return fam/ids, don't right now
    return(
        list(
            X = X, # genotype matrix
            bim = bim,
            kinship = kinship, # kinship matrix of last generation (founders if G=1)
            admix_proportions = admix_proportions, # admixture proportions of last generation (founders if G=1)
            trait = trait, # trait vector
            causal_indexes = causal_indexes, # randomly-picked causal locus index
            causal_coeffs = causal_coeffs # locus effect size vector (for causal loci only, rest are zero)
        )
    )
}
