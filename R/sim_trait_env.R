#' Simulate a complex trait
#'
#' Mostly a wrapper around [simtrait::sim_trait()], but can optionally add a structured environment effect.
#' However, in that case the heritability is no longer as specified.
#'
#' @param X The genotype matrix.
#' @param p_anc The vector of ancestral allele frequencies.
#' @param k_subpops The number of intermediate subpopulations for admixture model.
#' Used for `env = 'gcat'` only.
#' @inheritParams sim_gen_phen
#'
#' @return A list containing the following elements:
#' - `trait`: The simulated trait vector.
#' - `causal_indexes`: The vector of randomly selected causal loci indexes.
#' - `causal_coeffs`: The vector of simulated regression coefficients for causal loci.
#'
#' @examples
#' # draw population parameters using `sim_pop` first
#' # a small example
#' k_subpops <- 3
#' data <- sim_pop( n_ind = 50, G = 3, k_subpops = k_subpops )
#'
#' # then draw genotypes
#' m_loci <- 100
#' data2 <- sim_geno(
#'     data$admix_proportions_1,
#'     data$inbr_subpops,
#'     data$fam,
#'     data$ids
#' )
#'
#' # now draw trait!
#' m_causal <- 5
#' herit <- 0.5
#' data3 <- sim_trait_env(
#'     data2$X,
#'     data2$p_anc,
#'     m_causal,
#'     herit,
#'     env = 'gcat',
#'     k_subpops = k_subpops
#' )
#'
#' # the simulated trait
#' data3$trait
#' # indexes of randomly-selected causal loci
#' data3$causal_indexes
#' # coefficients of causal loci
#' data3$causal_coeffs
#'
#' @export
sim_trait_env <- function(
                          X,
                          p_anc,
                          m_causal,
                          herit,
                          env = NA,
                          env_var = 1,
                          k_subpops = NA,
                          fes = FALSE,
                          verbose = TRUE
                          ) {

    if (verbose) message('simtrait')

    # parameters of simulation

    # create simulated trait and associated data
    # use version for known p_anc (prefered, only applicable to simulated data)
    # let sim_trait complain if its params are missing
    obj <- simtrait::sim_trait(
        X = X,
        m_causal = m_causal,
        herit = herit,
        p_anc = p_anc,
        fes = fes
    )
    
    if (!is.na(env) && env == 'gcat') {
        if ( is.na( k_subpops ) )
            stop( '`k_subpops` is required if `env = "gcat"`!' )
        
        # construct non-genetic environment effect
        # this is similar, in spirit, to the GCAT paper simulations
        n_ind <- length( obj$trait )
        # in this simulation, cluster membership is given by existing index:
        # the number is used as the environmental effect
        env <- ceiling( ( 1 : n_ind ) / n_ind * k_subpops )
        # rescale with sample variance, to have unit variance
        #        env <- env / sqrt( var(env) )
        # center and scale
        env <- as.numeric( scale( env ) )
        # scale to get desired variance 
        env <- env * sqrt( env_var )
        # simple combine
        obj$trait <- obj$trait + env
    }

    # return what we want from trait simulation
    return(
        list(
            trait = obj$trait,
            causal_indexes = obj$causal_indexes,
            causal_coeffs = obj$causal_coeffs
        )
    )
}
