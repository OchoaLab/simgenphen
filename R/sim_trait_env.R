sim_trait_env <- function(
                          X,
                          p_anc,
                          m_causal,
                          herit,
                          env = NA,
                          env_var = 1,
                          k_subpops = NA, # used for env='gcat'
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
