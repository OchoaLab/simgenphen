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
    admix_proportions <- obj$admix_proportions
    inbr_subpops <- obj$inbr_subpops
    fam <- obj$fam
    ids <- obj$ids
    kinship <- obj$kinship
    
    ####################
    ### 1: GENOTYPES ###
    ####################

    obj <- sim_geno(
        admix_proportions = admix_proportions,
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
    # - we never use coancestry outside, only kinship
    # - inbr_subpops, p_anc also not used outside
    # - could return fam/ids, don't right now
    return(
        list(
            X = X, # genotype matrix
            bim = bim,
            kinship = kinship, # kinship matrix
            admix_proportions = admix_proportions, # admixture proportions
            trait = trait, # trait vector
            causal_indexes = causal_indexes, # randomly-picked causal locus index
            causal_coeffs = causal_coeffs # locus effect size vector (for causal loci only, rest are zero)
        )
    )
}
