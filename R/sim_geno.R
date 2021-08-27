sim_geno <- function(
                     admix_proportions,
                     inbr_subpops,
                     fam = NULL,
                     ids = NULL,
                     m_loci = 100000,
                     p_anc_input = NULL,
                     verbose = TRUE
                     ) {
    if ( missing( admix_proportions ) )
        stop( '`admix_proportions` is required!' )
    if ( missing( inbr_subpops ) )
        stop( '`inbr_subpops` is required!' )
    
    # draw allele freqs and genotypes
    if (verbose)
        message('draw_all_admix')
    out <- bnpsd::draw_all_admix(admix_proportions, inbr_subpops, m_loci, p_anc = p_anc_input)
    X <- out$X # genotypes
    p_anc <- out$p_anc # ancestral AFs

    if ( !is.null( fam ) ) {
        # expect family structure now!
        if ( is.null( ids ) )
            stop( '`ids` must be non-`NULL` if `fam` is not `NULL`' )
        G <- length( ids )
    
        # `simfam` requires names for `X`
        colnames( X ) <- ids[[ 1 ]]
        
        if (verbose)
            message('geno_last_gen')
        # final genotypes
        X <- simfam::geno_last_gen(X, fam, ids)

        # NOTE: p_anc doesn't change, this is accurate
        
        # handle fixed loci (a big pain!)
        fixed_loci_indexes <- bnpsd::fixed_loci(X)
        m_loci_fixed <- sum( fixed_loci_indexes )
        while (m_loci_fixed > 0) { # draw things anew, over and over until nothing was fixed
            # draw allele freqs and genotypes
            out <- bnpsd::draw_all_admix(admix_proportions, inbr_subpops, m_loci_fixed)
            # overwrite fixed loci with redrawn polymorphic data
            #            X[fixed_loci_indexes, ] <- out$X # genotypes
            p_anc[fixed_loci_indexes] <- out$p_anc # ancestral AFs
            X_redrawn <- out$X # renamed for clarity
            # `simfam` requires names for `X_redrawn`
            colnames( X_redrawn ) <- ids[[ 1 ]]

            if (verbose)
                message('geno_last_gen (redrawn)')
            # repeat children draws through generations, then
            # overwrite fixed loci with redrawn (hopefully) polymorphic data
            X[fixed_loci_indexes, ] <- simfam::geno_last_gen( X_redrawn, fam, ids )
            
            # look for remaining fixed loci (to continue or stop loop)
            fixed_loci_indexes <- bnpsd::fixed_loci(X)
            m_loci_fixed <- sum( fixed_loci_indexes )
        }
    }

    return(
        list(
            X = X,
            p_anc = p_anc
        )
    )
}

