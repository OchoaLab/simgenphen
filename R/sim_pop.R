sim_pop <- function(
                    n_ind = 1000,
                    k_subpops = 3,
                    bias_coeff = 0.5,
                    fst = 0.3,
                    G = 1,
                    verbose = TRUE
                    ) {
    # FST values for admixture subpopulations
    inbr_subpops <- 1 : k_subpops
    # high admixture FSTs require subpopulations to all be more differentiated, this works when k==3
    if ( fst >= 0.5 )
        inbr_subpops <- inbr_subpops + 1
    
    # define population structure
    if (verbose)
        message('admix_prop_1d_linear')
    obj <- bnpsd::admix_prop_1d_linear(
        n_ind = n_ind,
        k_subpops = k_subpops,
        bias_coeff = bias_coeff,
        coanc_subpops = inbr_subpops,
        fst = fst
    )
    # in this case return value is a named list with three items:
    admix_proportions <- obj$admix_proportions # admixture proportions
    inbr_subpops <- obj$coanc_subpops # rescaled Fst vector for intermediate subpops

    # get true kinship matrix
    if (verbose)
        message('coanc_admix, coanc_to_kinship')
    coancestry <- bnpsd::coanc_admix(admix_proportions, inbr_subpops) # the coancestry matrix
    kinship <- bnpsd::coanc_to_kinship(coancestry) # kinship matrix

    # should be null when there's no family structure
    fam <- NULL
    ids <- NULL
    if (G > 1) {
        # simulate realistic generations of families

        # simulate pedigree first
        if (verbose)
            message('sim_pedigree')
        data_simfam <- simfam::sim_pedigree( n_ind, G )
        fam <- data_simfam$fam # the pedigree itself
        ids <- data_simfam$ids # to filter generations later

        # prune fam table now, to not simulate unnecessary individuals without descendants
        if (verbose)
            message('prune_fam')
        fam <- simfam::prune_fam( fam, ids[[ G ]] )

        if (verbose)
            message('kinship_last_gen')
        # label founders in matrix
        rownames( kinship ) <- ids[[ 1 ]]
        colnames( kinship ) <- ids[[ 1 ]]
        # now calculate total kinship in final generation
        kinship <- simfam::kinship_last_gen( kinship, fam, ids )

        if (verbose)
            message('admix_last_gen')
        # label founders in matrix
        rownames( admix_proportions ) <- ids[[ 1 ]]
        # get correct admixture proportions for these children!
        admix_proportions <- simfam::admix_last_gen( admix_proportions, fam, ids )
    }
    
    return(
        list(
            admix_proportions = admix_proportions,
            inbr_subpops = inbr_subpops,
            fam = fam,
            ids = ids,
            kinship = kinship
        )
    )
}

