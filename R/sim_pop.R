#' Simulate population structure
#'
#' This function defines admixture and optionally family parameters establishing how individuals are related to each other.
#' Admixture proportions are deterministically constructed from the linear 1D admixture model [bnpsd::admix_prop_1d_linear()].
#' If `fst < 0.5` the intermediate inbreeding coefficients are proportional to their index (`1:k_subpops`), otherwise they are proportional to `(1:k_subpops) + 1`.
#' If `G > 1`, a random family structure is drawn using [simfam::sim_pedigree()].
#' No per-locus parameters are constructed or simulated with this function.
#'
#' @inheritParams sim_gen_phen
#'
#' @return A list containing the following elements:
#' - `admix_proportions_1`: The true admixture proportions of the founders (required to draw genotypes for admixed founders and seed the family structure if `G > 1`)
#' - `admix_proportions`: The true admixture proportions of the joint admixed family model (final generation only).
#' - `inbr_subpops`: The vector of intermediate subpopulation inbreeding/FST values for admixture model.
#' - `fam`: If `G > 1`, the pedigree structure as a plink FAM table, pruned with [simfam::prune_fam()] to include only ancestors of the final generation; otherwise `NULL`.
#' - `ids`: If `G > 1`, a list of IDs from `fam` above split into non-overlapping generations; otherwise `NULL`.
#' - `kinship`: The true kinship of the joint admixed family model (final generation only).
#'
#' @examples
#' # a small example
#' data <- sim_pop( n_ind = 50, G = 3 )
#' 
#' # parameters for admixture model of founders
#' data$admix_proportions_1
#' data$inbr_subpops
#' 
#' # parameters of family structure
#' data$fam
#' data$ids
#' 
#' # parameters for final generation
#' data$admix_proportions
#' data$kinship
#' 
#' @seealso
#' The `bnpsd` and `simfam` packages.
#'
#' @export
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
    admix_proportions_1 <- obj$admix_proportions # admixture proportions
    inbr_subpops <- obj$coanc_subpops # rescaled Fst vector for intermediate subpops

    # get true kinship matrix
    if (verbose)
        message('coanc_admix, coanc_to_kinship')
    coancestry <- bnpsd::coanc_admix( admix_proportions_1, inbr_subpops ) # the coancestry matrix
    kinship <- bnpsd::coanc_to_kinship( coancestry ) # kinship matrix
    
    # should be null when there's no family structure
    fam <- NULL
    ids <- NULL
    # let these be the final ones either way (if G == 1 then use *_1)
    admix_proportions_G <- admix_proportions_1
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
        rownames( admix_proportions_1 ) <- ids[[ 1 ]]
        # get correct admixture proportions for these children!
        admix_proportions_G <- simfam::admix_last_gen( admix_proportions_1, fam, ids )
    }
    
    return(
        list(
            admix_proportions_1 = admix_proportions_1,
            admix_proportions = admix_proportions_G,
            inbr_subpops = inbr_subpops,
            fam = fam,
            ids = ids,
            kinship = kinship
        )
    )
}

