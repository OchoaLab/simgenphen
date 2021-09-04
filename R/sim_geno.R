#' Simulate genotypes from an admixed family structure
#'
#' This function draws genotypes from an admixture model and optionally propagates them through a pedigree.
#' The admixture genotypes are calculated using [bnpsd::draw_all_admix()] and simply returned if `fam` is `NULL`.
#' If `fam` is not `NULL`, genotypes are propagated through the pedigree using [simfam::geno_last_gen()].
#' The admixture model by default redraws fixed loci from the model.
#' However, loci may become fixed after propagation through the pedigree.
#' This function identifies these fixed loci (using [bnpsd::fixed_loci()]) and redraws them first from the admixture model and then through the pedigree, iterating until no loci are fixed.
#'
#' @param admix_proportions_1 The admixture proportions matrix of the founders (required to draw the founder genotypes from the admixture model).
#' @param inbr_subpops The vector of intermediate subpopulation inbreeding/FST values for admixture model.
#' @param fam The pedigree structure as a plink FAM table, or `NULL` if no family structure is to be simulated.
#' @param ids The list of IDs from `fam` above split into non-overlapping generations.
#' Ignored if `fam` is `NULL`.
#' @param p_anc The desired ancestral allele frequencies (scalar or length-`m_loci` vector passed to [bnpsd::draw_all_admix()]).
#' By default (`NULL`), ancestral allele frequencies are drawn randomly.
#' @inheritParams sim_gen_phen
#'
#' @return A list containing the following elements:
#' - `X`: The simulated genotype matrix.
#' - `p_anc`: The vector of ancestral allele frequencies (required to simulate traits with desired heritability).
#'
#' @examples
#' # draw population parameters using `sim_pop` first
#' # a small example
#' data_pop <- sim_pop( n_ind = 50, G = 3 )
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
#' # ancestral allele frequencies
#' data_geno$p_anc
#' 
#' @export
sim_geno <- function(
                     admix_proportions_1,
                     inbr_subpops,
                     fam = NULL,
                     ids = NULL,
                     m_loci = 100000,
                     p_anc = NULL,
                     verbose = TRUE
                     ) {
    if ( missing( admix_proportions_1 ) )
        stop( '`admix_proportions_1` is required!' )
    if ( missing( inbr_subpops ) )
        stop( '`inbr_subpops` is required!' )

    # rename internally to avoid confusion with output p_anc, which gets edited for fixed loci (so we need to refer to input again in those cases)
    p_anc_input <- p_anc
    
    # draw allele freqs and genotypes
    if (verbose)
        message('draw_all_admix')
    out <- bnpsd::draw_all_admix(admix_proportions_1, inbr_subpops, m_loci, p_anc = p_anc_input)
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
            out <- bnpsd::draw_all_admix(admix_proportions_1, inbr_subpops, m_loci_fixed)
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

