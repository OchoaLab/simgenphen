#' Construct a variant table in plink BIM format with several chromosomes.
#'
#' Purpose is to build an artificial BIM table with several chromosomes, which some methods require, such as leave-one-chromosome-out (LOCO) estimation.
#' This function produces chromosomes with approximately equal numbers of variants.
#' Positions are spaced uniformly by a given size.
#' 
#' @param pos_gap The number of basepairs to space variants.
#' @inheritParams sim_gen_phen
#'
#' @return A `tibble` with these columns:
#' - `chr`: The constructed chromosome positions.
#' - `id`: Locus ID values (autocompleted by [genio::make_bim()] to match row numbers).
#' - `posg`: Genetic distance (autocompleted by [genio::make_bim()] to `0`, all missing).
#' - `pos`: locus positions, starting at `pos_gap` for each chromosome and increasing in increments of `pos_gap`.
#' - `ref`: Reference allele (autocompleted by [genio::make_bim()] to `1`).
#' - `alt`: Alternative allele (autocompleted by [genio::make_bim()] to `2`).
#'
#' @examples
#' # simulate a small table
#' m_loci <- 100
#' bim <- sim_bim( m_loci )
#' 
#' @export
sim_bim <- function(
                    m_loci = 100000,
                    n_chr = 22,
                    pos_gap = 1000
                    ) {
    bim <- tibble::tibble(
        # chromosomes are all roughly equal length in this case
        chr = ceiling( (1 : m_loci ) / m_loci * n_chr ),
        pos = 0 # initialize for now
    )
    # add positions in this loop, contiguous per chr, restarting after each case
    for ( chr in 1 : n_chr ) {
        indexes <- which( bim$chr == chr )
        bim$pos[ indexes ] <- ( 1 : length( indexes ) ) * pos_gap # space them
    }
    # autocomplete the rest
    bim <- genio::make_bim( bim )
    # return that tibble
    return( bim )
}
