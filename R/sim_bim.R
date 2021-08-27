# construct reasonable BIM table for BOLT, which relies on this data more than other approaches
# main thing is to have more than one chr
sim_bim <- function(
                    m_loci = 100000,
                    n_chr = 22
                    ) {
    bim <- tibble::tibble(
        # chromosomes are all roughly equal length in this case
        chr = ceiling( (1 : m_loci ) / m_loci * n_chr ),
        pos = 0 # initialize for now
    )
    # add positions in this loop, contiguous per chr, restarting after each case
    for ( chr in 1 : n_chr ) {
        indexes <- which( bim$chr == chr )
        bim$pos[ indexes ] <- ( 1 : length( indexes ) ) * 1000 # space them a lot
    }
    # autocomplete the rest
    bim <- genio::make_bim( bim )
    # return that tibble
    return( bim )
}
