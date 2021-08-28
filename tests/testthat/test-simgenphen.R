### CHECKERS ###

# write generic checkers that are applied repeatedly

check_kinship <- function( kinship, n_ind ) {
    expect_true( is.matrix( kinship ) )
    expect_equal( nrow( kinship ), n_ind )
    expect_true( isSymmetric( kinship ) )
    expect_true( min( kinship ) >= 0 )
    expect_true( max( kinship ) <= 1 )
}

check_admix_prop <- function( admix, n_ind, k_subpops ) {
    expect_true( is.matrix( admix ) )
    expect_equal( nrow( admix ), n_ind )
    expect_equal( ncol( admix ), k_subpops )
    admix_row_sums <- rowSums( admix ) # calculate
    names( admix_row_sums ) <- NULL # don't want names to be part of comparison
    expect_equal( admix_row_sums, rep.int( 1, n_ind ) )
    expect_true( min( admix ) >= 0 )
    expect_true( max( admix ) <= 1 )
}

check_sim_pop_data <- function ( data, n_ind, k_subpops, fst, G ) {
    # check overall object
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('admix_proportions_1', 'admix_proportions', 'inbr_subpops', 'fam', 'ids', 'kinship') )
    
    # check admixture proportions matrix of founders
    check_admix_prop( data$admix_proportions_1, n_ind, k_subpops )
    
    # check inbr_subpops vector
    inbr_subpops <- data$inbr_subpops
    expect_true( is.numeric( inbr_subpops ) )
    expect_equal( length( inbr_subpops ), k_subpops )
    expect_true( min( inbr_subpops ) >= 0 )
    expect_true( max( inbr_subpops ) <= 1 )
    
    # check kinship of last generation (founders if G=1)
    check_kinship( data$kinship, n_ind )

    # check admixture proportions matrix of last generation (founders if G=1)
    check_admix_prop( data$admix_proportions, n_ind, k_subpops )
    
    # check fst, for admixture params of founders
    admix <- data$admix_proportions_1
    expect_equal( mean( diag( bnpsd::coanc_admix( admix, inbr_subpops ) ) ), fst )
    
    # checks dependent on family structure or not
    if ( G == 1 ) {
        # because there's no family structure here, the last two must be null
        expect_true( is.null( data$fam ) )
        expect_true( is.null( data$ids ) )
    } else {
        # check fam
        fam <- data$fam
        expect_true( is.data.frame( fam ) )
        expect_equal( names( fam ), c('fam', 'id', 'pat', 'mat', 'sex', 'pheno') )
        # number of rows varies because of pruning step, but we do have a clear maximum size
        expect_true( nrow( fam ) <= n_ind * G )
        
        # check ids
        ids <- data$ids
        expect_true( is.list( ids ) )
        expect_equal( length( ids ), G )
    }
}

check_geno <- function( X, n_ind, m_loci ) {
    expect_true( is.matrix( X ) )
    expect_true( all( X %in% c(0L, 1L, 2L) ) )
    expect_equal( nrow( X ), m_loci )
    expect_equal( ncol( X ), n_ind )
    # only special requirement is that no sites are fixed!
    x_bar <- rowMeans( X )
    expect_true( min( x_bar ) > 0 )
    expect_true( max( x_bar ) < 2 )
}

check_sim_geno_data <- function( data, n_ind, m_loci ) {
    # check overall object
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('X', 'p_anc') )

    # check genotypes
    check_geno( data$X, n_ind, m_loci )
    
    # check p_anc
    p_anc <- data$p_anc
    expect_true( is.numeric( p_anc ) )
    expect_equal( length( p_anc ), m_loci )
    expect_true( min( p_anc ) >= 0 )
    expect_true( max( p_anc ) <= 1 )
}

check_sim_bim <- function( bim, m_loci, n_chr ) {
    # validate output
    expect_true( is.data.frame( bim ) )
    expect_equal( names( bim ), c('chr', 'id', 'posg', 'pos', 'ref', 'alt') )
    expect_equal( nrow( bim ), m_loci )
    # validate main non-trivial feature, that there are various chromosomes!
    expect_true( all( bim$chr %in% 1 : n_chr ) )
    # ids are unique
    expect_equal( length( bim$id ), length( unique( bim$id ) ) )
    # no positions are zero or negative
    expect_true( all( bim$pos > 0 ) )
    # the rest are defaults from genio::make_bim, meh
}

check_sim_trait_env_data <- function( data, n_ind, m_loci, m_causal ) {
    # check overall object
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('trait', 'causal_indexes', 'causal_coeffs') )

    # check trait
    trait <- data$trait
    expect_true( is.numeric( trait ) )
    expect_equal( length( trait ), n_ind )

    # check causal_indexes
    indexes <- data$causal_indexes
    expect_true( is.numeric( indexes ) )
    expect_equal( length( indexes ), m_causal )
    expect_true( min( indexes ) >= 1 )
    expect_true( max( indexes ) <= m_loci )

    # check causal_coeffs
    coeffs <- data$causal_coeffs
    expect_true( is.numeric( coeffs ) )
    expect_equal( length( coeffs ), m_causal )
}

### TESTS ###

# now actually run specific tests on toy simulated cases

test_that( "sim_pop works", {
    # there are no mandatory params
    
    # default values, but have them defined here to check all dimensions
    n_ind <- 1000
    k_subpops <- 3
    fst <- 0.3
    G <- 1
    expect_silent( 
        data <- sim_pop(
            n_ind = n_ind,
            k_subpops = k_subpops,
            fst = fst,
            G = G,
            verbose = FALSE
        )
    )
    check_sim_pop_data( data, n_ind, k_subpops, fst, G )

    # try higher K
    n_ind <- 1000
    k_subpops <- 10
    fst <- 0.1 # must lower FST to make it all work
    G <- 1
    expect_silent( 
        data <- sim_pop(
            n_ind = n_ind,
            k_subpops = k_subpops,
            fst = fst,
            G = G,
            verbose = FALSE
        )
    )
    check_sim_pop_data( data, n_ind, k_subpops, fst, G )

    # a special case in the code is for high FST
    # (but lower k)
    n_ind <- 1000
    k_subpops <- 3
    fst <- 0.5
    G <- 1
    expect_silent( 
        data <- sim_pop(
            n_ind = n_ind,
            k_subpops = k_subpops,
            fst = fst,
            G = G,
            verbose = FALSE
        )
    )
    check_sim_pop_data( data, n_ind, k_subpops, fst, G )

    # now add generations!
    # return other params to defaults
    n_ind <- 1000
    k_subpops <- 3
    fst <- 0.3
    G <- 3
    expect_silent( 
        data <- sim_pop(
            n_ind = n_ind,
            k_subpops = k_subpops,
            fst = fst,
            G = G,
            verbose = FALSE
        )
    )
    check_sim_pop_data( data, n_ind, k_subpops, fst, G )
    
})

test_that( "sim_geno works", {
    # always have to use sim_pop first, same as before but much smaller sample sizes
    # will check those outputs again here
    n_ind <- 100
    m_loci <- 50
    k_subpops <- 3
    fst <- 0.3
    G <- 1
    p_anc <- NULL
    expect_silent( 
        data <- sim_pop(
            n_ind = n_ind,
            k_subpops = k_subpops,
            fst = fst,
            G = G,
            verbose = FALSE
        )
    )
    check_sim_pop_data( data, n_ind, k_subpops, fst, G )
    admix_proportions_1 <- data$admix_proportions_1
    inbr_subpops <- data$inbr_subpops
    fam <- data$fam
    ids <- data$ids
    kinship <- data$kinship
    # sim_geno has required params without defaults, check that it dies here only (one time only)
    expect_error( sim_geno() )
    expect_error( sim_geno( admix_proportions_1 = admix_proportions_1 ) )
    expect_error( sim_geno( inbr_subpops = inbr_subpops ) )
    # now a successful run
    expect_silent (
        data <- sim_geno(
            admix_proportions_1,
            inbr_subpops,
            fam = fam,
            ids = ids,
            m_loci = m_loci,
            p_anc = p_anc,
            verbose = FALSE
        )
    )
    check_sim_geno_data( data, n_ind, m_loci )

    # change ancestral allele frequencies only
    # this is the only value we've used in practice
    p_anc <- 0.5
    expect_silent (
        data <- sim_geno(
            admix_proportions_1,
            inbr_subpops,
            fam = fam,
            ids = ids,
            m_loci = m_loci,
            p_anc = p_anc,
            verbose = FALSE
        )
    )
    check_sim_geno_data( data, n_ind, m_loci )

    # all the same except with family structure!
    n_ind <- 100
    m_loci <- 50
    k_subpops <- 3
    fst <- 0.3
    G <- 3
    p_anc <- NULL
    expect_silent( 
        data <- sim_pop(
            n_ind = n_ind,
            k_subpops = k_subpops,
            fst = fst,
            G = G,
            verbose = FALSE
        )
    )
    check_sim_pop_data( data, n_ind, k_subpops, fst, G )
    admix_proportions_1 <- data$admix_proportions_1
    inbr_subpops <- data$inbr_subpops
    fam <- data$fam
    ids <- data$ids
    kinship <- data$kinship
    expect_silent (
        data <- sim_geno(
            admix_proportions_1,
            inbr_subpops,
            fam = fam,
            ids = ids,
            m_loci = m_loci,
            p_anc = p_anc,
            verbose = FALSE
        )
    )
    check_sim_geno_data( data, n_ind, m_loci )

    # change ancestral allele frequencies only
    # this is the only value we've used in practice
    p_anc <- 0.5
    expect_silent (
        data <- sim_geno(
            admix_proportions_1,
            inbr_subpops,
            fam = fam,
            ids = ids,
            m_loci = m_loci,
            p_anc = p_anc,
            verbose = FALSE
        )
    )
    check_sim_geno_data( data, n_ind, m_loci )
})

test_that( "sim_bim works", {
    # there are no mandatory parameters
    # pass params for checks
    m_loci <- 100
    n_chr <- 22
    expect_silent( bim <- sim_bim( m_loci, n_chr ) )
    # validate output
    check_sim_bim( bim, m_loci, n_chr )
})

test_that( "sim_trait_env works", {
    # simulate a small X, same as before, but no validations
    # no particular need for family structure or other non-default choices, so keep it simple
    n_ind <- 100
    m_loci <- 50
    k_subpops <- 3
    expect_silent( 
        data <- sim_pop(
            n_ind = n_ind,
            k_subpops = k_subpops,
            verbose = FALSE
        )
    )
    expect_silent (
        data <- sim_geno(
            data$admix_proportions_1,
            data$inbr_subpops,
            m_loci = m_loci,
            verbose = FALSE
        )
    )
    X <- data$X
    p_anc <- data$p_anc
    
    # set other params
    m_causal <- 5
    herit <- 0.5
    env <- NA
    env_var <- 1
    fes <- FALSE
    expect_silent( 
        data <- sim_trait_env(
            X = X,
            p_anc = p_anc,
            m_causal = m_causal,
            herit = herit,
            env = env,
            env_var = env_var,
            fes = fes,
            verbose = FALSE
        )
    )
    check_sim_trait_env_data( data, n_ind, m_loci, m_causal )

    # test FES version
    fes <- TRUE
    expect_silent( 
        data <- sim_trait_env(
            X = X,
            p_anc = p_anc,
            m_causal = m_causal,
            herit = herit,
            env = env,
            env_var = env_var,
            fes = fes,
            verbose = FALSE
        )
    )
    check_sim_trait_env_data( data, n_ind, m_loci, m_causal )

    # test env version, only gcat is implemented
    env <- 'gcat'
    # revert FES (shouldn't matter either way)
    fes <- FALSE
    expect_silent( 
        data <- sim_trait_env(
            X = X,
            p_anc = p_anc,
            m_causal = m_causal,
            herit = herit,
            env = env,
            env_var = env_var,
            k_subpops = k_subpops, # required for `env='gcat'` only
            fes = fes,
            verbose = FALSE
        )
    )
    check_sim_trait_env_data( data, n_ind, m_loci, m_causal )
})

test_that( "sim_gen_phen works", {
    # no parameters are strictly required
    # test a toy case
    n_ind <- 100
    m_loci <- 50
    k_subpops <- 3
    m_causal <- 5
    n_chr <- 22
    # run!
    expect_silent(
        data <- sim_gen_phen(
            n_ind = n_ind,
            m_loci = m_loci,
            k_subpops = k_subpops,
            m_causal = m_causal,
            n_chr = n_chr,
            verbose = FALSE
        )
    )
    # test overall object
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('X', 'bim', 'kinship', 'admix_proportions', 'trait', 'causal_indexes', 'causal_coeffs') )
    # check X
    check_geno( data$X, n_ind, m_loci )
    # check bim
    check_sim_bim( data$bim, m_loci, n_chr )
    # check kinship
    check_kinship( data$kinship, n_ind )
    # check admixture proportions matrix
    check_admix_prop( data$admix_proportions, n_ind, k_subpops )
    # check trait components, simple subset of given data
    check_sim_trait_env_data(
        data[ c('trait', 'causal_indexes', 'causal_coeffs') ],
        n_ind, m_loci, m_causal
    )
})
