% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_gen_phen.R
\name{sim_gen_phen}
\alias{sim_gen_phen}
\title{Simulate genotypes and phenotypes}
\usage{
sim_gen_phen(
  n_ind = 1000,
  m_loci = 1e+05,
  k_subpops = 3,
  fst = 0.3,
  bias_coeff = 0.5,
  G = 1,
  m_causal = 100,
  herit = 0.8,
  env = NA,
  env_var = 1,
  fes = FALSE,
  n_chr = 22,
  beta = NA,
  verbose = TRUE
)
}
\arguments{
\item{n_ind}{The number of individuals.}

\item{m_loci}{The number of loci.}

\item{k_subpops}{The number of intermediate subpopulations for admixture model.}

\item{fst}{The FST of the admixed individuals (for the founders if \code{G > 1}).}

\item{bias_coeff}{The bias coefficient of the admixed individuals (for the founders if \code{G > 1}).}

\item{G}{The number of generations for random family.
The \code{G = 1} case is no family structure (just admixture).
For \code{G > 1} first admixed founders are simulated, then the family structure is simulated from them.}

\item{m_causal}{The number of causal loci for the trait, selected randomly from among the simulated loci.}

\item{herit}{The trait heritability.}

\item{env}{A string describing environment model.
Only \code{NA} (no environment) or "gcta" are accepted.}

\item{env_var}{The variance of the environment effect.
Ignored if \code{env} is \code{NA}.}

\item{fes}{If \code{FALSE} (default) constructs trait from Random Coefficients (RC) model.
If \code{TRUE}, the Fixed Effect Sizes (FES) model is used instead.}

\item{n_chr}{Number of chromosomes to simulate.
Chromosome assignments are not biologically meaningful, as all loci are drawn independently (no LD).}

\item{beta}{Shape parameter for a symmetric Beta for ancestral allele frequencies \code{p_anc}.
If \code{NA} (default), \code{p_anc} is uniform with range in [0.01, 0.5].
Otherwise, \code{p_anc} has a symmetric Beta distribution with range in [0, 1].}

\item{verbose}{If \code{TRUE} reports progress, otherwise it is silent.}
}
\value{
A list containing the following elements:
\itemize{
\item \code{X}: The simulated genotype matrix.
\item \code{bim}: The variant info table.
\item \code{p_anc}: The ancestral allele frequency vector (aligned with rows of \code{X}).
\item \code{kinship}: The true kinship of the joint admixed family model (final generation only).
\item \code{admix_proportions}: The true admixture proportions of the joint admixed family model (final generation only).
\item \code{trait}: The simulated trait vector.
\item \code{causal_indexes}: The vector of randomly selected causal loci indexes.
\item \code{causal_coeffs}: The vector of simulated regression coefficients for causal loci.
}
}
\description{
Simulate entire parameters and random variables, particularly random genotypes and phenotypes, with desired dimensions and other base parameters.
}
\examples{
# run with default values, except smaller
data <- sim_gen_phen( n_ind = 50, G = 3, m_loci = 100, m_causal = 5 )

# main objects of interest:
# genotype matrix
data$X
# trait vector
data$trait

}
\seealso{
This function is a wrapper around \code{\link[=sim_pop]{sim_pop()}}, \code{\link[=sim_geno]{sim_geno()}}, \code{\link[=sim_bim]{sim_bim()}}, and \code{\link[=sim_trait_env]{sim_trait_env()}}, see those for more details.
}
