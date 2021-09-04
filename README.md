# simgenphen

`simgenphen` (Simulate Genotypes and Phenotypes) is a wrapper package around several more specialized packages: `bnpsd` for admixture, `simfam` for family structure, and `simtrait` for complex traits.
This package assumes reasonable defaults and generally saves on typing for some recurrent cases.
Key non-trivial functionality is drawing genotypes for admixed families whose loci are all polymorphic.

## Installation

The current development version can be installed from the GitHub repository using `devtools`:
```R
install.packages("devtools") # if needed
library(devtools)
install_github('OchoaLab/simgenphen')
```


## Examples

There are essentially two recommended ways of using this code, shown in the next two subsections.
``` r
library(simgenphen)
```

### One-liner version

Here is a simple way of generating genotypes that feature both admixture and family structure, along with a phenotype!
It's run mostly with default values, except data dimensions are smaller:
```r
data <- sim_gen_phen( n_ind = 50, G = 3, m_loci = 100, m_causal = 5 )

# main objects of interest:
# genotype matrix
data$X
# trait vector
data$trait
```

### Version with separate steps

This approach requires more commands, but is more logical and efficient in loops where population structure is constant but want newly-drawn genotype and phenotype replicates.
These more explicit commands also return additional values.

First construct population parameters (relatedness but no loci):
```r
k_subpops <- 3
data_pop <- sim_pop( n_ind = 50, G = 3, k_subpops = k_subpops )
```
Now draw genotypes!
```r
m_loci <- 100
data_geno <- sim_geno(
    data_pop$admix_proportions_1,
    data_pop$inbr_subpops,
    data_pop$fam,
    data_pop$ids,
    m_loci
)

# genotype matrix
data_geno$X
```
Lastly, draw trait!
In this case we specify a non-default heritability and an environment effect "gcat" to showcase that version (those parameters can also be passed to the "one-liner" function `sim_gen_phen`).
```r
m_causal <- 5
herit <- 0.5
data_trait <- sim_trait_env(
    data_geno$X,
    data_geno$p_anc,
    m_causal,
    herit,
    env = 'gcat',
    k_subpops = k_subpops
)

# the simulated trait
data_trait$trait
# indexes of randomly-selected causal loci
data_trait$causal_indexes
# coefficients of causal loci
data_trait$causal_coeffs
```
