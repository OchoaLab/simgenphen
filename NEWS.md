# simgenphen 0.0.0.9000 (2021-08-27)

- First commit of public version
- Includes full tests for all functions
- Nothing exported or documented yet

# simgenphen 0.0.1.9000 (2021-08-27)

- Documented and exported all functions!
- Function `sim_pop` now also returns `admix_proportions_1`, the admixture proportions of the founders, distinguished from `admix_proportions` which is the admixture proportions of the final generation (different if `G > 1`).
- Function `sim_geno`:
  - First parameter is renamed from `admix_proportions` to `admix_proportions_1` to emphasize that it should be the admixture proportions of the founders!
  - Optional parameter `p_anc_input` is renamed to `p_anc`.
- Function `sim_gen_phen` debugged to use founder admixture proportions to draw genotypes for founders.
  The previous version had the bug because `sim_geno` did not return founder admixture proportions, only last-generation values, but the bug only manifested if `G > 1`.
  Bug was very recently introduced (compared to previous non-public code in a separate repository).
- Function `sim_bim` added parameter `pos_gap`.

# simgenphen 0.0.2.9000 (2021-09-03)

- Added package documentation with runnable examples.
- Also updated examples in individual functions to better match the package entry examples.

# simgenphen 0.0.3.9000 (2021-09-03)

- Added README with examples from package documentation.

# simgenphen 0.0.4.9000 (2022-02-21)

- Function `sim_geno` fixed a bug that fixed loci were redrawn without respecting `p_anc` if provided (no bug for `p_anc = NULL`).
- Functions `sim_geno` and `sim_gen_phen` added option `beta` for a Beta distribution for ancestral allele frequencies.

# simgenphen 0.0.5.9000 (2023-02-01)

- Function `sim_gen_phen` added element `p_anc` to return list, needed in order to reuse the simulated genotypes to simulate new traits for them (for example, with different heritabilities).
- Updated unit tests for function `sim_bim` to expect column `alt` before `ref` in return tibble (a change inherited from dependency `genio::make_bim()` made in `genio` version 1.1.0.9000 (2022-04-20))
