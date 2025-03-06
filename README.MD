# A lack of repeated evolution of *ben-1* mediated benzimidazole resistance in *Caenorhabditis* nematodes

Data, analysis scripts, and results from the article "A lack of repeated evolution of *ben-1* mediated benzimidazole resistance in *Caenorhabditis* nematodes"

## Abstract
Benzimidazoles (BZs), a widely used class of anthelmintic drugs, target beta-tubulin proteins, disrupt microtubule formation, and cause nematode death. In parasitic nematode species, mutations in beta-tubulin genes (e.g., isotype-1 beta-tubulin) are predicted to inhibit BZ binding and are associated with BZ resistance. Similarly, in the free-living nematode Caenorhabditis elegans, mutations in an isotype-1 beta-tubulin ortholog, ben-1, are the primary drivers of BZ resistance. The recurrent association of BZ resistance with alleles predicted to impact beta-tubulin function across nematode species suggests that BZ resistance repeatedly evolves through variants in conserved beta-tubulin genes and provides compelling evidence to predict the emergence of BZ resistance across nematode species. To evaluate the repeated evolution of BZ resistance mediated by conserved beta-tubulin genes, we identified predicted BZ resistance alleles in beta-tubulin genes across wild strains from three Caenorhabditis species: C. elegans, Caenorhabditis briggsae, and Caenorhabditis tropicalis. We hypothesized that if these species experienced similar selective pressures that they would evolve resistance to BZs through mutations in three beta-tubulin genes (ben-1, tbb-1, and tbb-2), which are expressed in multiple tissues and implicated in tubulin polymerization dynamics. Using high-throughput imaging assays, we tested the association of predicted BZ resistance alleles with BZ resistance. We found that a heterogenous set of variants identified in Cel-ben-1 were highly associated with BZ resistance. In C. briggsae, only two detected variants in Cbr-ben-1 (W21stop and Q134H) were associated with BZ resistance. In C. tropicalis, two missense variants were identified in ben-1, but neither was associated with BZ resistance. To determine whether C. briggsae and C. tropicalis evolved BZ resistance through other beta-tubulin genes, we tested whether predicted BZ resistance alleles in tbb-1 and tbb-2 were associated with BZ resistance. Variants in tbb-1 or tbb-2 in C. briggsae and C. tropicalis were not associated with BZ resistance. Our findings reveal a lack of repeated evolution of BZ resistance across the three Caenorhabditis species and highlight the importance of defining BZ resistance mechanisms outside the beta-tubulin gene family. 


## Code
### Required software: 


# Scripts
`expression_variation` - folder containing scripts to generate figures for 
    - S1. The realationship between ben-1 expression levels and ABZ responses 
    in C. elegans strains.
    - S2. The relationship between tbb-1 (S2A) and tbb-2 (S2B) expression levels and ABZ responses
`phylogeographic_distribution` - folder containing scripts to generate figures for 5, S14, and S15
`hta_res_threshold` - folder containing scripts to set the resistance threshold for HTA data collected in paper. 
`pull_bg.R` - Script to examine Blossum and Grantham scores for beta-tubulin missense variants. Generates figures S10, S11, S12, and S13. 

# Figures

Figure S1:  The relationship between ben-1 expression levels and albendazole response in C. elegans strains.
  - `scripts/expression_variation/beta_tub_expression_plots.R`

Figure S2: The relationship between C. elegans beta-tubulin expression and albendazole response
  - `scripts/expression_variation/beta_tub_expression_plots.R`


# Data

- `blast_summary`: Summary of the recipoval best hits for beta-tubulin genes in C. elegans, C. briggsae, and C. tropicalis.
- `caendr_annotation`: Contains annotation files downloaded from CaeNDR
- `hta_summaries`: Summary of ABZ responses for strains assyed in this paper.