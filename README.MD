# Independent mechanisms of benzimidazole resistance across *Caenorhabditis* nematodes

Data, analysis scripts, and results from the article "Independent mechanisms of benzimidazole resistance across *Caenorhabditis* nematodes"

## Abstract
Benzimidazoles (BZs), a widely used class of anthelmintic drugs, target beta-tubulin proteins, disrupt microtubule formation, and cause nematode death. In parasitic nematode species, mutations in beta-tubulin genes (*e.g.*, isotype-1 beta-tubulin) are predicted to inhibit BZ binding and are associated with BZ resistance. Similarly, in the free-living nematode *Caenorhabditis elegans*, mutations in an isotype-1 beta-tubulin ortholog, *ben-1*, are the primary drivers of BZ resistance. The recurrent association of BZ resistance with beta-tubulins suggests that BZ resistance is repeatedly caused by mutations in beta-tubulin genes, an example of repeated evolution of drug resistance across nematode species. To evaluate the hypothesis of repeated evolution of BZ resistance mediated by beta-tubulin, we identified predicted resistance alleles in beta-tubulin genes across wild strains from three *Caenorhabditis* species: *C. elegans*, *Caenorhabditis briggsae*, and *Caenorhabditis tropicalis*. We hypothesized that, if these species experienced similar selective pressures, they would evolve resistance to BZs by mutations in any of three beta-tubulin genes (*ben-1*, *tbb-1*, and *tbb-2*). Using high-throughput development assays, we tested the association of predicted beta-tubulin alleles with BZ resistance. We found that a heterogeneous set of variants identified in *C. elegans ben-1* were associated with BZ resistance. In *C. briggsae*, only two variants in *ben-1*, predicted to encode a premature stop codon (W21stop) and a missense substitution (Q134H), were associated with BZ resistance. In *C. tropicalis*, two missense variants were identified in *ben-1*, but neither was associated with BZ resistance. *C. briggsae* and *C. tropicalis* might have evolved BZ resistance by mutations in other beta-tubulin genes, but we found that variants in *tbb-1* or *tbb-2* in these species were not associated with BZ resistance. Our findings reveal a lack of repeated evolution of BZ resistance across the three *Caenorhabditis* species and highlight the importance of defining BZ resistance mechanisms outside of beta-tubulins.


<!-- ## Code
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
- `hta_summaries`: Summary of ABZ responses for strains assyed in this paper. -->
