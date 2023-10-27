Repository contaning data and codes to reproduce the comparison of sample-wise methods in the rRoma article.

# Pathway database

GMT file: h.all.v2023.1.Hs.symbols.gmt

# Scenario #1: comparison of CF and healthy control samples

- gene expression data: Rehman_CF_NCF_RNA_exp_matrix_norm.rda
- samples information : Rehman_CF_NCF_sample_labels.rda

- rROMA script: rROMA vignette (https://sysbio-curie.github.io/rROMA/index.html)
- comparison with the other state of the art methods: Rehman_CF_NCF_other_methods.R

# Scenario #2: comparison of healthy control samples before and after treatment with IL17+TNFalpha

- gene expression data: Rehman_IL17_TNFa_RNA_exp_matrix_norm.rda
- samples information : Rehman_IL17_TNFa_sample_labels.rda

- exploratory analysis and normalisation of data matrix: Rehman_IL17_TNFa_RNA_exp_matrix_norm.R
- rROMA script: Rehman_TNF_IL17_rROMA.R
- comparison with the other state of the art methods: Rehman_TNF_IL17_other_methods.R

