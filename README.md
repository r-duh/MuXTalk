# MuXTalk

The basic workflow of MuXTalk consists of the following functions:
- **process_data(...):** Prepare the KEGG signaling pathway, protein-protein interaction (PPI) and gene regulatory network (GRN) data to build the signaling-regulatory multilayer network. This step involves parsing data, matching identifiers, reindexing, creating edge and node dictionaries, and so on, to be suitable for representation as a multilayer network.
- **sparse_layers(...):** Generate sparse adjacency matrices for the signaling and regulatory layer.
- **randomize_\*(...):** Create (or read from file) ensembles of randomized versions of the sparse adjacency matrices. Please note that the interaction-specific sparse matrices KEGG_e in particular is a large file (28GB). The user may choose to download these using the links below or create them locally. The latter option takes ~10hrs on a typical laptop.
- **between/shortest_paths_multilink_counts_discovery(...):** Count multilinks for the actual layers.
- **between/shortest_paths_multilink_counts_rand_discovery(...):** Count multilinks for the randomized layers.
- **between/shortest_paths_multilink_zscores_pvals_discovery(...):** Calculate multilink statistics.
- **get_ranked_pathway_pairs_discovery(...):** Calculate MuXTalk scores to prioritize pathway pairs by their propensity to crosstalk with each other based on the multilink statistics. The resulting file is a ranked list of all pathway pairs in .csv (default) or .parquet format.



## Running MuXTalk locally (recommended for user-defined GRNs)

The user has two options to run MuXTalk end-to-end. 


While the ensemble of randomized network layers can be generated locally

Download the ensembles of randomized networks from 
https://www.dropbox.com/s/rbwi4qo2rsx2tqg/A_KEGG_e_sparr_rand_dict_500runs.pickle?dl=0 [KEGG_e -- interaction-specific KEGG]
https://www.dropbox.com/s/pknam5ok2rhccre/A_KEGGPPI_sparr_rand_dict_500runs.pickle?dl=0 [KEGG + PPI combined]
https://www.dropbox.com/s/p1xq59s44rgzxof/HumanGRN10e4_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-4]
https://www.dropbox.com/s/hd49mgdle9uhc33/HumanGRN10e5_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-5]
https://www.dropbox.com/s/z14qywvqsa3mkqw/HumanGRN10e6_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-6]


## Running the MuXTalk Streamit app (recommended to visually explore the default GRNs)
