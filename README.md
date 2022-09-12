# MuXTalk

The basic workflow of MuXTalk is as follows:
- **process_data(...):** Prepare the KEGG signaling pathway, protein-protein interaction (PPI) and gene regulatory network (GRN) data to build the signaling-regulatory multilayer network. This step involves parsing data, matching identifiers, reindexing, creating edge and node dictionaries, and so on, to be suitable for representation as a multilayer network.
- **sparse_layers(...):** Generate sparse adjacency matrices for the signaling and regulatory layer.
- **randomize_\*(...):** Create (or read from file) ensembles of randomized versions of the sparse adjacency matrices. Please note that the interaction-specific sparse matrices KEGG_e in particular is a large file (28GB). The user may choose to download these using the links below or create them locally. The latter option takes ~10hrs on a typical laptop.
- **between/shortest_paths_multilink_counts_discovery(...):** Count multilinks for the actual layers.
- **between/shortest_paths_multilink_counts_rand_discovery(...):** Count multilinks for the randomized layers.
- **between/shortest_paths_multilink_zscores_pvals_discovery(...):** Calculate multilink statistics.
- **get_ranked_pathway_pairs_discovery(...):** Calculate MuXTalk scores to prioritize pathway pairs by their propensity to crosstalk with each other based on the multilink statistics.



## Running MuXTalk locally (recommended for custom GRNs)

The user has two options to run MuXTalk end-to-end. 


While the ensemble of randomized network layers can be generated locally



## Running the MuXTalk Streamit app (recommended to visually explore the default GRNs)
