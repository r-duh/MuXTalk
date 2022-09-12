# MuXTalk

The basic workflow of MuXTalk is as follows:
- Prepare the KEGG signaling pathway, protein-protein interaction (PPI) and gene regulatory network (GRN) data to build the signaling-regulatory multilayer network. This step involves parsing data, matching identifiers, reindexing, creating edge and node dictionaries, and so on, to be suitable for representation as a multilayer network. [Responsible function: process_data()] 
- Generate sparse adjacency matrices for the signaling and regulatory layer. function: sparse_layers()
- Create ensembles of randomized versions of the sparse adjacency matrices. function: 
- Count multilinks and calculate multilink statistics.

## Running MuXTalk locally (recommended for custom GRNs)

The user has two options to run MuXTalk end-to-end. 


While the ensemble of randomized network layers can be generated locally



## Running the MuXTalk Streamit app (recommended to visually explore the default GRNs)
