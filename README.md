# MuXTalk

## What MuXTalk does, in a nutshell:
The basic workflow of MuXTalk consists of the following functions:
- **process_data(...):** Prepare the KEGG signaling pathway, protein-protein interaction (PPI) and gene regulatory network (GRN) data to build the signaling-regulatory multilayer network. This step involves parsing data, matching identifiers, reindexing, creating edge and node dictionaries, and so on, to be suitable for representation as a multilayer network.
- **sparse_layers(...):** Generate sparse adjacency matrices for the signaling and regulatory layer.
- **randomize_\*(...) or randomize_\*\_npz(...):** Create (or read from file) ensembles of randomized versions of the sparse adjacency matrices. Please note that the first version of the randomize function saves the entire ensemble in a .pickle file as a single dictionary for fast access and, therefore, the interaction-specific sparse matrices KEGG_e in particular is a large file (28GB). In contrast, the \_npz version saves each randomized sparse adjacency matrix as an individual .npz (a compressed NumPy array format) file. The overall size of .npz files is much smaller compared to .pickle (e.g., ~5GB for KEGG_e), however, reading .npz files in real-time takes longer than accessing the elements of the dictionary. This is why we use this npz version in the Docker image (detailed below). In either case, the user may choose to download the network ensembles directly using the links below or create them from scratch locally. The latter option takes ~10hrs on a typical laptop.
- **between/shortest_paths_multilink_counts_discovery(...):** Count multilinks for the actual layers.
- **between/shortest_paths_multilink_counts_rand_discovery(...):** Count multilinks for the randomized layers.
- **between/shortest_paths_multilink_zscores_pvals_discovery(...):** Calculate multilink statistics.
- **get_ranked_pathway_pairs_discovery(...):** Calculate MuXTalk scores to prioritize pathway pairs by their propensity to crosstalk with each other based on multilink statistics. The resulting file is a ranked list of all pathway pairs in .csv (default) or .parquet (used for the Streamlit app) format.



## Running MuXTalk locally (recommended for user-defined GRNs)

There are two options to run MuXTalk end-to-end: 
1) Run MuXTalk 

2) Run the MuXTalk script directly (requires 


While the ensemble of randomized network layers can be generated locally

Download the ensembles of randomized networks (N=500) using the below links:
https://www.dropbox.com/s/rbwi4qo2rsx2tqg/A_KEGG_e_sparr_rand_dict_500runs.pickle?dl=0 [KEGG_e -- interaction-specific KEGG]
https://www.dropbox.com/s/pknam5ok2rhccre/A_KEGGPPI_sparr_rand_dict_500runs.pickle?dl=0 [KEGG + PPI combined]
https://www.dropbox.com/s/p1xq59s44rgzxof/HumanGRN10e4_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-4]
https://www.dropbox.com/s/hd49mgdle9uhc33/HumanGRN10e5_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-5]
https://www.dropbox.com/s/z14qywvqsa3mkqw/HumanGRN10e6_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-6]

## Running MuXTalk as a Docker container (Requires >10GB of RAM)

We recommending using Docker as the most user-friendly option of running MuXTalk. 

To run MuXTalk as a Docker container, Docker must be installed. Please see https://docs.docker.com/get-docker/ for installation instructions. Once Docker is installed, the MuXTalk image can either be pulled from Docker Hub or run locally using the .tar file.

Next, follow the below steps to run the MuXTalk container in Docker:
- Download the file (...) to be used as a mounted volume.
- Then simply type in the below command in the terminal:
docker run -it --memory 12g --memory-swap -1 -v /path/to/MuXTalk_Docker_mount/:/MuXTalk_app/ --rm muxtalk-docker-app-slim --proj_path=/MuXTalk_app/ --input_GRN=HumanGRN10e6 --MuXTalk_method=MuXTalk_shortest --get_n=150 --get_randomly=True --sp_threshold=1 --parquet=False

For user-defined GRNs, use "customGRN" as the input GRN name, i.e. --input_GRN=customGRN.




## Running MuXTalk with custom GRNs
- Make GRN edgelist file with the name "input_GRN_edges.csv" where input_GRN is the name of the custom GRN, e.g. "HumanGRN10e6". No column names.




## Running the MuXTalk Streamit app (recommended to visually explore both the default GRNs and custom GRNs)

.parquet files that are output by the MuXTalk script can be used as input to the Streamlit app. 



