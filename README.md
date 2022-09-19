# MuXTalk

## What MuXTalk does, in a nutshell:
The basic workflow of MuXTalk consists of the following functions:
- **process_data(...):** Prepare the KEGG signaling pathway, protein-protein interaction (PPI) and gene regulatory network (GRN) data to build the signaling-regulatory multilayer network. This step involves parsing data, matching identifiers, reindexing, creating edge and node dictionaries, and so on, to be suitable for representation as a multilayer network.
- **sparse_layers(...):** Generate sparse adjacency matrices for the signaling and regulatory layer.
- **randomize_\*(...) or randomize_\*\_npz(...):** Create (or read from file) ensembles of randomized versions of the sparse adjacency matrices. Please note that the first version of the randomize function saves the entire ensemble in a .pickle file as a single dictionary for fast access and, therefore, the interaction-specific sparse matrices KEGG_e in particular is a large file (28GB). In contrast, the \_npz version saves each randomized sparse adjacency matrix as an individual .npz (a compressed NumPy array format) file. The overall size of .npz files is much smaller compared to .pickle (e.g., ~5GB for KEGG_e), however, reading .npz files in real-time takes much longer than accessing the elements of the dictionary. This is why we use the npz version in the Docker image (detailed below). In either case, the user may choose to download the network ensembles directly using the links below or create them from scratch locally. The latter option takes ~10hrs on a typical laptop. If the randomized networks are available in the **proj_path** directory, the randomize function will skip to just reading them.
- **between/shortest_paths_multilink_counts_discovery(...):** Count multilinks for the actual layers.
- **between/shortest_paths_multilink_counts_rand_discovery(...):** Count multilinks for the randomized layers.
- **between/shortest_paths_multilink_zscores_pvals_discovery(...):** Calculate multilink statistics.
- **get_ranked_pathway_pairs_discovery(...):** Calculate MuXTalk scores to prioritize pathway pairs by their propensity to crosstalk with each other based on multilink statistics. The resulting file is a ranked list of all pathway pairs in .csv (default) or .parquet (used for the Streamlit app) format.



## Running MuXTalk locally (recommended for user-defined GRNs)

There are two options to run MuXTalk end-to-end: 
1) Run MuXTalk as a Docker container (requires Docker installation)
This is perhaps the most straightforward choice for the general user and only requires Docker to be installed. Please see https://docs.docker.com/get-docker/ for installation instructions. 

- Once Docker is installed, the MuXTalk image "muxtalk-docker-app-slim" can either be pulled from Docker Hub or run locally using the .tar file in the link below:
[docker pull -----]
[docker load -----]

- [Download](https://www.dropbox.com/sh/ztlc8spxyvu5cgn/AABVSaaTLQQUrs3_SwLo-B8ca?dl=0) the MuXTalk folders to be mounted as a volume to the Docker container. This local folder (i.e., located in the user's machine) will act as the main folder in which the MuXTalk container will read and write files.

- Simply type in the below command in the terminal to run the MuXTalk image as a container:
```
docker run -it --memory 12g --memory-swap -1 -v /path/to/MuXTalk_Docker_mount/:/MuXTalk_app/ --rm muxtalk-docker-app-slim --proj_path=/MuXTalk_app/ --input_GRN=HumanGRN10e6 --MuXTalk_method=MuXTalk_shortest --get_n=150 --get_randomly=True --sp_threshold=1 --parquet=False
```

2) Run the MuXTalk script directly (requires familiarity with Python environments)
- Create a new environment for MuXTalk using conda: 
- conda create -n myenv python=3.9


While the ensemble of randomized network layers can be generated locally

Download the ensembles of randomized networks (N=500) using the below links:
https://www.dropbox.com/s/rbwi4qo2rsx2tqg/A_KEGG_e_sparr_rand_dict_500runs.pickle?dl=0 [KEGG_e -- interaction-specific KEGG]
https://www.dropbox.com/s/pknam5ok2rhccre/A_KEGGPPI_sparr_rand_dict_500runs.pickle?dl=0 [KEGG + PPI combined]
https://www.dropbox.com/s/p1xq59s44rgzxof/HumanGRN10e4_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-4]
https://www.dropbox.com/s/hd49mgdle9uhc33/HumanGRN10e5_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-5]
https://www.dropbox.com/s/z14qywvqsa3mkqw/HumanGRN10e6_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-6]


For user-defined GRNs, use "customGRN" as the input GRN name, i.e. --input_GRN=customGRN.



The options for the MuXTalk script can be accessed by 
```
python3 run_MuXTalk.py --help
```
which will return
```
usage: run_MuXTalk.py [-h] [--proj_path PROJ_PATH] [--input_GRN INPUT_GRN]
                      [--MuXTalk_method MUXTALK_METHOD] [--get_n GET_N]
                      [--get_randomly GET_RANDOMLY]
                      [--sp_threshold SP_THRESHOLD] [--parquet PARQUET]

optional arguments:
  -h, --help            show this help message and exit
  --proj_path PROJ_PATH
                        Path to the project file where all the input/output
                        files will be located.
  --input_GRN INPUT_GRN
                        Name of the input gene regulatory network (GRN).
                        Please note that this must be the same as the name
                        preceding '_edges.csv' in the GRN edgelist.
  --MuXTalk_method MUXTALK_METHOD
                        MuXTalk method to be used: The two valid options are
                        'MuXTalk_between' and 'MuXTalk_shortest'.
  --get_n GET_N         Number of randomized instances to be retrieved from
                        the ensemble. For memory efficiency, get_n=100 is
                        recommended.
  --get_randomly GET_RANDOMLY
                        Option to retrieve randomized networks in order or
                        randomly. Default value is True.
  --sp_threshold SP_THRESHOLD
                        The shortest path threshold to be used in MuXTalk. The
                        valid options are 'None', '1', or '2'.
  --parquet PARQUET     Option to output the MuXTalk results as a .parquet
                        file. False (default) outputs .csv files.
```


## Running MuXTalk with custom GRNs
- Make GRN edgelist file with the name "input_GRN_edges.csv" where input_GRN is the name of the custom GRN, e.g. "HumanGRN10e6". No column names.




## Running the MuXTalk Streamit app (recommended to visually explore both the default GRNs and custom GRNs)

.parquet files that are output by the MuXTalk script can be used as input to the Streamlit app. 



## Troubleshooting

In the event that the Docker container quits with or without errors, the following to mitigate the potentially memory-related issues.
- Use the --memory and --memory-swap flags to increase the memory and swap partition allocated to the container
-   e.g. - Then simply type in the below command in the terminal: 
```
docker run -it --memory 12g --memory-swap -1 -v /path/to/MuXTalk_Docker_mount/:/MuXTalk_app/ --rm muxtalk-docker-app-slim --proj_path=/MuXTalk_app/ --input_GRN=HumanGRN10e6 --MuXTalk_method=MuXTalk_shortest --get_n=150 --get_randomly=True --sp_threshold=1 --parquet=False
```
