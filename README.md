<!-- To-do:
- Hover functionality for the pyvis visualization -- https://stackoverflow.com/questions/70546467/further-editing-for-items-in-the-pyvis-tooltip
- Option to limit the # of sps for sp>1 for visual clarity
-  print out terminal in the progress status bar for details
-   -->


# MuXTalk: Detecting and dissecting signaling crosstalk via multilayer networks. 

A network-based statistical framework to infer and explore crosstalk between signaling pathways.

**Preprint URL:** https://www.biorxiv.org/content/10.1101/2022.09.29.510183v1

![Overview of the MuXTalk framework](https://github.com/r-duh/MuXTalk/blob/main/MuXTalk_Overview.png?raw=true)

## Running MuXTalk locally (recommended for user-defined GRNs)

### There are two options to run MuXTalk end-to-end: 
#### 1) Run MuXTalk as a Docker container (requires Docker installation)
This is perhaps the most straightforward choice for the general user and only requires Docker to be installed. Please see https://docs.docker.com/get-docker/ for installation instructions. Note: The Docker option is somewhat slower to run but has the advantage of not depending on the specific package environment and operating system of the user.

- Once Docker is installed, the MuXTalk image "rduh/muxtalk:slim" can either be pulled from Docker Hub
```
docker pull rduh/muxtalk:slim
```
or loaded from the [.tar file](https://www.dropbox.com/scl/fi/15m7nuob3d325h4sgelcb/muxtalk_slim.tar?dl=0&rlkey=jl3wi9loats4sz4ezeyj114ok) using
```
docker load --input muxtalk_slim.tar
```

- [Download](https://www.dropbox.com/sh/jhir7f3srikhjxj/AAC_Mjk2tzJ6OiMdJp7PEg9Ua?dl=0) the MuXTalk folder to be mounted as a volume to the Docker container. This local folder (i.e., located in the user's machine), named /MuXTalk_Docker_mount/, will act as the main folder in which the MuXTalk container will read and write files.

- As an example use case, type in the below command in the terminal to run the MuXTalk image as a container. /path/to/MuXTalk_Docker_mount/ is where the folder you downloaded is located in your computer. Details about the MuXTalk parameters can be found in the following sections.
```
docker run -it -v /path/to/MuXTalk_Docker_mount/:/MuXTalk_app/ --rm rduh/muxtalk:slim --proj_path=/MuXTalk_app/ --input_GRN=HumanGRN10e6 --input_PPI=PPI_Cheng_NatComms2019 --input_PPI_ID_format=Entrez --MuXTalk_method=MuXTalk_shortest --sp_threshold=1
```

#### 2) Run the MuXTalk script directly (requires conda to be installed). 
This is the faster option but requires familiarity with creating environments and running scripts. As the first step, [install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) on your system . Next, type the following commands in your command prompt, or terminal, in the following order to set up and run MuXTalk.
- Using Conda, create a new environment named "muxtalk" (or any name of your choosing) for MuXTalk and install in it the dependencies needed by MuXTalk using the MuXTalk_requirements.yml file in /path/to/MuXTalk_Docker_mount/ (also available in the project GitHub [page](https://github.com/r-duh/MuXTalk/blob/main/MuXTalk_environment.yml)): 
```
conda env create -n muxtalk -f /path/to/MuXTalk_Docker_mount/MuXTalk_environment.yml
```
- Activate the newly created Conda environment:
```
conda activate muxtalk
```
- You have now created a new conda environment and installed in it all the packages MuXTalk needs to run. The only remaining step is to run MuXTalk. First, navigate to /path/to/MuXTalk_Docker_mount/
```
cd /path/to/MuXTalk_Docker_mount/
 ```
You can then inspect the options for the MuXTalk script by 
```
python run_MuXTalk_v2.py --help
```
which will return
```
usage: run_MuXTalk_v2.py [-h] [--proj_path PROJ_PATH] [--input_GRN INPUT_GRN]
                         [--input_GRN_ID_format INPUT_GRN_ID_FORMAT]
                         [--input_PPI INPUT_PPI]
                         [--input_PPI_ID_format INPUT_PPI_ID_FORMAT]
                         [--MuXTalk_method MUXTALK_METHOD] [--get_n [GET_N]]
                         [--get_randomly [GET_RANDOMLY]]
                         [--sp_threshold [SP_THRESHOLD]] [--parquet [PARQUET]]

optional arguments:
  -h, --help            show this help message and exit
  --proj_path PROJ_PATH
                        Path to the project file where all the input/output
                        files will be located.
  --input_GRN INPUT_GRN
                        Name of the input gene regulatory network (GRN).
                        Please note that this must be the same as the name
                        preceding '_edges.csv' in the GRN edgelist.
  --input_GRN_ID_format INPUT_GRN_ID_FORMAT
                        ID format for the input GRN. Specify 'Gene_Symbol' or
                        'Entrez' depending on the input data. Default is
                        'Gene_Symbol'.
  --input_PPI INPUT_PPI
                        Name of the input protein-protein interaction (PPI)
                        network. Please note that this must be the same as the
                        name preceding '_edges.csv' in the PPI edgelist.
  --input_PPI_ID_format INPUT_PPI_ID_FORMAT
                        ID format for the input PPI. Specify 'Gene_Symbol' or
                        'Entrez' depending on the input data. Default is
                        'Gene_Symbol'.
  --MuXTalk_method MUXTALK_METHOD
                        MuXTalk method to be used: The two valid options are
                        'MuXTalk_between' and 'MuXTalk_shortest'.
  --get_n [GET_N]       Number of randomized instances to be retrieved from
                        the ensemble. For memory efficiency, the default value
                        is set to 100.
  --get_randomly [GET_RANDOMLY]
                        Option to retrieve randomized networks in order or
                        randomly. Default value is True.
  --sp_threshold [SP_THRESHOLD]
                        The shortest path threshold to be used in MuXTalk. The
                        valid options are 'None', '1', or '2'. Default is '1'.
  --parquet [PARQUET]   Option to output the MuXTalk results as a .parquet
                        file. False (default) outputs .csv files.
```

The only required argument is --proj_path, which should be set to /path/to/MuXTalk_Docker_mount/. Set --parquet=True to generate .parquet files to be used with the Streamlit visualization (detailed below).


For example, if we wanted to run MuXTalk_shortest on "HumanGRN106" as the input gene regulatory network (GRN) with a shortest path threshold of 1, we would run MuXTalk using
```
python run_MuXTalk_v2.py --proj_path=/path/to/MuXTalk_Docker_mount/ --input_GRN=HumanGRN10e6 --input_PPI=PPI_Cheng_NatComms2019 --input_PPI_ID_format=Entrez --MuXTalk_method=MuXTalk_shortest --sp_threshold=1
```

Running the above script will result in the main output of MuXTalk, i.e., the two files named "input_GRN_input_PPI_sp_threshold_MuXTalk_method_detected_discovery.csv" and "input_GRN_input_PPI_sp_threshold_MuXTalk_method_discovery.csv" where input_GRN, input_PPI, sp_threshold and MuXTalk_method are the user-defined parameters. These files contain the MuXTalk-ranked lists of pathway pairs for the "detected" (meaning pathway pairs that were assessed by MuXTalk, excluding the pathway pairs that had no shortest paths between them in the PPI network) and complete set, respectively.

## Running MuXTalk with custom GRNs and PPI networks
- Once MuXTalk is set up, we can also run it with user-defined GRNs and PPI networks. For this, we need to first add into our local MuXTalk folder (/path/to/MuXTalk_Docker_mount/) an edgelist file for the custom GRN and/or PPI network named "customGRN_edges.csv" and/or "customPPI_edges.csv" where "customGRN" and "customPPI" can be user-defined names. This file must have two columns, without headers: For GRN: the first one for the source gene (or transcription factor) and the second one for the target gene, as shown below, or for PPI, the first one for Protein A and the second one for Protein B. "customGRN" and "customPPI" will also be the name of the input_GRN and input_PPI variable, respectively. Genes can have Gene Symbols or Entrez IDs as identifiers; MuXTalk will take care of all the ID conversions.

|  |  |
| --- | --- |
| Source Gene A | Target Gene B |
| Source Gene C | Target Gene D |
| Source Gene E | Target Gene F |
| ... | ... |
 
- When we run MuXTalk with --input_GRN=customGRN with customGRN_edges.csv (and/or --input_PPI=customPPI with customPPI_edges.csv) in /path/to/MuXTalk_Docker_mount/, MuXTalk will create the randomized versions of the custom GRN and signaling networks and store them in the /customGRN_customPPI_A_GRN_sparr_rand_npz_files/,  /customPPI_A_KEGG_e_sparr_rand_npz_files/ and /customPPI_A_KEGGPPI_sparr_rand_npz_files/ folders. These folders are created and populated as MuXTalk runs. This step will have to be only done once per each new GRN/PPI combination. Here is an example that uses the same parameters but on a custom GRN and PPI provided by the user. **Note:** When running MuXTalk with custom PPI networks, the user might want to ensure there is enough disk space as the randomized ensembles of the signaling layer (i.e., the A_KEGG_e and A_KEGGPPI files) may take up about 20-30GB in total.
```
docker run -it -v /path/to/MuXTalk_Docker_mount/:/MuXTalk_app/ --rm rduh/muxtalk:slim --proj_path=/MuXTalk_app/ --input_GRN=customGRN --input_PPI=customPPI --input_PPI_ID_format=Gene_Symbol --MuXTalk_method=MuXTalk_shortest --get_n=150 --get_randomly=True --sp_threshold=1 --parquet=True
```


## What MuXTalk does, in a nutshell:
The basic workflow of MuXTalk consists of the following functions:
- **process_data(...):** Prepare the KEGG signaling pathway, protein-protein interaction (PPI) and gene regulatory network (GRN) data to build the signaling-regulatory multilayer network. This step involves parsing data, matching identifiers, reindexing, creating edge and node dictionaries, and so on, to be suitable for representation as a multilayer network.
- **sparse_layers(...):** Generate sparse adjacency matrices for the signaling and regulatory layer.
- **randomize_\*(...) or randomize_\*\_npz(...):** Create (or read from file) ensembles of randomized versions of the sparse adjacency matrices. Please note that the first version of the randomize function saves the entire ensemble in a .pickle file as a single dictionary for fast access and, therefore, the interaction-specific sparse matrices KEGG_e in particular is a large file (~30GB). In contrast, the \_npz version saves each randomized sparse adjacency matrix as an individual .npz (a compressed NumPy array format) file. The overall size of .npz files is much smaller compared to .pickle (e.g., ~5GB for KEGG_e), however, reading .npz files in real-time takes much longer than accessing the elements of the dictionary. This is why we use the npz version in the Docker image (detailed below). In either case, the user may choose to download the network ensembles directly using the links below or create them from scratch locally. The latter option takes ~10hrs on a typical laptop. If the randomized networks are available in the **proj_path** directory, the randomize function will skip to just reading them.
- **between/shortest_paths_multilink_counts_discovery(...):** Count multilinks for the actual layers.
- **between/shortest_paths_multilink_counts_rand_discovery(...):** Count multilinks for the randomized layers.
- **between/shortest_paths_multilink_zscores_pvals_discovery(...):** Calculate multilink statistics.
- **get_ranked_pathway_pairs_discovery(...):** Calculate MuXTalk scores to prioritize pathway pairs by their propensity to crosstalk with each other based on multilink statistics. The resulting file is a ranked list of all pathway pairs in .csv (default) or .parquet (used for the Streamlit app) format.


## Running the MuXTalk Streamit app (recommended to visually explore both the default GRNs and custom GRNs)

The steps to run the MuXTalk Streamlit app are very similar to those of the MuXTalk algorithm above. We have two options: 

#### 1) Run MuXTalk Streamlit app as a Docker container (requires [Docker installation](https://docs.docker.com/get-docker/))
This option is somewhat slower to run but has the advantage of not depending on the specific package environment and operating system of the user.

- Once Docker is installed, the MuXTalk image "rduh/muxtalk-streamlit:slim" can either be pulled from Docker Hub
```
docker pull rduh/muxtalk-streamlit:slim
```
or loaded from the [.tar file](https://www.dropbox.com/s/r1yd6rcl9pyeyrc/muxtalk_streamlit_slim.tar?dl=0) using
```
docker load --input muxtalk_streamlit_slim.tar
```

- [Download](https://www.dropbox.com/sh/671hr100hqpymro/AADNet8iHQCdckMyMSnPMlzma?dl=0) the MuXTalk folder to be mounted as a volume to the Docker container. This local folder (i.e., located in the user's machine), named /MuXTalk_Streamlit_Docker_mount/, will act as the main folder in which MuXTalk Streamlit app's container will read and write files.

- To run the MuXTalk Streamlit app as a Docker container, type in the below command in the terminal. /path/to/MuXTalk_Streamlit_Docker_mount/ is where the folder you downloaded above is located in your computer. Details about the MuXTalk visualization parameters can be found in the readme within the app.
```
docker run -it -v /path/to/MuXTalk_Streamlit_Docker_mount/:/MuXTalk_Streamlit_app/ -p 8501:8501 rduh/muxtalk-streamlit:slim
```

- Finally, to see the Streamlit app, go to your browser and enter the address that appears in your terminal (It looks like this: "You can now view your Streamlit app in your browser. URL: `http://0.0.0.0:8501` "). So, typically `http://0.0.0.0:8501`. If you have more than one Streamlit instance running, this can be `http://0.0.0.0:8502`, `http://0.0.0.0:8503`, and so on.


#### 2) Run the MuXTalk Streamlit app script directly (requires conda to be installed). 
This is the faster option but requires familiarity with creating environments and running scripts.
- If not already done so, install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) on your system. Next, type the following commands in your command prompt, or terminal, in the following order to set up and run MuXTalk.
- Using Conda, create a new environment named "muxtalk_streamlit" (or any name of your choosing) for the MuXTalk Streamlit app and install in it the dependencies needed by the app using the MuXTalk_Streamlit_requirements.yml file in /path/to/MuXTalk_Streamlit_Docker_mount/ (also available in the project GitHub [page](https://github.com/r-duh/MuXTalk/blob/main/Streamlit/MuXTalk_Streamlit_environment.yml)): 
```
conda env create -n muxtalk_streamlit -f /path/to/MuXTalk_Streamlit_Docker_mount/MuXTalk_Streamlit_environment.yml
```
- Activate the newly created Conda environment:
```
conda activate muxtalk_streamlit
```
- You have now created a new conda environment and installed in it all the packages the MuXTalk Streamlit app needs to run. The only remaining step is to run it. First, navigate to /path/to/MuXTalk_Streamlit_Docker_mount/
```
cd /path/to/MuXTalk_Streamlit_Docker_mount/
```
Then, run the streamlit run by typing the below command (note that it has the additional --proj_path flag, which needs to be set to the MuXTalk Streamlit app directory
```
streamlit run MuXTalk_viz_Streamlit_Docker.py -- --proj_path=/path/to/MuXTalk_Streamlit_Docker_mount/
```

- Finally, to see the Streamlit app, go to your browser and enter the address that appears in your terminal (It looks like this: "You can now view your Streamlit app in your browser. URL: `http://0.0.0.0:8501` "). So, typically `http://0.0.0.0:8501`. If you have more than one Streamlit instance running, this can be `http://0.0.0.0:8502`, `http://0.0.0.0:8503`, and so on.


### Visualizing custom GRNs using the MuXTalk Streamlit app
The .parquet files that are output by the MuXTalk script (details described above) can be used as input to the Streamlit app by a simple drag-and-drop on the app.


## Troubleshooting

In the event that the Docker container quits with or without errors, the following measures might help to mitigate the potentially memory-related issues.
- Use the --memory and --memory-swap flags to increase the memory and swap partition allocated to the container
-   e.g. - Then simply type in the below command in the terminal: 
```
docker run -it --memory 12g --memory-swap -1 -v /path/to/MuXTalk_Docker_mount/:/MuXTalk_app/ --rm rduh/muxtalk:slim --proj_path=/MuXTalk_app/ --input_GRN=HumanGRN10e6 --MuXTalk_method=MuXTalk_shortest --get_n=150 --get_randomly=True --sp_threshold=1 --parquet=False
```

## Manuscript analyses
The analyses presented in our paper can be viewed via the link below:

https://nbviewer.org/github/r-duh/MuXTalk/blob/main/manuscript_analysis/R2_MuXTalk_manuscript_analyses.ipynb

<br>
<br>
<br>

<!-- Additional considerations:

While the ensemble of randomized network layers can be generated locally

Download the ensembles of randomized networks (N=500) using the below links:
https://www.dropbox.com/s/rbwi4qo2rsx2tqg/A_KEGG_e_sparr_rand_dict_500runs.pickle?dl=0 [KEGG_e -- interaction-specific KEGG]
https://www.dropbox.com/s/pknam5ok2rhccre/A_KEGGPPI_sparr_rand_dict_500runs.pickle?dl=0 [KEGG + PPI combined]
https://www.dropbox.com/s/p1xq59s44rgzxof/HumanGRN10e4_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-4]
https://www.dropbox.com/s/hd49mgdle9uhc33/HumanGRN10e5_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-5]
https://www.dropbox.com/s/z14qywvqsa3mkqw/HumanGRN10e6_A_GRN_sparr_rand_dict_500runs.pickle?dl=0 [GRN, p<10e-6] -->

# Contact:
Created and maintained by Arda Halu. For requests for assistance, email arda.halu@channing.harvard.edu.



