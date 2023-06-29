import streamlit as st
import pandas as pd
import numpy as np
from itertools import product
from tqdm import tqdm
from scipy import sparse
import json

# Input must be in the form of GRN edgelist (tab-separated) with two columns for source (TF) and target, respectively, and with no headers. 
@st.cache(show_spinner=False, allow_output_mutation=True)
def process_GRN_data(proj_path, input_GRN, input_filenames_dict, ID_format='Gene_Symbol'):

	print('Processing GRN edgelist...')
	
	if ID_format=='Entrez':
		GRN_edges_df_entrez = pd.read_csv(proj_path + input_GRN + '_edges.csv', sep='\t', header=None)
		GRN_edges_df_entrez = GRN_edges_df_entrez.rename(columns={0: 'Source_Entrez_ID', 1: 'Target_Entrez_ID'})
			
		return GRN_edges_df_entrez

	elif ID_format=='Gene_Symbol':
	
		HUGO_symb_entrez_uniprot = pd.read_csv(proj_path + input_filenames_dict['HUGO_symb_entrez_uniprot'], sep='\t')
		HUGO_symb_entrez_uniprot = HUGO_symb_entrez_uniprot[~pd.isnull(HUGO_symb_entrez_uniprot['NCBI Gene ID(supplied by NCBI)'])]
		HUGO_symb_entrez_uniprot = HUGO_symb_entrez_uniprot.astype({'NCBI Gene ID(supplied by NCBI)': int}) 
		
		GRN_edges_df = pd.read_csv(proj_path + input_GRN + '_edges.csv', sep='\t', header=None)
		GRN_edges_df = GRN_edges_df.rename(columns={0: 'Source', 1: 'Target'})
		GRN_edges_df_entrez = pd.merge(pd.merge(GRN_edges_df, HUGO_symb_entrez_uniprot, left_on='Source', right_on='Approved symbol'), 
										HUGO_symb_entrez_uniprot, 
										left_on='Target', 
										right_on='Approved symbol')[['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y']]
		GRN_edges_df_entrez.columns = ['Source_Entrez_ID', 'Target_Entrez_ID']\
		
		return GRN_edges_df_entrez
		
	else:
	
		print("ID_format must either be 'Gene_Symbol' or 'Entrez'")

@st.cache(show_spinner=False, allow_output_mutation=True)
# Input must be in the form of PPI edgelist (tab-separated) with two columns for Gene A and Gene B, with no headers. 
def process_PPI_data(proj_path, input_PPI, input_filenames_dict, ID_format='Gene_Symbol'):

	print('Processing PPI edgelist...')
	
	if ID_format=='Entrez':
	
		PPI_edges_df_entrez = pd.read_csv(proj_path + input_PPI + '_edges.csv', sep='\t', header=None)
		PPI_edges_df_entrez = PPI_edges_df_entrez.rename(columns={0: 'Gene_A_Entrez_ID', 1: 'Gene_B_Entrez_ID'})
		
		return PPI_edges_df_entrez
		
	elif ID_format=='Gene_Symbol':
		
		HUGO_symb_entrez_uniprot = pd.read_csv(proj_path + input_filenames_dict['HUGO_symb_entrez_uniprot'], sep='\t')
		HUGO_symb_entrez_uniprot = HUGO_symb_entrez_uniprot[~pd.isnull(HUGO_symb_entrez_uniprot['NCBI Gene ID(supplied by NCBI)'])]
		HUGO_symb_entrez_uniprot = HUGO_symb_entrez_uniprot.astype({'NCBI Gene ID(supplied by NCBI)': int}) 	  
	
		PPI_edges_df = pd.read_csv(proj_path + input_PPI + '_edges.csv', sep='\t', header=None)
		PPI_edges_df = PPI_edges_df.rename(columns={0: 'Gene_A', 1: 'Gene_B'})	
		PPI_edges_df_entrez = pd.merge(pd.merge(PPI_edges_df, HUGO_symb_entrez_uniprot, left_on='Gene_A', right_on='Approved symbol'), 
										HUGO_symb_entrez_uniprot, 
										left_on='Gene_B', 
										right_on='Approved symbol')[['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y']]
		PPI_edges_df_entrez.columns = ['Gene_A_Entrez_ID', 'Gene_B_Entrez_ID']    
		
		return PPI_edges_df_entrez
		
	else:

		print("ID_format must either be 'Gene_Symbol' or 'Entrez'")

# preprocess the KEGG, PPI and GRN data
@st.cache(show_spinner=False, allow_output_mutation=True)
def process_data(proj_path, input_PPI, input_PPI_ID_format, input_GRN, input_GRN_ID_format, input_filenames_dict):
    
    ### import the most recent conversions from HUGO (more reliable than my gene-info-stripped since it's automatically updated)
    print('Loading Gene Symbol-Entrez ID-Uniprot ID mappings...')
    HUGO_symb_entrez_uniprot = pd.read_csv(proj_path + input_filenames_dict['HUGO_symb_entrez_uniprot'], sep='\t')
    # here we take Entrez IDs as the basis and only the part for which Entrez IDs are known (the NaN rows are usually the "withdrawn" entries)
    HUGO_symb_entrez_uniprot = HUGO_symb_entrez_uniprot[~pd.isnull(HUGO_symb_entrez_uniprot['NCBI Gene ID(supplied by NCBI)'])]
    HUGO_symb_entrez_uniprot = HUGO_symb_entrez_uniprot.astype({'NCBI Gene ID(supplied by NCBI)': int}) 
    
    ### import the protein-protein interaction network
    print('Loading PPI edges...')
    PPI_edges_df_entrez = process_PPI_data(proj_path, input_PPI, input_filenames_dict, ID_format=input_PPI_ID_format)
  
    ### import the GRN
    print('Loading GRN edges...')
    GRN_edges_df_entrez = process_GRN_data(proj_path, input_GRN, input_filenames_dict, ID_format=input_GRN_ID_format)

    ### import the KEGG signaling network, remove the GErel (gene regulatory) type
    print('Loading KEGG signaling network edges...')
    KEGG_all_nodes_df = pd.read_csv(proj_path + input_filenames_dict['KEGG_all_nodes_df'])
    KEGG_all_edges_df = pd.read_csv(proj_path + input_filenames_dict['KEGG_all_edges_df'])
    # Remove regulatory interactions
    KEGG_all_edges_df = KEGG_all_edges_df[KEGG_all_edges_df['Edge_type']!='GErel']
    

    
    ### make a master conversion file by taking the union of all three layers in terms of Entrez IDs
    print('Cross-mapping Gene Symbols, Entrez IDs and Uniprot IDs...')
    PPI_allnodes_entrez = set(PPI_edges_df_entrez['Gene_A_Entrez_ID']) | set(PPI_edges_df_entrez['Gene_B_Entrez_ID'])

    GRN_allnodes_entrez = set(GRN_edges_df_entrez['Source_Entrez_ID']) | set(GRN_edges_df_entrez['Target_Entrez_ID'])

    KEGG_allnodes_entrez = set([int(x.split('hsa:')[1]) for x in 
                                list(KEGG_all_nodes_df[KEGG_all_nodes_df['Node_type']=='gene']['KEGG_name(s)_expanded'].values)])

    all_layers_nodes_entrez = sorted(list(PPI_allnodes_entrez | GRN_allnodes_entrez | KEGG_allnodes_entrez))


    entrez_conversion_df = pd.DataFrame()
    entrez_conversion_df['Entrez ID'] = all_layers_nodes_entrez
    entrez_conversion_df = pd.merge(entrez_conversion_df, HUGO_symb_entrez_uniprot, 
                                    left_on='Entrez ID', right_on='NCBI Gene ID(supplied by NCBI)')[['Approved symbol', 
                                                                                                     'NCBI Gene ID(supplied by NCBI)']]
    entrez_conversion_df['KEGG_ID'] = 'hsa:' + entrez_conversion_df['NCBI Gene ID(supplied by NCBI)'].astype(str)
    
    ### Convert all edges to Entrez ID and build multilayer 
    KEGG_all_edges_entrez = pd.merge(pd.merge(KEGG_all_edges_df,  entrez_conversion_df, left_on='KEGG_name(s)_expanded_1', right_on='KEGG_ID'), 
                                     entrez_conversion_df, left_on='KEGG_name(s)_expanded_2', right_on='KEGG_ID')

    ### Add PPI edges to the KEGG signaling layer
    
    print('Combining KEGG and PPI edges...')
    # Designate it with the type "PPI". Add directed edges in both directions for all undirected edges.    
    PPI_all_edges_entrez = pd.DataFrame(columns=KEGG_all_edges_entrez.columns)
    PPI_all_edges_entrez['NCBI Gene ID(supplied by NCBI)_x'] = PPI_edges_df_entrez['Gene_A_Entrez_ID']
    PPI_all_edges_entrez['NCBI Gene ID(supplied by NCBI)_y'] = PPI_edges_df_entrez['Gene_B_Entrez_ID']
    PPI_all_edges_entrez['Path_label'] = 'ppi'
    PPI_all_edges_entrez['Edge_type'] = 'ppi'
    PPI_all_edges_entrez['Edge_subtype'] = 'ppi'

    KEGG_PPI_all_edges_entrez = pd.concat([PPI_all_edges_entrez, KEGG_all_edges_entrez])
    # remove the few edges with edge_subtype=NaN:
    KEGG_PPI_all_edges_entrez = KEGG_PPI_all_edges_entrez[~pd.isnull(KEGG_PPI_all_edges_entrez['Edge_subtype'])]
    KEGG_PPI_allnodes_entrez = (set(KEGG_PPI_all_edges_entrez['NCBI Gene ID(supplied by NCBI)_x']) |
                            set(KEGG_PPI_all_edges_entrez['NCBI Gene ID(supplied by NCBI)_y']))

    ### Get subgraph of GRN based on PPI+KEGG's nodes
    GRN_KEGG_PPI_edges = GRN_edges_df_entrez[(GRN_edges_df_entrez['Source_Entrez_ID'].isin(KEGG_PPI_allnodes_entrez)) & 
                                             (GRN_edges_df_entrez['Target_Entrez_ID'].isin(KEGG_PPI_allnodes_entrez))]
    
    ### Dataframe of Entrez IDs and consecutive numbers to be used as indices in the adjacency matrices 
    KEGG_PPI_allnodes_entrez_df = pd.DataFrame(sorted(list(KEGG_PPI_allnodes_entrez)))
    KEGG_PPI_allnodes_entrez_df['ix'] = pd.DataFrame(sorted(list(KEGG_PPI_allnodes_entrez))).index
    KEGG_PPI_allnodes_entrez_df = pd.merge(KEGG_PPI_allnodes_entrez_df, entrez_conversion_df, 
                                           left_on=0, right_on='NCBI Gene ID(supplied by NCBI)', how='left')[[0, 'ix','Approved symbol']]
    
    ### Edge and node dictionaries of KEGG pathways
    print('Generating edge and node dictionaries of KEGG pathways...')
    # Get all pathway names with the exception of "ppi", which means those edges just belong to the PPI and not a specific signaling pathway.
    KEGG_all_paths = sorted(list(set(sorted(KEGG_PPI_all_edges_entrez['Path_label'].unique())) - set(['ppi'])))

    
    print('Generating KEGG pathway-specific dataframes...')
    KEGG_path_edges_df_dict = {}
    KEGG_path_nodes_df_dict = {}
    KEGG_path_nodes_dict = {}
    KEGG_path_nodes_entrez_dict = {}

    for p in KEGG_all_paths:
        KEGG_path_edges_df_dict[p] = KEGG_PPI_all_edges_entrez[KEGG_PPI_all_edges_entrez['Path_label']==p]
        # remove the few edges with edge_subtype=NaN:
        KEGG_path_edges_df_dict[p] = KEGG_path_edges_df_dict[p][~pd.isnull(KEGG_path_edges_df_dict[p]['Edge_subtype'])]

        KEGG_path_nodes_dict[p] = set(KEGG_path_edges_df_dict[p]['KEGG_name(s)_expanded_1']) | \
                                                    set(KEGG_path_edges_df_dict[p]['KEGG_name(s)_expanded_2'])
        KEGG_path_nodes_entrez_dict[p] = set(KEGG_path_edges_df_dict[p]['NCBI Gene ID(supplied by NCBI)_x']) | \
                                                    set(KEGG_path_edges_df_dict[p]['NCBI Gene ID(supplied by NCBI)_y'])    
        KEGG_path_nodes_df_dict[p] = KEGG_all_nodes_df[(KEGG_all_nodes_df['Path_label']==p) & 
                                                       (KEGG_all_nodes_df['KEGG_name(s)_expanded'].isin(KEGG_path_nodes_dict[p]))]   
        
        
    
    KEGG_interaction_types_dict = {key: value for key, value in zip(KEGG_PPI_all_edges_entrez['Edge_subtype'].sort_values().unique()[1:],
                                                                    np.arange(1, len(KEGG_PPI_all_edges_entrez['Edge_subtype'].value_counts())))}
    all_motif_types = ['%s%s' %(i, j) for i, j in list(product(*[np.arange(len(KEGG_interaction_types_dict)+1), [-1, 0, 1]]))]
    all_motif_types_list = np.array([(i, j) for i, j in list(product(*[np.arange(len(KEGG_interaction_types_dict)+1), [-1, 0, 1]]))])  
    
    
    return (KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, PPI_all_edges_entrez, KEGG_PPI_all_edges_entrez, KEGG_PPI_allnodes_entrez_df, 
            KEGG_interaction_types_dict, KEGG_all_edges_entrez, KEGG_all_paths, KEGG_path_nodes_entrez_dict, 
            all_motif_types, all_motif_types_list)

# generate sparse arrays of the KEGG, PPI and GRN networks for memory-efficient use in downstream calculations
@st.cache(show_spinner=False, allow_output_mutation=True)
def sparse_layers(KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, PPI_all_edges_entrez, KEGG_all_edges_entrez, KEGG_PPI_all_edges_entrez, 
                  KEGG_interaction_types_dict):
    
    # remove self edges of all sorts while filling the dataframes (these are mostly in PPI, but KEGG layer also has some, e.g. "state change" edges) 
    # so that the motif counts are symmetric as well.
    
    print('Generating sparse array for the  GRN...')
    A_GRN = pd.DataFrame(0, index=sorted(list(KEGG_PPI_allnodes_entrez)), columns=sorted(list(KEGG_PPI_allnodes_entrez)))
    for s, t in tqdm(GRN_KEGG_PPI_edges[['Source_Entrez_ID', 'Target_Entrez_ID']].drop_duplicates().values, 
                     position=0, leave=True):
        if s != t:
            A_GRN.at[s, t] = 1
            A_GRN.at[t, s] = -1
    
    print('Generating sparse array for the  PPI...')
    A_PPI = pd.DataFrame(0, index=sorted(list(KEGG_PPI_allnodes_entrez)), columns=sorted(list(KEGG_PPI_allnodes_entrez)))
    for s, t in tqdm(PPI_all_edges_entrez[['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y']].values, position=0, leave=True):
        if s != t:
            A_PPI.at[s, t] = 9
            A_PPI.at[t, s] = 9

    A_GRN_sparr = sparse.csr_matrix(A_GRN.to_numpy())
    A_PPI_sparr = sparse.csr_matrix(A_PPI.to_numpy())
    
    print('Generating sparse array for KEGG and PPI combined...')    
    # This is the collapsed KEGG+PPI network that includes all KEGG edges and PPI edges (for the counting of 0X motifs) and does not differentiate 
    # between the different edge types in the KEGG/PPI layer. Everything is coded '13'
    A_KEGGPPI = pd.DataFrame(0, index=sorted(list(KEGG_PPI_allnodes_entrez)), columns=sorted(list(KEGG_PPI_allnodes_entrez)))
    for s, t in tqdm(KEGG_PPI_all_edges_entrez[['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y']].drop_duplicates().values, 
                     position=0, leave=True):
        if s != t:
            A_KEGGPPI.at[s, t] = 13
            A_KEGGPPI.at[t, s] = 13

    A_KEGGPPI_sparr = sparse.csr_matrix(A_KEGGPPI.to_numpy()) 

    print('Generating edge type-specific sparse arrays for KEGG...')     
    A_KEGG_e_sparr_dict = {}
    for e in tqdm(KEGG_interaction_types_dict.keys(), position=0, leave=True):    
        KEGG_e_edges_entrez = KEGG_all_edges_entrez[KEGG_all_edges_entrez['Edge_subtype']==e]

        # dataframe for the KEGG signaling pathways layer with the edge type e: 
        # 0: no edge, 1-11 as described above in "KEGG_interaction_types_dict"
        A_KEGG_e = A_PPI.copy()
        for s, t in KEGG_e_edges_entrez[['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y']].values:
            if (e == 'Within_group') | (e == 'binding/association'):
                A_KEGG_e.at[s, t] = KEGG_interaction_types_dict['binding/association']
                A_KEGG_e.at[t, s] = KEGG_interaction_types_dict['binding/association']
            else:    
                A_KEGG_e.at[s, t] = KEGG_interaction_types_dict[e]
                A_KEGG_e.at[t, s] = KEGG_interaction_types_dict[e] * -1

        A_KEGG_e_sparr_dict[e] = sparse.csr_matrix(A_KEGG_e.to_numpy())
        
        
    return (A_GRN_sparr, A_PPI_sparr, A_KEGGPPI_sparr, A_KEGG_e_sparr_dict)


@st.cache(show_spinner=False, allow_output_mutation=True, hash_funcs={'_json.Scanner': hash})
def return_between_multilink_edges(p1, p2, multilink_type, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, A_GRN_sparr, 
                                  A_KEGGPPI_sparr, KEGG_interaction_types_dict, A_KEGG_e_sparr_dict):

    multilink_counts = {}   

    path_ix1 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p1]))), how='inner', 
                   left_on=0, right_on=0)['ix'].values
    path_ix2 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p2]))), how='inner', 
                   left_on=0, right_on=0)['ix'].values

    A_2_arr = A_GRN_sparr[tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray() 

    if '-' in multilink_type:
        sig_ix = int(multilink_type[:len(multilink_type)-2])
        grn_ix = int(multilink_type[-2:])
    else:
        sig_ix = int(multilink_type[:len(multilink_type)-1])
        grn_ix = int(multilink_type[-1:])  

    if multilink_type in ['0-1', '00', '01']:       
        A_1_arr = A_KEGGPPI_sparr[tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()
        A_mux = np.dstack([A_1_arr, A_2_arr])

    else:
        e = list(KEGG_interaction_types_dict.keys())[list(KEGG_interaction_types_dict.values()).index(sig_ix)]
        A_1_arr = A_KEGG_e_sparr_dict[e][tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()
        A_mux = np.dstack([A_1_arr, A_2_arr])


    sig_edges_ix = [(i, j) for i, j in zip(np.setdiff1d(path_ix1, path_ix2)[np.argwhere((A_mux == np.array([sig_ix,  grn_ix])).all(-1))[:, 0]],
                                           np.setdiff1d(path_ix2, path_ix1)[np.argwhere((A_mux == np.array([sig_ix,  grn_ix])).all(-1))[:, 1]])]
    between_multilink_edges_df = pd.DataFrame(index=np.arange(len(sig_edges_ix)), 
                                              columns=['Node 1 ix', 'Node 2 ix', 'Node 1 Gene Symbol', 'Node 2 Gene Symbol', 
                                                       'Node 1 Entrez ID', 'Node 2 Entrez ID'])
    
    for n, (i, j) in enumerate(sig_edges_ix):
        between_multilink_edges_df.at[n, 'Node 1 ix'] = i
        between_multilink_edges_df.at[n, 'Node 2 ix'] = j
        between_multilink_edges_df.at[n, 'Node 1 Gene Symbol'] = KEGG_PPI_allnodes_entrez_df[KEGG_PPI_allnodes_entrez_df['ix']==i]['Approved symbol'].values[0]
        between_multilink_edges_df.at[n, 'Node 2 Gene Symbol'] = KEGG_PPI_allnodes_entrez_df[KEGG_PPI_allnodes_entrez_df['ix']==j]['Approved symbol'].values[0]
        between_multilink_edges_df.at[n, 'Node 1 Entrez ID'] = KEGG_PPI_allnodes_entrez_df[KEGG_PPI_allnodes_entrez_df['ix']==i][0].values[0]
        between_multilink_edges_df.at[n, 'Node 2 Entrez ID'] = KEGG_PPI_allnodes_entrez_df[KEGG_PPI_allnodes_entrez_df['ix']==j][0].values[0]

    return between_multilink_edges_df


@st.cache(show_spinner=False, allow_output_mutation=True, hash_funcs={'_json.Scanner': hash})
def return_shortest_multilink_edges_forStreamlit(proj_path, sp_threshold, input_PPI, p1, p2, multilink_type, KEGG_PPI_allnodes_entrez_df, 
												KEGG_path_nodes_entrez_dict, A_GRN_sparr, 
                                  				A_KEGGPPI_sparr, KEGG_interaction_types_dict, A_KEGG_e_sparr_dict):

    f = open(proj_path + '%s_sp%s_discovery_shortest_path_dicts/%s.json' % (input_PPI, sp_threshold, p1))
    temp_dict = json.load(f)
    
    multilink_counts = {}   

    A_2_arr = np.asarray(A_GRN_sparr[np.array(temp_dict[p2])[:, 0], 
                                     np.array(temp_dict[p2])[:, 1]])[0]

    if '-' in multilink_type:
        sig_ix = int(multilink_type[:len(multilink_type)-2])
        grn_ix = int(multilink_type[-2:])
    else:
        sig_ix = int(multilink_type[:len(multilink_type)-1])
        grn_ix = int(multilink_type[-1:])  

    if multilink_type in ['0-1', '00', '01']:       
        A_1_arr = np.asarray(A_KEGGPPI_sparr[np.array(temp_dict[p2])[:, 0], 
                                             np.array(temp_dict[p2])[:, 1]])[0]
        A_mux = np.dstack([A_1_arr, A_2_arr])

    else:
        e = list(KEGG_interaction_types_dict.keys())[list(KEGG_interaction_types_dict.values()).index(sig_ix)]
        A_1_arr = np.asarray(A_KEGG_e_sparr_dict[e][np.array(temp_dict[p2])[:, 0], 
                                                    np.array(temp_dict[p2])[:, 1]])[0]
        A_mux = np.dstack([A_1_arr, A_2_arr])


    sig_edges_ix = np.array(temp_dict[p2])[np.argwhere((A_mux[0] == np.array([sig_ix,  grn_ix])).all(-1)).flatten()]    
    shortest_multilink_edges_df = pd.DataFrame(index=np.arange(len(sig_edges_ix)), 
                                               columns=['Node 1 ix', 'Node 2 ix', 'Node 1 Gene Symbol', 'Node 2 Gene Symbol', 
                                                       'Node 1 Entrez ID', 'Node 2 Entrez ID'])
    
    for n, (i, j) in enumerate(sig_edges_ix):
        shortest_multilink_edges_df.at[n, 'Node 1 ix'] = i
        shortest_multilink_edges_df.at[n, 'Node 2 ix'] = j
        shortest_multilink_edges_df.at[n, 'Node 1 Gene Symbol'] = KEGG_PPI_allnodes_entrez_df[KEGG_PPI_allnodes_entrez_df['ix']==i]['Approved symbol'].values[0]
        shortest_multilink_edges_df.at[n, 'Node 2 Gene Symbol'] = KEGG_PPI_allnodes_entrez_df[KEGG_PPI_allnodes_entrez_df['ix']==j]['Approved symbol'].values[0]
        shortest_multilink_edges_df.at[n, 'Node 1 Entrez ID'] = KEGG_PPI_allnodes_entrez_df[KEGG_PPI_allnodes_entrez_df['ix']==i][0].values[0]
        shortest_multilink_edges_df.at[n, 'Node 2 Entrez ID'] = KEGG_PPI_allnodes_entrez_df[KEGG_PPI_allnodes_entrez_df['ix']==j][0].values[0]

    return shortest_multilink_edges_df         
        
