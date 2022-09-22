import sys
import numpy as np
import pandas as pd
import random
from itertools import product
from itertools import permutations
from tqdm import tqdm
import pickle
from scipy import sparse
from numba import jit
import os.path
import networkx as nx
from sklearn.metrics import auc
import scipy.stats as st
import statsmodels.stats.multitest as multi
from itertools import islice
import matplotlib.pyplot as plt


# This version is made to incorporate the opposite signs in the KEGG signaling and GRN networks. It is meant for symmetric adj. matrices like the 
# above, but after the rewiring is done, replaces the opposite matrix element with the exact previous corresponding value, i.e. AB_val and BA_val 
# (and similarly CD_val and DC_val are separately kept track of).
@jit(nopython=True)
def edge_swap_undir_sign(adj, ntry):
    
    adj_rand = np.copy(adj)
    nrew = 0
    ix1 = np.where(adj_rand>0)[0]
    ix2 = np.where(adj_rand>0)[1]
    aux = np.where(ix1 > ix2)
    ix1 = ix1[aux]
    ix2 = ix2[aux]

    # swap ntry * # of edges (upper triangle)
    for i in np.arange(ntry*len(ix1)):
        
        # choose two edges at random
        swap_ix1 = int(random.random() * len(ix1)) #random.sample(set(np.arange(len(ix1))), 1)[0]
        swap_ix2 = int(random.random() * len(ix1)) #random.sample(set(np.arange(len(ix1))), 1)[0]
        A = ix1[swap_ix1]
        B = ix2[swap_ix1]
        C = ix1[swap_ix2]
        D = ix2[swap_ix2]

        if (len(set([A, B, C, D])) == 4): #(A != C) & (A != D) & (B != D) & (B != C)
            if random.random() > 0.5:
                if ((adj_rand[A, C] == 0) & (adj_rand[B, D] == 0)):
                    
                    AB_val = adj_rand[A, B]
                    CD_val = adj_rand[C, D]
                    BA_val = adj_rand[B, A]
                    DC_val = adj_rand[D, C]
                    
                    adj_rand[A, B] = 0
                    adj_rand[C, D] = 0
                    adj_rand[B, A] = 0
                    adj_rand[D, C] = 0
                    adj_rand[A, C] = AB_val
                    adj_rand[B, D] = CD_val
                    adj_rand[C, A] = BA_val
                    adj_rand[D, B] = DC_val
                    nrew = nrew + 1
                    ix1[swap_ix1] = A
                    ix2[swap_ix1] = C
                    ix1[swap_ix2] = B
                    ix2[swap_ix2] = D 
            
            else:
                E = C
                C = D
                D = E
                #del E
                if ((adj_rand[A, C] == 0) & (adj_rand[B, D] == 0)):

                    AB_val = adj_rand[A, B]
                    CD_val = adj_rand[C, D]
                    BA_val = adj_rand[B, A]
                    DC_val = adj_rand[D, C]
                    
                    adj_rand[A, B] = 0
                    adj_rand[C, D] = 0
                    adj_rand[B, A] = 0
                    adj_rand[D, C] = 0
                    adj_rand[A, C] = AB_val
                    adj_rand[B, D] = CD_val
                    adj_rand[C, A] = BA_val
                    adj_rand[D, B] = DC_val
                    nrew = nrew + 1
                    ix1[swap_ix1] = A
                    ix2[swap_ix1] = C
                    ix1[swap_ix2] = B
                    ix2[swap_ix2] = D 
                                    
    return (adj_rand, nrew)      

# Input must be in the form of GRN edgelist (tab-separated) with two columns for source (TF) and target, respectively, and with no headers. 
def process_GRN_data(proj_path, input_GRN):
    
    print('Processing GRN edgelist...')
    GRN_edges_df = pd.read_csv(proj_path + input_GRN + '_edges.csv', sep='\t', header=None)
    GRN_edges_df = GRN_edges_df.rename(columns={0: 'Source', 1: 'Target'})
    
    return GRN_edges_df

# Function to convert the human GRN refseq identifiers at different p-value cutoffs to the Source-Target format to be read by the 
# "process_GRN_data" function.
def process_GRN_cisbp_data(proj_path, input_refseq,  output_GRN):
    
    if not os.path.exists(proj_path + output_GRN + '_edges.csv'):

        print('Creating GRN edgelist by mapping RefSeq identifiers...')
        df_motinf = pd.read_csv(proj_path + input_filenames_dict['df_motinf'], sep='\t')
        df_motinf['motif_id'] = df_motinf['V1'].str.split('_', expand=True)[0]

        motifDict = {motifID.split('_')[0]: gene for motifID, gene in zip(df_motinf['V1'].values, df_motinf['V2'].values)}

        df_edges = pd.read_csv(proj_path + input_refseq, sep='\t', header=None)
        GRN_pvals = pd.merge(df_edges, df_motinf, left_on=0, right_on='motif_id', how='left').sort_values([0, 1])
        GRN_edges_df = GRN_pvals[['V2', 1]].rename(columns={'V2': 'Source', 1: 'Target'})

        GRN_edges_df.to_csv(proj_path + output_GRN + '_edges.csv', sep='\t', index=False, header=False)

# preprocess the KEGG, PPI and GRN data
def process_data(proj_path, input_GRN, input_filenames_dict):
    
    ### import the most recent conversions from HUGO (more reliable than my gene-info-stripped since it's automatically updated)
    print('Loading Gene Symbol-Entrez ID-Uniprot ID mappings...')
    HUGO_symb_entrez_uniprot = pd.read_csv(proj_path + input_filenames_dict['HUGO_symb_entrez_uniprot'], sep='\t')
    # here we take Entrez IDs as the basis and only the part for which Entrez IDs are known (the NaN rows are usually the "withdrawn" entries)
    HUGO_symb_entrez_uniprot = HUGO_symb_entrez_uniprot[~pd.isnull(HUGO_symb_entrez_uniprot['NCBI Gene ID(supplied by NCBI)'])]
    HUGO_symb_entrez_uniprot = HUGO_symb_entrez_uniprot.astype({'NCBI Gene ID(supplied by NCBI)': int}) 
    
    ### import the protein-protein interaction network from: Cheng, Feixiong, István A. Kovács, and Albert-László Barabási. 
    # "Network-based prediction of drug combinations." Nature communications 10.1 (2019): 1-11.
    print('Loading the PPI edges from Cheng et al. Nature communications 10.1 (2019): 1-11...')
    PPI_Cheng_2019_data = pd.read_csv(proj_path + input_filenames_dict['PPI_Cheng_2019_data'])

    PPI_Cheng_2019_data_symb = pd.merge(pd.merge(PPI_Cheng_2019_data, HUGO_symb_entrez_uniprot, 
                                                 left_on='Protein_A_Entrez_ID', 
                                                 right_on='NCBI Gene ID(supplied by NCBI)'), 
                                        HUGO_symb_entrez_uniprot, 
                                        left_on='Protein_B_Entrez_ID', 
                                        right_on='NCBI Gene ID(supplied by NCBI)')[['Approved symbol_x', 'Approved symbol_y', 'data_source(s)']]
    PPI_Cheng_2019_data_symb.columns = ['Protein_A_Gene_Symbol', 'Protein_B_Gene_Symbol', 'data_source(s)']
  
    ### import the GRN
    print('Loading GRN edges...')
    GRN_edges_df = process_GRN_data(proj_path, input_GRN)

    ### import the KEGG signaling network, remove the GErel (gene regulatory) type
    print('Loading KEGG signaling network edges...')
    KEGG_all_nodes_df = pd.read_csv(proj_path + input_filenames_dict['KEGG_all_nodes_df'])
    KEGG_all_edges_df = pd.read_csv(proj_path + input_filenames_dict['KEGG_all_edges_df'])
    # Remove regulatory interactions
    KEGG_all_edges_df = KEGG_all_edges_df[KEGG_all_edges_df['Edge_type']!='GErel']
    

    
    ### make a master conversion file by taking the union of all three layers in terms of Entrez IDs
    print('Cross-mapping Gene Symbols, Entrez IDs and Uniprot IDs...')
    PPI_Cheng_2019_allnodes_entrez = set(PPI_Cheng_2019_data['Protein_A_Entrez_ID']) | set(PPI_Cheng_2019_data['Protein_B_Entrez_ID'])

    GRN_allnodes_symb = set(GRN_edges_df['Source']) | set(GRN_edges_df['Target'])
    GRN_allnodes_symb_df = pd.DataFrame(sorted(list(GRN_allnodes_symb)))
    GRN_allnodes_entrez = set(pd.merge(GRN_allnodes_symb_df, HUGO_symb_entrez_uniprot, 
                                       left_on=0, right_on='Approved symbol')['NCBI Gene ID(supplied by NCBI)'])

    KEGG_allnodes_entrez = set([int(x.split('hsa:')[1]) for x in 
                                list(KEGG_all_nodes_df[KEGG_all_nodes_df['Node_type']=='gene']['KEGG_name(s)_expanded'].values)])

    all_layers_nodes_entrez = sorted(list(PPI_Cheng_2019_allnodes_entrez | GRN_allnodes_entrez | KEGG_allnodes_entrez))


    entrez_conversion_df = pd.DataFrame()
    entrez_conversion_df['Entrez ID'] = all_layers_nodes_entrez
    entrez_conversion_df = pd.merge(entrez_conversion_df, HUGO_symb_entrez_uniprot, 
                                    left_on='Entrez ID', right_on='NCBI Gene ID(supplied by NCBI)')[['Approved symbol', 
                                                                                                     'NCBI Gene ID(supplied by NCBI)']]
    entrez_conversion_df['KEGG_ID'] = 'hsa:' + entrez_conversion_df['NCBI Gene ID(supplied by NCBI)'].astype(str)
    
    ### Convert all edges to Entrez ID and build multilayer (PPI already in Entrez in "PPI_Cheng_2019_data")
    KEGG_all_edges_entrez = pd.merge(pd.merge(KEGG_all_edges_df,  entrez_conversion_df, left_on='KEGG_name(s)_expanded_1', right_on='KEGG_ID'), 
                                     entrez_conversion_df, left_on='KEGG_name(s)_expanded_2', right_on='KEGG_ID')
    GRN_edges_entrez = pd.merge(pd.merge(GRN_edges_df, entrez_conversion_df, left_on='Source', right_on='Approved symbol'), 
                                entrez_conversion_df, left_on='Target', right_on='Approved symbol')[['NCBI Gene ID(supplied by NCBI)_x', 
                                                                                                    'NCBI Gene ID(supplied by NCBI)_y']]

    ### Add PPI edges to the KEGG signaling layer
    
    print('Combining KEGG and PPI edges...')
    # Designate it with the type "PPI". Add directed edges in both directions for all undirected edges.    
    PPI_all_edges_entrez = pd.DataFrame(columns=KEGG_all_edges_entrez.columns)
    PPI_all_edges_entrez['NCBI Gene ID(supplied by NCBI)_x'] = PPI_Cheng_2019_data['Protein_A_Entrez_ID']
    PPI_all_edges_entrez['NCBI Gene ID(supplied by NCBI)_y'] = PPI_Cheng_2019_data['Protein_B_Entrez_ID']
    PPI_all_edges_entrez['Path_label'] = 'ppi'
    PPI_all_edges_entrez['Edge_type'] = 'ppi'
    PPI_all_edges_entrez['Edge_subtype'] = 'ppi'

    KEGG_PPI_all_edges_entrez = pd.concat([PPI_all_edges_entrez, KEGG_all_edges_entrez])
    # remove the few edges with edge_subtype=NaN:
    KEGG_PPI_all_edges_entrez = KEGG_PPI_all_edges_entrez[~pd.isnull(KEGG_PPI_all_edges_entrez['Edge_subtype'])]
    KEGG_PPI_allnodes_entrez = (set(KEGG_PPI_all_edges_entrez['NCBI Gene ID(supplied by NCBI)_x']) |
                            set(KEGG_PPI_all_edges_entrez['NCBI Gene ID(supplied by NCBI)_y']))

    ### Get subgraph of GRN based on PPI+KEGG's nodes
    GRN_KEGG_PPI_edges = GRN_edges_entrez[(GRN_edges_entrez['NCBI Gene ID(supplied by NCBI)_x'].isin(KEGG_PPI_allnodes_entrez)) & 
                                                    (GRN_edges_entrez['NCBI Gene ID(supplied by NCBI)_y'].isin(KEGG_PPI_allnodes_entrez))]
    
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
def sparse_layers(KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, PPI_all_edges_entrez, KEGG_all_edges_entrez, KEGG_PPI_all_edges_entrez, 
                  KEGG_interaction_types_dict):
    
    # remove self edges of all sorts while filling the dataframes (these are mostly in PPI, but KEGG layer also has some, e.g. "state change" edges) 
    # so that the motif counts are symmetric as well.
    
    print('Generating sparse array for the  GRN...')
    A_GRN = pd.DataFrame(0, index=sorted(list(KEGG_PPI_allnodes_entrez)), columns=sorted(list(KEGG_PPI_allnodes_entrez)))
    for s, t in tqdm(GRN_KEGG_PPI_edges[['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y']].drop_duplicates().values, 
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

def XTalk_DB(KEGG_all_paths, proj_path, input_filenames_dict):
    
    print('Processing XTalkDB data for benchmarking...')
    XTalk_DB = pd.read_csv(proj_path + input_filenames_dict['XTalk_DB'])
    
    # Two names do not match exactly with the latest version of the KEGG database: 'Insulin secretion pathway' --> 'Insulin secretion'
    # and 'Jak-STAT signaling pathway' --> 'JAK-STAT signaling pathway'.
    XTalk_DB['Pathway A'] = XTalk_DB['Pathway A'].replace('Insulin secretion pathway', 'Insulin secretion')
    XTalk_DB['Pathway B'] = XTalk_DB['Pathway B'].replace('Insulin secretion pathway', 'Insulin secretion')

    XTalk_DB['Pathway A'] = XTalk_DB['Pathway A'].replace('Jak-STAT signaling pathway', 'JAK-STAT signaling pathway')
    XTalk_DB['Pathway B'] = XTalk_DB['Pathway B'].replace('Jak-STAT signaling pathway', 'JAK-STAT signaling pathway')
    
    # set of pathways common between our 32 pathways and the 26 pathways in XTalkDB
    common_crosstalk_paths = set(KEGG_all_paths) & set(XTalk_DB['Pathway A'])

    # pathway pairs from the common set of pathways
    common_perms_str = ['%s-%s'%(i, j) for i, j in permutations(common_crosstalk_paths, 2)]
    common_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]    

    XTalk_DB_positives = XTalk_DB[XTalk_DB['Crosstalk']=='yes'][['Pathway A', 'Pathway B']].drop_duplicates()

    # true positive common cross-talking pathways
    XTalk_DB_positives_common = XTalk_DB_positives[XTalk_DB_positives['Pathway A'].isin(common_crosstalk_paths) & 
                                                   XTalk_DB_positives['Pathway B'].isin(common_crosstalk_paths)]

    # make pathways into strings separated by dash to match the format of the multilink results
    XTalk_DB_positives_common_str = XTalk_DB_positives_common['Pathway A'].str.cat(XTalk_DB_positives_common['Pathway B'], sep='-').values
    
    return (common_crosstalk_paths, XTalk_DB_positives_common, XTalk_DB_positives_common_str)

def randomize_GRN(A_GRN_sparr, proj_path, input_GRN, N_runs=500, N_swap=10, return_output=False):
    
    if not os.path.exists(proj_path + input_GRN + '_A_GRN_sparr_rand_dict_%sruns.pickle' % N_runs):
    
        print('Computing ensemble of randomized GRN sparse  matrices...')
        A_GRN_sparr_rand_dict = {}

        for nrand in tqdm(np.arange(N_runs), position=0, leave=True):
            A_GRN_arr_rand, A_GRN_nrew_rand = edge_swap_undir_sign(A_GRN_sparr.toarray(), N_swap)
            A_GRN_sparr_rand_dict[nrand] = sparse.csr_matrix(A_GRN_arr_rand)
   
        with open(proj_path + input_GRN + '_A_GRN_sparr_rand_dict_%sruns.pickle' % N_runs, 'wb') as fp:
            pickle.dump(A_GRN_sparr_rand_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
        
        if return_output==True:        
            return A_GRN_sparr_rand_dict
        
    else:
        
        print('Reading ensemble of randomized GRN sparse  matrices...')
        with open(proj_path + input_GRN + '_A_GRN_sparr_rand_dict_500runs.pickle', 'rb') as fp:
            A_GRN_sparr_rand_dict = pickle.load(fp)        

        if return_output==True:        
            return A_GRN_sparr_rand_dict

def randomize_KEGG_e(A_KEGG_e_sparr_dict, KEGG_interaction_types_dict, proj_path, input_GRN, N_runs=500, N_swap=10, return_output=False):

    if not os.path.exists(proj_path + 'A_KEGG_e_sparr_rand_dict_%sruns.pickle' % N_runs):    
    
        print('Computing ensemble of randomized edge type-specific KEGG sparse  matrices...')
        A_KEGG_e_sparr_rand_dict = {i: {} for i in np.arange(N_runs)}

        for nrand in tqdm(np.arange(N_runs), position=0, leave=True):
            for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):
                A_KEGG_e_arr_rand, A_KEGG_e_nrew_rand = edge_swap_undir_sign(A_KEGG_e_sparr_dict[e].toarray(), N_swap)
                A_KEGG_e_sparr_rand_dict[nrand][e] = sparse.csr_matrix(A_KEGG_e_arr_rand)   
         
        with open(proj_path + 'A_KEGG_e_sparr_rand_dict_%sruns.pickle' % N_runs, 'wb') as fp:
            pickle.dump(A_KEGG_e_sparr_rand_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
        
        if return_output==True:
            return A_KEGG_e_sparr_rand_dict
        
    else:
        
        print('Reading ensemble of randomized edge type-specific KEGG sparse  matrices... (This may take 5-10 mins.)')
        with open(proj_path + 'A_KEGG_e_sparr_rand_dict_%sruns.pickle' % N_runs, 'rb') as fp:
            A_KEGG_e_sparr_rand_dict = pickle.load(fp)        
        
        if return_output==True:
            return A_KEGG_e_sparr_rand_dict        

def randomize_KEGGPPI(KEGG_PPI_allnodes_entrez, KEGG_interaction_types_dict, A_KEGG_e_sparr_rand_dict_dir, proj_path, input_GRN, 
                      return_output=False):
    
    if not os.path.exists(proj_path + 'A_KEGGPPI_sparr_rand_dict_%sruns.pickle' % A_KEGG_e_sparr_rand_dict_dir.split('runs')[0].split('_')[-1]):  
        
        print('Reading ensemble of randomized edge type-specific KEGG sparse  matrices...')
        with open(A_KEGG_e_sparr_rand_dict_dir, 'rb') as fp:
            A_KEGG_e_sparr_rand_dict = pickle.load(fp)  
    
        print('Computing ensemble of randomized KEGG+PPI combined sparse  matrices...')
        A_KEGGPPI_sparr_rand_dict = {}

        for nrand in tqdm(np.arange(len(A_KEGG_e_sparr_rand_dict)), position=0, leave=True):

            temp_sparse = sparse.csr_matrix((len(KEGG_PPI_allnodes_entrez), len(KEGG_PPI_allnodes_entrez)))

            # get only the KEGG interaction types from the A_KEGG_e adjacency matrices (which include both the KEGG interaction type and ppi (9))
            for e in sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))):
                temp_sparse = temp_sparse + (abs(A_KEGG_e_sparr_rand_dict[nrand][e])==KEGG_interaction_types_dict[e]).astype(int)

            # add at the end the ppi type
            temp_sparse = temp_sparse + (abs(A_KEGG_e_sparr_rand_dict[nrand]['ppi'])==KEGG_interaction_types_dict['ppi']).astype(int)      

            A_KEGGPPI_sparr_rand_dict[nrand] = (temp_sparse!=0).astype(int)
  
        with open(proj_path + 'A_KEGGPPI_sparr_rand_dict_%sruns.pickle' % len(A_KEGG_e_sparr_rand_dict), 'wb') as fp:
            pickle.dump(A_KEGGPPI_sparr_rand_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)    
        
        if return_output==True:
            return A_KEGGPPI_sparr_rand_dict
        
    else:
        
        print('Reading ensemble of randomized KEGG+PPI combined sparse  matrices...')
        with open(proj_path + 'A_KEGGPPI_sparr_rand_dict_%sruns.pickle' % A_KEGG_e_sparr_rand_dict_dir.split('runs')[0].split('_')[-1], 'rb') as fp:
            A_KEGGPPI_sparr_rand_dict = pickle.load(fp)
        
        if return_output==True:
            return A_KEGGPPI_sparr_rand_dict     

def overall_edge_overlap(A_GRN_sparr, A_KEGGPPI_sparr, A_GRN_sparr_rand_dict, A_KEGGPPI_sparr_rand_dict, proj_path, input_GRN):
    
    mat_product = np.multiply(A_GRN_sparr.toarray(), A_KEGGPPI_sparr.toarray())
    overall_overlap_edge_counts = len(np.triu(mat_product)[np.triu(mat_product)!=0])                   
                         
    if os.path.exists(proj_path + input_GRN + '_overall_overlap_edge_counts_rand.csv'):
        print('Reading the overlapping edges between the signaling and regulatory layer for the entire multilayer network...')
        overall_overlap_edge_counts_rand = pd.read_csv(proj_path + input_GRN + '_overall_overlap_edge_counts_rand.csv', index_col=0) 
    else:
        print('Calculating the overlapping edges between the signaling and regulatory layer for the entire multilayer network...')
        overall_overlap_edge_counts_rand = np.zeros(len(A_GRN_sparr_rand_dict))
        for i in tqdm(np.arange(len(A_GRN_sparr_rand_dict)), position=0, leave=True):
            mat_product = np.multiply(A_GRN_sparr_rand_dict[i].toarray(), A_KEGGPPI_sparr_rand_dict[i].toarray())
            overall_overlap_edge_counts_rand[i] = len(np.triu(mat_product)[np.triu(mat_product)!=0])               
        pd.DataFrame(overall_overlap_edge_counts_rand).to_csv(proj_path + input_GRN + '_overall_overlap_edge_counts_rand.csv')
    
    overall_overlap_zscore = (overall_overlap_edge_counts - np.mean(overall_overlap_edge_counts_rand)) / np.std(overall_overlap_edge_counts_rand)
    overall_overlap_emp_pval = (overall_overlap_edge_counts_rand >= overall_overlap_edge_counts).sum() / len(overall_overlap_edge_counts_rand)
                         
    return (overall_overlap_edge_counts, overall_overlap_edge_counts_rand, overall_overlap_zscore, overall_overlap_emp_pval)

def pathway_edge_overlap(KEGG_all_paths, A_GRN_sparr_rand_dict, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, A_KEGGPPI_sparr, 
                         A_GRN_sparr, A_KEGGPPI_sparr_rand_dict, proj_path, input_GRN):
    
    if (os.path.exists(proj_path + input_GRN + '_pathway_overlap_edge_stats.csv') & 
        os.path.exists(proj_path + input_GRN + '_pathway_overlap_edge_counts_rand.csv')):
        print('Reading the overlapping edges between the signaling and regulatory layer within each KEGG pathway...')
        pathway_overlap_edge_stats = pd.read_csv(proj_path + input_GRN + '_pathway_overlap_edge_stats.csv', index_col=0)     
        pathway_overlap_edge_counts_rand = pd.read_csv(proj_path + input_GRN + '_pathway_overlap_edge_counts_rand.csv', index_col=0)   
    
    else:
        print('Calculating the overlapping edges between the signaling and regulatory layer within each KEGG pathway...')
        pathway_overlap_edge_stats = pd.DataFrame(index=KEGG_all_paths, columns=['Edge overlap', 'emp. pval', 'z-score'])  
        pathway_overlap_edge_counts_rand = pd.DataFrame(index=KEGG_all_paths, columns=np.arange(len(A_GRN_sparr_rand_dict)))

        for p in tqdm(KEGG_all_paths, position=0, leave=True): 

            path_ix = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values

            mat_product = np.multiply(A_KEGGPPI_sparr[tuple([path_ix])][:, path_ix].toarray(), A_GRN_sparr[tuple([path_ix])][:, path_ix].toarray())    
            pathway_overlap_edge_stats.at[p, 'Edge overlap'] = len(np.triu(mat_product)[np.triu(mat_product)!=0])

            for i in np.arange(len(A_GRN_sparr_rand_dict)):

                mat_product = np.multiply(A_KEGGPPI_sparr_rand_dict[i][tuple([path_ix])][:, path_ix].toarray(), 
                                       A_GRN_sparr_rand_dict[i][tuple([path_ix])][:, path_ix].toarray())    
                pathway_overlap_edge_counts_rand.at[p, i] = len(np.triu(mat_product)[np.triu(mat_product)!=0])  

            pathway_overlap_edge_stats.at[p, 'z-score'] = ((pathway_overlap_edge_stats.at[p, 'Edge overlap'] - 
                                                            pathway_overlap_edge_counts_rand.loc[p].mean()) /
                                                           pathway_overlap_edge_counts_rand.loc[p].std())

            if pathway_overlap_edge_stats.at[p, 'z-score'] > 0:
                pathway_overlap_edge_stats.at[p, 'emp. pval'] = ((pathway_overlap_edge_stats.at[p, 'Edge overlap'] <= 
                                                                  pathway_overlap_edge_counts_rand.loc[p]).sum() / 
                                                                 pathway_overlap_edge_counts_rand.shape[1])
            elif pathway_overlap_edge_stats.at[p, 'z-score'] < 0:
                pathway_overlap_edge_stats.at[p, 'emp. pval'] = ((pathway_overlap_edge_stats.at[p, 'Edge overlap'] >= 
                                                                  pathway_overlap_edge_counts_rand.loc[p]).sum() /
                                                                 pathway_overlap_edge_counts_rand.shape[1])


        pathway_overlap_edge_counts_rand.to_csv(proj_path + input_GRN + '_pathway_overlap_edge_counts_rand.csv')
        pathway_overlap_edge_stats.to_csv(proj_path + input_GRN + '_pathway_overlap_edge_stats.csv')
        
    return (pathway_overlap_edge_stats, pathway_overlap_edge_counts_rand)

# "hash" function to assign unique integers to each motif (no overlapping values ie the number of hashed integers is the same as the 
# number of motifs)
@jit(nopython=True)
def custom_hash(a):
    a_hashed = (a[:, 0]*5+100) + a[:, 1]
    return a_hashed

# fast function to count multilink motifs of large graphs using the hash function above
def count_motifs(A_mux, all_motif_types_list):
    A_mux_reshaped = A_mux.reshape(int(len(A_mux.ravel())/2), 2)
    all_motifs_hashed = custom_hash(A_mux_reshaped)
    all_motifs_hashed_bincounts = np.bincount(all_motifs_hashed)
    
    unique_motifs_hashed = custom_hash(all_motif_types_list)
    motif_hash_dict = {k: '%s%s'%(i, j) for (i, j), k in zip(all_motif_types_list, unique_motifs_hashed)}
    
    motif_counts_hashed = {motif_hash_dict[i]: j for i, j in zip(np.where(all_motifs_hashed_bincounts != 0)[0], 
                                       all_motifs_hashed_bincounts[np.where(all_motifs_hashed_bincounts != 0)]) 
                           if i in unique_motifs_hashed} # the if statement ensures that there are no negative edge types in the KEGG+PPI layer
    
    return motif_counts_hashed

# the 1D-array version of the function above for use with the shortest paths approach
def count_motifs_1D(A_mux, all_motif_types_list):
    A_mux_reshaped = A_mux[0]
    all_motifs_hashed = custom_hash(A_mux_reshaped)
    all_motifs_hashed_bincounts = np.bincount(all_motifs_hashed)
    
    unique_motifs_hashed = custom_hash(all_motif_types_list)
    motif_hash_dict = {k: '%s%s'%(i, j) for (i, j), k in zip(all_motif_types_list, unique_motifs_hashed)}
    
    motif_counts_hashed = {motif_hash_dict[i]: j for i, j in zip(np.where(all_motifs_hashed_bincounts != 0)[0], 
                                       all_motifs_hashed_bincounts[np.where(all_motifs_hashed_bincounts != 0)]) 
                           if i in unique_motifs_hashed} # the if statement ensures that there are no negative edge types in the KEGG+PPI layer
    
    return motif_counts_hashed

# multilink counts of the entire KEGG+PPI and GRN multiplex
def overall_multilink_counts(A_KEGGPPI_sparr, A_GRN_sparr, A_KEGG_e_sparr_dict, all_motif_types, all_motif_types_list, 
                             KEGG_interaction_types_dict, proj_path, input_GRN):

    if not os.path.exists(proj_path + input_GRN + '_overall_multilink_counts.pickle'):
        
        print('Calculating multilink counts in the entire multilayer network...')
        overall_multilink_counts_df_dict = {}
        for p in ['overall_multiplex']: 

            multilink_counts = {}   

            # only for [0X]
            A_mux = np.dstack([A_KEGGPPI_sparr.toarray(), A_GRN_sparr.toarray()])   
            e_multilink_counts = count_motifs(A_mux, all_motif_types_list)

            for motif in e_multilink_counts.keys():
                if motif in ['0-1', '00', '01']:
                    multilink_counts[motif] = e_multilink_counts[motif]    

            # for the remaining motifs (count ppi layer last since it's included in all KEGG layers)
            for e in tqdm(np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]), position=0, leave=True):       

                A_mux = np.dstack([A_KEGG_e_sparr_dict[e].toarray(), A_GRN_sparr.toarray()])
                e_multilink_counts = count_motifs(A_mux, all_motif_types_list)

                # disregard the motifs with no links in the signaling network (0X type motifs) since these are counted separately above
                # (since we don't want to multiple-count these for every interaction type e)
                for motif in e_multilink_counts.keys():
                    if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                        multilink_counts[motif] = e_multilink_counts[motif]

            multilink_counts_df = pd.DataFrame(0, index=all_motif_types, columns=['motif_counts'])
            for key in multilink_counts.keys():
                if key in all_motif_types:
                    multilink_counts_df.at[key, 'motif_counts'] = multilink_counts[key]


            overall_multilink_counts_df_dict[p] = multilink_counts_df  

        with open(proj_path + input_GRN + '_overall_multilink_counts.pickle', 'wb') as fp:
            pickle.dump(overall_multilink_counts_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)

    else:
        
        print('Reading multilink counts in the entire multilayer network...')
        with open(proj_path + input_GRN + '_overall_multilink_counts.pickle', 'rb') as fp:
            overall_multilink_counts_df_dict = pickle.load(fp)        
        
    return overall_multilink_counts_df_dict

# multilink counts of the randomized versions of the entire KEGG+PPI and GRN multiplex
# "rand_start" is the index of the random pair of layers to be started from. This is to both allow for a deterministic randomization that can be
# traced back and also to avoid different methods sharing the same randomized instances.
def overall_multilink_counts_rand(A_GRN_sparr_rand_dict, all_motif_types, A_KEGGPPI_sparr_rand_dict, all_motif_types_list, 
                                  KEGG_interaction_types_dict, A_KEGG_e_sparr_rand_dict, proj_path, input_GRN,
                                  rand_start = 0, N_rand = 100, N_swap = 10):

    if not os.path.exists(proj_path + input_GRN + '_overall_multilink_counts_rand_%srands.pickle' % N_rand):
    
    
        print('Calculating multilink counts in the entire multilayer network for randomized multilayer network ensebles...')
        overall_multilink_counts_rand_df_dict = {}
        all_pairs = [(i, j) for i, j in product(list(np.arange(len(A_GRN_sparr_rand_dict))), list(np.arange(len(A_GRN_sparr_rand_dict))))]

        for p in ['overall_multiplex']: 

            multilink_counts_rand_df = pd.DataFrame(0, index=all_motif_types, columns=np.arange(N_rand))

            for nrand in tqdm(np.arange(N_rand), position=0, leave=True):    

                # for easy reproducibility, deterministically iterate over the randomized ensembles in order
                rand_ix = all_pairs[nrand]

                multilink_counts_rand = {}   

                # slice GRN to get the pathway nodes only
                A_2_arr_rand = A_GRN_sparr_rand_dict[rand_ix[1]].toarray()
                A_1_arr_rand = A_KEGGPPI_sparr_rand_dict[rand_ix[0]].toarray()

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs [only for 0X]
                A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                e_multilink_counts_rand = count_motifs(A_mux_rand, all_motif_types_list)

                for motif in e_multilink_counts_rand.keys():
                    if motif in ['0-1', '00', '01']:
                        multilink_counts_rand[motif] = e_multilink_counts_rand[motif]    


                for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):       

                    # slice KEGG(e)+PPI to get the pathway nodes only
                    A_1_arr_rand = A_KEGG_e_sparr_rand_dict[rand_ix[0]][e].toarray()

                    # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                    A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                    e_multilink_counts_rand = count_motifs(A_mux_rand, all_motif_types_list)

                    # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                    # (since we don't want to multiple-count these for every interaction type e)
                    for motif in e_multilink_counts_rand.keys():
                        if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                            multilink_counts_rand[motif] = e_multilink_counts_rand[motif]


                for key in multilink_counts_rand.keys():
                    if key in all_motif_types:
                        multilink_counts_rand_df.at[key, nrand] = multilink_counts_rand[key]

            overall_multilink_counts_rand_df_dict[p] = multilink_counts_rand_df     

        with open(proj_path + input_GRN  + '_overall_multilink_counts_rand_%srands.pickle' % N_rand, 'wb') as fp:
            pickle.dump(overall_multilink_counts_rand_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)   
            
    else:
        
        print('Reading multilink counts in the entire multilayer network for randomized multilayer network ensebles...')
        with open(proj_path + input_GRN + '_overall_multilink_counts_rand_%srands.pickle' % N_rand, 'rb') as fp:
            overall_multilink_counts_rand_df_dict = pickle.load(fp)                
            
    return overall_multilink_counts_rand_df_dict

# multilink z-scores of the entire KEGG+PPI and GRN multiplex
def overall_multilink_zscores_pvals(KEGG_all_paths, overall_multilink_counts_rand_df_dict, overall_multilink_counts_df_dict, all_motif_types, 
                                    proj_path, input_GRN, N_rand=100):

    if not (os.path.exists(proj_path + input_GRN + '_overall_multilink_counts_stats.pickle') & 
            os.path.exists(proj_path + input_GRN + '_overall_rand_zscores.pickle')):
       
        print('Calculating multilink statistics for the entire multilayer network...')
        overall_multilink_counts_stats_df_dict = {}
        overall_rand_zscores_df_dict = {}

        for p in ['overall_multiplex']: 

            overall_multilink_counts_stats_df = pd.DataFrame(index=all_motif_types)
            overall_multilink_counts_stats_df['rand_mean'] = overall_multilink_counts_rand_df_dict[p].mean(axis=1)
            overall_multilink_counts_stats_df['rand_std'] = overall_multilink_counts_rand_df_dict[p].std(axis=1)
            overall_multilink_counts_stats_df['actual_motif_counts'] = overall_multilink_counts_df_dict[p]['motif_counts']
            nostd_ix = overall_multilink_counts_stats_df[overall_multilink_counts_stats_df['rand_std'] != 0].index
            overall_multilink_counts_stats_df['z-score'] = ((overall_multilink_counts_stats_df.loc[nostd_ix]['actual_motif_counts'] - 
                                                             overall_multilink_counts_stats_df.loc[nostd_ix]['rand_mean']) / 
                                                            overall_multilink_counts_stats_df.loc[nostd_ix]['rand_std'])

            pvals_df = pd.DataFrame(index=all_motif_types, columns=['p-value'])
            for m in overall_multilink_counts_stats_df.index:
                if overall_multilink_counts_stats_df['z-score'].loc[m] > 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(overall_multilink_counts_rand_df_dict[p].loc[m] >= 
                                                     overall_multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif overall_multilink_counts_stats_df['z-score'].loc[m] < 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(overall_multilink_counts_rand_df_dict[p].loc[m] <= 
                                                     overall_multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif pd.isnull(overall_multilink_counts_stats_df['z-score'].loc[m]) & (overall_multilink_counts_stats_df['actual_motif_counts'].loc[m] >
                                                                                       overall_multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0
                elif pd.isnull(overall_multilink_counts_stats_df['z-score'].loc[m]) & (overall_multilink_counts_stats_df['actual_motif_counts'].loc[m] <
                                                                                       overall_multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0

            overall_multilink_counts_stats_df['p-value'] = pvals_df['p-value'] 

            overall_rand_zscores_df = overall_multilink_counts_rand_df_dict[p].subtract(overall_multilink_counts_rand_df_dict[p].mean(axis=1), 
                                                                                        axis=0)\
                                                            .div(overall_multilink_counts_rand_df_dict[p].std(axis=1), axis=0)

            overall_multilink_counts_stats_df_dict[p] = overall_multilink_counts_stats_df
            overall_rand_zscores_df_dict[p] = overall_rand_zscores_df  

        with open(proj_path + input_GRN + '_overall_multilink_counts_stats.pickle', 'wb') as fp:
            pickle.dump(overall_multilink_counts_stats_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)  
        with open(proj_path + input_GRN + '_overall_rand_zscores.pickle', 'wb') as fp:
            pickle.dump(overall_rand_zscores_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        print('Reading multilink statistics for the entire multilayer network...')
        with open(proj_path + input_GRN + '_overall_multilink_counts_stats.pickle', 'rb') as fp:
            overall_multilink_counts_stats_df_dict = pickle.load(fp)
        with open(proj_path + input_GRN + '_overall_rand_zscores.pickle', 'rb') as fp:
            overall_rand_zscores_df_dict = pickle.load(fp)   
                        
    return (overall_multilink_counts_stats_df_dict, overall_rand_zscores_df_dict)

def plot_zscores_overall(overall_rand_zscores_df_dict, overall_multilink_counts_stats_df_dict, all_motif_types, proj_path, input_GRN, 
                         save_fig=False):

    # suppress the '00' types of multilink that throw off the ylim
    overall_rand_zscores_df_dict['overall_multiplex'].loc['00'] = np.nan
    overall_multilink_counts_stats_df_dict['overall_multiplex'].loc['00'] = np.nan         
    
    fig = plt.figure()
    fig.set_size_inches(10, 5)
    sns.set(font_scale=1.5)
    with sns.axes_style('ticks'):

        for p in ['overall_multiplex']:
            
            N_rand = overall_rand_zscores_df_dict[p].shape[1]
            plt.subplot(111)
            for i in np.arange(N_rand):
                plt.scatter(x=overall_rand_zscores_df_dict[p].index, y=overall_rand_zscores_df_dict[p][i], marker='+', s=150, c='grey', alpha=0.5)

            shapes = []
            colors = []
            for i in overall_multilink_counts_stats_df_dict[p].index:
                if ((overall_multilink_counts_stats_df_dict[p].loc[i]['z-score'] > 0) & 
                    (overall_multilink_counts_stats_df_dict[p].loc[i]['p-value'] <= 0.05)):
                    shapes.append('^')
                    colors.append('tomato')
                elif ((overall_multilink_counts_stats_df_dict[p].loc[i]['z-score'] < 0) & 
                      (overall_multilink_counts_stats_df_dict[p].loc[i]['p-value'] <= 0.05)):
                    shapes.append('v')
                    colors.append('dodgerblue')
                elif (overall_multilink_counts_stats_df_dict[p].loc[i]['p-value'] > 0.05):
                    shapes.append('o')
                    colors.append('k')
                else:
                    shapes.append('')
                    colors.append('w')     

            for xp, yp, m, c in zip(overall_rand_zscores_df_dict[p].index, overall_multilink_counts_stats_df_dict[p]['z-score'], shapes, colors):
                plt.scatter(x=xp, y=yp, marker=m, s=150, c=c)

            plt.plot([-1, len(all_motif_types)+1], [0, 0], '-', color='k')
            plt.xlim(-1, len(all_motif_types))
            plt.xlabel('Multilink motif')
            plt.ylabel('z-score')
            plt.xticks(rotation=90)
            plt.title('%s' % p)
            ax = plt.gca()
            ax.grid(which='major', axis='x', linestyle='--')
            plt.tight_layout()
        
        if save_fig == True:
            plt.savefig(proj_path + input_GRN + '_overall_multilink_zscores.pdf', format='pdf')
        
        plt.show()

# multilink counts of the individual KEGG signaling pathways
def multilink_counts(KEGG_all_paths, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, A_GRN_sparr, A_KEGGPPI_sparr, 
                     all_motif_types_list, all_motif_types, KEGG_interaction_types_dict, A_KEGG_e_sparr_dict, proj_path, input_GRN):

    if not os.path.exists(proj_path + input_GRN + '_multilink_counts.pickle'):
        
        print('Calculating multilink counts in each KEGG signaling pathway...')
        multilink_counts_df_dict = {}
        for p in tqdm(KEGG_all_paths, position=0, leave=True): 

            multilink_counts = {}   

            path_ix = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values

            # slice GRN to get the pathway nodes only
            A_2_arr = A_GRN_sparr[tuple([path_ix])][:, path_ix].toarray()
            A_1_arr = A_KEGGPPI_sparr[tuple([path_ix])][:, path_ix].toarray()

            # only for [0X]
            A_mux = np.dstack([A_1_arr, A_2_arr])    
            e_multilink_counts = count_motifs(A_mux, all_motif_types_list)

            for motif in e_multilink_counts.keys():
                if motif in ['0-1', '00', '01']:
                    multilink_counts[motif] = e_multilink_counts[motif]    

            # for the remaining motifs (count ppi layer last since it's included in all KEGG layers)
            for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):    

                # slice KEGG(e)+PPI to get the pathway nodes only
                A_1_arr = A_KEGG_e_sparr_dict[e][tuple([path_ix])][:, path_ix].toarray()

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                A_mux = np.dstack([A_1_arr, A_2_arr])

                e_multilink_counts = count_motifs(A_mux, all_motif_types_list)

                # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                # (since we don't want to multiple-count these for every interaction type e)
                for motif in e_multilink_counts.keys():
                    if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                        multilink_counts[motif] = e_multilink_counts[motif]


            multilink_counts_df = pd.DataFrame(0, index=all_motif_types, columns=['motif_counts'])
            for key in multilink_counts.keys():
                if key in all_motif_types:
                    multilink_counts_df.at[key, 'motif_counts'] = multilink_counts[key]


            multilink_counts_df_dict[p] = multilink_counts_df  

        with open(proj_path + input_GRN + '_multilink_counts.pickle', 'wb') as fp:
            pickle.dump(multilink_counts_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        print('Reading multilink counts in each KEGG signaling pathway...')
        with open(proj_path + input_GRN + '_multilink_counts.pickle', 'rb') as fp:
            multilink_counts_df_dict = pickle.load(fp)        
        
            
    return multilink_counts_df_dict

# multilink counts of randomized versions of individual KEGG signaling pathways
def multilink_counts_rand(A_GRN_sparr_rand_dict,KEGG_all_paths, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, 
                          all_motif_types, A_KEGGPPI_sparr_rand_dict, all_motif_types_list, KEGG_interaction_types_dict, 
                          A_KEGG_e_sparr_rand_dict, proj_path, input_GRN, rand_start = 0, N_rand = 100, N_swap = 10):

    if not os.path.exists(proj_path + input_GRN + '_multilink_counts_rand_%srands.pickle' % N_rand):    
    
        print('Calculating multilink counts in each KEGG signaling pathway for randomized multilayer network ensebles...')
        multilink_counts_rand_df_dict = {}
        all_pairs = [(i, j) for i, j in product(list(np.arange(len(A_GRN_sparr_rand_dict))), list(np.arange(len(A_GRN_sparr_rand_dict))))]

        for nn, p in enumerate(tqdm(KEGG_all_paths, position=0, leave=True)): 

            path_ix = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values

            multilink_counts_rand_df = pd.DataFrame(0, index=all_motif_types, columns=np.arange(N_rand))

            for nrand in np.arange(N_rand):    

                # for easy reproducibility, deterministically iterate over the randomized ensembles in order
                rand_ix = all_pairs[rand_start + (nn*N_rand+nrand)]

                multilink_counts_rand = {}   

                # slice GRN to get the pathway nodes only
                A_2_arr_rand = A_GRN_sparr_rand_dict[rand_ix[1]][tuple([path_ix])][:, path_ix].toarray()
                A_1_arr_rand = A_KEGGPPI_sparr_rand_dict[rand_ix[0]][tuple([path_ix])][:, path_ix].toarray()

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs [only for 0X]
                A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                e_multilink_counts_rand = count_motifs(A_mux_rand, all_motif_types_list)

                for motif in e_multilink_counts_rand.keys():
                    if motif in ['0-1', '00', '01']:
                        multilink_counts_rand[motif] = e_multilink_counts_rand[motif]    


                for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):       

                    # slice KEGG(e)+PPI to get the pathway nodes only
                    A_1_arr_rand = A_KEGG_e_sparr_rand_dict[rand_ix[0]][e][tuple([path_ix])][:, path_ix].toarray()

                    # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                    A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                    e_multilink_counts_rand = count_motifs(A_mux_rand, all_motif_types_list)

                    # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                    # (since we don't want to multiple-count these for every interaction type e)
                    for motif in e_multilink_counts_rand.keys():
                        if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                            multilink_counts_rand[motif] = e_multilink_counts_rand[motif]


                for key in multilink_counts_rand.keys():
                    if key in all_motif_types:
                        multilink_counts_rand_df.at[key, nrand] = multilink_counts_rand[key]

            multilink_counts_rand_df_dict[p] = multilink_counts_rand_df     

        with open(proj_path + input_GRN + '_multilink_counts_rand_%srands.pickle' % N_rand, 'wb') as fp:
            pickle.dump(multilink_counts_rand_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL) 
            
    else:
        
        print('Reading multilink counts in each KEGG signaling pathway for randomized multilayer network ensebles...')
        with open(proj_path + input_GRN + '_multilink_counts_rand_%srands.pickle' % N_rand, 'rb') as fp:
            multilink_counts_rand_df_dict = pickle.load(fp)               
            
    return multilink_counts_rand_df_dict

# multilink z-scores of individual KEGG signaling pathways
def multilink_zscores_pvals(KEGG_all_paths, multilink_counts_rand_df_dict, multilink_counts_df_dict, all_motif_types, proj_path, input_GRN,
                            N_rand=100):

    if not (os.path.exists(proj_path + input_GRN + '_multilink_counts_stats.pickle') &
            os.path.exists(proj_path + input_GRN + '_rand_zscores.pickle')):        
    
        print('Calculating multilink statistics for each KEGG signaling pathway...')
        multilink_counts_stats_df_dict = {}
        rand_zscores_df_dict = {}

        for p in tqdm(KEGG_all_paths, position=0, leave=True): 

            multilink_counts_stats_df = pd.DataFrame(index=all_motif_types)
            multilink_counts_stats_df['rand_mean'] = multilink_counts_rand_df_dict[p].mean(axis=1)
            multilink_counts_stats_df['rand_std'] = multilink_counts_rand_df_dict[p].std(axis=1)
            multilink_counts_stats_df['actual_motif_counts'] = multilink_counts_df_dict[p]['motif_counts']
            nostd_ix = multilink_counts_stats_df[multilink_counts_stats_df['rand_std'] != 0].index
            multilink_counts_stats_df['z-score'] = ((multilink_counts_stats_df.loc[nostd_ix]['actual_motif_counts'] - 
                                                                  multilink_counts_stats_df.loc[nostd_ix]['rand_mean']) / 
                                                                  multilink_counts_stats_df.loc[nostd_ix]['rand_std'])

            pvals_df = pd.DataFrame(index=all_motif_types, columns=['p-value'])
            for m in multilink_counts_stats_df.index:
                if multilink_counts_stats_df['z-score'].loc[m] > 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(multilink_counts_rand_df_dict[p].loc[m] >= 
                                                     multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif multilink_counts_stats_df['z-score'].loc[m] < 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(multilink_counts_rand_df_dict[p].loc[m] <= 
                                                     multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif pd.isnull(multilink_counts_stats_df['z-score'].loc[m]) & (multilink_counts_stats_df['actual_motif_counts'].loc[m] >
                                                                               multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0
                elif pd.isnull(multilink_counts_stats_df['z-score'].loc[m]) & (multilink_counts_stats_df['actual_motif_counts'].loc[m] <
                                                                               multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0

            multilink_counts_stats_df['p-value'] = pvals_df['p-value'] 

            rand_zscores_df = multilink_counts_rand_df_dict[p].subtract(multilink_counts_rand_df_dict[p].mean(axis=1), axis=0)\
                                                            .div(multilink_counts_rand_df_dict[p].std(axis=1), axis=0)

            multilink_counts_stats_df_dict[p] = multilink_counts_stats_df
            rand_zscores_df_dict[p] = rand_zscores_df  

        with open(proj_path + input_GRN + '_multilink_counts_stats.pickle', 'wb') as fp:
            pickle.dump(multilink_counts_stats_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)  
        with open(proj_path + input_GRN + '_rand_zscores.pickle', 'wb') as fp:
            pickle.dump(rand_zscores_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        print('Reading multilink statistics for each KEGG signaling pathway...')        
        with open(proj_path + input_GRN + '_multilink_counts_stats.pickle', 'rb') as fp:
            multilink_counts_stats_df_dict = pickle.load(fp)
        with open(proj_path + input_GRN + '_rand_zscores.pickle', 'rb') as fp:
            rand_zscores_df_dict = pickle.load(fp)      
                       
    return (multilink_counts_stats_df_dict, rand_zscores_df_dict)

def plot_zscores(KEGG_all_paths, rand_zscores_df_dict, multilink_counts_stats_df_dict, all_motif_types, proj_path, input_GRN, save_fig=False):
        
    fig = plt.figure()
    fig.set_size_inches(48, 60)
    sns.set(font_scale=1.5)
    with sns.axes_style('ticks'):

        for nfig, p in enumerate(tqdm(KEGG_all_paths, position=0, leave=True)):
            
            # suppress the '00' types of multilink that throw off the ylim
            rand_zscores_df_dict[p].loc['00'] = np.nan
            multilink_counts_stats_df_dict[p].loc['00'] = np.nan     
            
            N_rand = rand_zscores_df_dict[p].shape[1]
            plt.subplot(13, 5, nfig+1)
            for i in np.arange(N_rand):
                plt.scatter(x=rand_zscores_df_dict[p].index, y=rand_zscores_df_dict[p][i], marker='+', s=150, c='grey', alpha=0.5)

            shapes = []
            colors = []
            for i in multilink_counts_stats_df_dict[p].index:
                if (multilink_counts_stats_df_dict[p].loc[i]['z-score'] > 0) & (multilink_counts_stats_df_dict[p].loc[i]['p-value'] <= 0.05):
                    shapes.append('^')
                    colors.append('tomato')
                elif (multilink_counts_stats_df_dict[p].loc[i]['z-score'] < 0) & (multilink_counts_stats_df_dict[p].loc[i]['p-value'] <= 0.05):
                    shapes.append('v')
                    colors.append('dodgerblue')
                elif (multilink_counts_stats_df_dict[p].loc[i]['p-value'] > 0.05):
                    shapes.append('o')
                    colors.append('k')
                else:
                    shapes.append('')
                    colors.append('w')     

            for xp, yp, m, c in zip(rand_zscores_df_dict[p].index, multilink_counts_stats_df_dict[p]['z-score'], shapes, colors):
                plt.scatter(x=xp, y=yp, marker=m, s=150, c=c)

            plt.plot([-1, len(all_motif_types)+1], [0, 0], '-', color='k')
            plt.xlim(-1, len(all_motif_types))
            plt.xlabel('Multilink motif')
            plt.ylabel('z-score')
            plt.xticks(rotation=90)
            plt.title('%s' % p)
            ax = plt.gca()
            ax.grid(which='major', axis='x', linestyle='--')
            plt.tight_layout()
        if save_fig == True:
            plt.savefig(proj_path + input_GRN + '_multilink_zscores.pdf', format='pdf')
        
        plt.show()

### multilink counts of the edges between pairs of KEGG signaling pathways
# note that the below is for the smaller subset of pathways that are in XTalkDB for cross-validation purposes
def between_paths_multilink_counts(common_crosstalk_paths, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, 
                                  A_GRN_sparr, A_KEGGPPI_sparr, all_motif_types_list, all_motif_types, KEGG_interaction_types_dict, 
                                  A_KEGG_e_sparr_dict, proj_path, input_GRN):

    if not os.path.exists(proj_path + input_GRN + '_between_paths_multilink_counts.pickle'):        
    
        print('Calculating multilink counts of direct edges between all pairs of KEGG signaling pathways...')
        between_paths_multilink_counts_df_dict = {}

        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]

        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            multilink_counts = {}   

            path_ix1 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p1]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values
            path_ix2 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p2]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values

            # slice GRN to get the pathway nodes only
            A_2_arr = A_GRN_sparr[tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()
            A_1_arr = A_KEGGPPI_sparr[tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()

            # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs [only for 0X]
            A_mux = np.dstack([A_1_arr, A_2_arr])
            e_multilink_counts = count_motifs(A_mux, all_motif_types_list)

            for motif in e_multilink_counts.keys():
                if motif in ['0-1', '00', '01']:
                    multilink_counts[motif] = e_multilink_counts[motif]    


            for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):       

                # slice KEGG(e)+PPI to get the pathway nodes only
                A_1_arr = A_KEGG_e_sparr_dict[e][tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                A_mux = np.dstack([A_1_arr, A_2_arr])
                e_multilink_counts = count_motifs(A_mux, all_motif_types_list)

                # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                # (since we don't want to multiple-count these for every interaction type e)
                for motif in e_multilink_counts.keys():
                    if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                        multilink_counts[motif] = e_multilink_counts[motif]


            multilink_counts_df = pd.DataFrame(0, index=all_motif_types, columns=['motif_counts'])
            for key in multilink_counts.keys():
                if key in all_motif_types:
                    multilink_counts_df.at[key, 'motif_counts'] = multilink_counts[key]


            between_paths_multilink_counts_df_dict[(p1, p2)] = multilink_counts_df  

        with open(proj_path + input_GRN + '_between_paths_multilink_counts.pickle', 'wb') as fp:
            pickle.dump(between_paths_multilink_counts_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        print('Reading multilink counts of direct edges between all pairs of KEGG signaling pathways...')
        with open(proj_path + input_GRN + '_between_paths_multilink_counts.pickle', 'rb') as fp:
            between_paths_multilink_counts_df_dict = pickle.load(fp)         
            
    return between_paths_multilink_counts_df_dict

def between_paths_multilink_counts_rand(A_GRN_sparr_rand_dict, common_crosstalk_paths, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, 
                                        all_motif_types,A_KEGGPPI_sparr_rand_dict, all_motif_types_list, KEGG_interaction_types_dict, 
                                        A_KEGG_e_sparr_rand_dict, proj_path, input_GRN, rand_start = 0, N_rand = 100, N_swap = 10):

    if not os.path.exists(proj_path + input_GRN + '_between_paths_multilink_counts_rand_%srands.pickle' % N_rand):     
        
        print('Calculating multilink counts of direct edges between all pairs of KEGG signaling pathways for randomized multilayer network ensebles...')
        between_paths_multilink_counts_rand_df_dict = {}
        all_pairs = [(i, j) for i, j in product(list(np.arange(len(A_GRN_sparr_rand_dict))), list(np.arange(len(A_GRN_sparr_rand_dict))))]
        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]

        for nn, (p1, p2) in enumerate(tqdm(all_perms, position=0, leave=True)): 

            path_ix1 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p1]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values
            path_ix2 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p2]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values

            multilink_counts_rand_df = pd.DataFrame(0, index=all_motif_types, columns=np.arange(N_rand))

            for nrand in np.arange(N_rand):    

                rand_ix = all_pairs[rand_start + (nn*N_rand+nrand)]

                multilink_counts_rand = {}   

                # slice GRN to get the pathway nodes only
                A_2_arr_rand = A_GRN_sparr_rand_dict[rand_ix[1]][tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()   
                A_1_arr_rand = A_KEGGPPI_sparr_rand_dict[rand_ix[0]][tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs [only for 0X]
                A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                e_multilink_counts_rand = count_motifs(A_mux_rand, all_motif_types_list)

                for motif in e_multilink_counts_rand.keys():
                    if motif in ['0-1', '00', '01']:
                        multilink_counts_rand[motif] = e_multilink_counts_rand[motif]    


                for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):       

                    # slice KEGG(e)+PPI to get the pathway nodes only
                    A_1_arr_rand = A_KEGG_e_sparr_rand_dict[rand_ix[0]][e][tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()

                    # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                    A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                    e_multilink_counts_rand = count_motifs(A_mux_rand, all_motif_types_list)

                    # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                    # (since we don't want to multiple-count these for every interaction type e)
                    for motif in e_multilink_counts_rand.keys():
                        if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                            multilink_counts_rand[motif] = e_multilink_counts_rand[motif]


                for key in multilink_counts_rand.keys():
                    if key in all_motif_types:
                        multilink_counts_rand_df.at[key, nrand] = multilink_counts_rand[key]

            between_paths_multilink_counts_rand_df_dict[(p1, p2)] = multilink_counts_rand_df     

        with open(proj_path + input_GRN + '_between_paths_multilink_counts_rand_%srands.pickle' % N_rand, 'wb') as fp:
            pickle.dump(between_paths_multilink_counts_rand_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)  
            
    else:
        
        print('Reading multilink counts of direct edges between all pairs of KEGG signaling pathways for randomized multilayer network ensebles...')
        with open(proj_path + input_GRN + '_between_paths_multilink_counts_rand_%srands.pickle' % N_rand, 'rb') as fp:
            between_paths_multilink_counts_rand_df_dict = pickle.load(fp)           
            
    return between_paths_multilink_counts_rand_df_dict

def between_paths_multilink_zscores_pvals(common_crosstalk_paths, between_paths_multilink_counts_rand_df_dict, 
                                          between_paths_multilink_counts_df_dict, all_motif_types, proj_path, input_GRN, N_rand=100):

    if not (os.path.exists(proj_path + input_GRN + '_between_paths_multilink_counts_stats.pickle') &
            os.path.exists(proj_path + input_GRN + '_between_paths_rand_zscores.pickle')):       
       
        print('Calculating multilink statistics for direct edges between all pairs of KEGG signaling pathways...')
        between_paths_multilink_counts_stats_df_dict = {}
        between_paths_rand_zscores_df_dict = {}

        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]

        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            between_paths_multilink_counts_stats_df = pd.DataFrame(index=all_motif_types)
            between_paths_multilink_counts_stats_df['rand_mean'] = between_paths_multilink_counts_rand_df_dict[(p1, p2)].mean(axis=1)
            between_paths_multilink_counts_stats_df['rand_std'] = between_paths_multilink_counts_rand_df_dict[(p1, p2)].std(axis=1)
            between_paths_multilink_counts_stats_df['actual_motif_counts'] = between_paths_multilink_counts_df_dict[(p1, p2)]['motif_counts']
            nostd_ix = between_paths_multilink_counts_stats_df[between_paths_multilink_counts_stats_df['rand_std'] != 0].index
            between_paths_multilink_counts_stats_df['z-score'] = ((between_paths_multilink_counts_stats_df.loc[nostd_ix]['actual_motif_counts'] - 
                                                                  between_paths_multilink_counts_stats_df.loc[nostd_ix]['rand_mean']) / 
                                                                  between_paths_multilink_counts_stats_df.loc[nostd_ix]['rand_std'])

            pvals_df = pd.DataFrame(index=all_motif_types, columns=['p-value'])
            for m in between_paths_multilink_counts_stats_df.index:
                if between_paths_multilink_counts_stats_df['z-score'].loc[m] > 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(between_paths_multilink_counts_rand_df_dict[(p1, p2)].loc[m] >= 
                                                     between_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif between_paths_multilink_counts_stats_df['z-score'].loc[m] < 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(between_paths_multilink_counts_rand_df_dict[(p1, p2)].loc[m] <= 
                                                     between_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif pd.isnull(between_paths_multilink_counts_stats_df['z-score'].loc[m]) & (between_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m] >
                                                                                            between_paths_multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0
                elif pd.isnull(between_paths_multilink_counts_stats_df['z-score'].loc[m]) & (between_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m] <
                                                                                            between_paths_multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0

            between_paths_multilink_counts_stats_df['p-value'] = pvals_df['p-value'] 

            between_paths_rand_zscores_df = between_paths_multilink_counts_rand_df_dict[(p1, p2)].subtract(between_paths_multilink_counts_rand_df_dict[(p1, p2)].mean(axis=1), axis=0).div(between_paths_multilink_counts_rand_df_dict[(p1, p2)].std(axis=1), axis=0)

            between_paths_multilink_counts_stats_df_dict[(p1, p2)] = between_paths_multilink_counts_stats_df
            between_paths_rand_zscores_df_dict[(p1, p2)] = between_paths_rand_zscores_df  


        with open(proj_path + input_GRN  + '_between_paths_multilink_counts_stats.pickle', 'wb') as fp:
            pickle.dump(between_paths_multilink_counts_stats_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)  
        with open(proj_path + input_GRN  + '_between_paths_rand_zscores.pickle', 'wb') as fp:
            pickle.dump(between_paths_rand_zscores_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)  
            
    else:
        
        print('Reading multilink statistics for direct edges between all pairs of KEGG signaling pathways...')
        with open(proj_path + input_GRN + '_between_paths_multilink_counts_stats.pickle', 'rb') as fp:
            between_paths_multilink_counts_stats_df_dict = pickle.load(fp)
        with open(proj_path + input_GRN + '_between_paths_rand_zscores.pickle', 'rb') as fp:
            between_paths_rand_zscores_df_dict = pickle.load(fp)            

    return (between_paths_multilink_counts_stats_df_dict, between_paths_rand_zscores_df_dict)

# gets the part of the shortest path between pathway1 and pathway2 that excludes edges that are completely in either pathway
def intermediary_sp(temp_sp, intermediary_nodes):

    inter_ind = np.where(np.in1d(np.array(temp_sp), intermediary_nodes))[0] # get the indices of these intermediary nodes in the shortest path array
    pre_ind = inter_ind[0] - 1 # get the index of the node preceding the first intermediary node
    post_ind = inter_ind[-1] + 1 # get the index of the node succeeding the last intermediary node
    inter_sp_ind = np.sort(np.hstack([pre_ind, inter_ind, post_ind])) # indices of the intermediary nodes and their connections to the pathways
    inter_sp = np.array(temp_sp)[inter_sp_ind] # the intermediary SP consisting of the intermediary nodes and the edges connecting them to the pathways

    return inter_sp

def get_shortest_paths(A_KEGGPPI_sparr, common_crosstalk_paths, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, proj_path, 
                       sp_threshold=None):

    if not (os.path.exists(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_edges_dict.pickle') &
            os.path.exists(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_intermediaries_dict.pickle')):     
    
        print('Calculating the shortest paths with threshold (%s) between all pairs of KEGG signaling pathways...' % str(sp_threshold))
        G_KEGGPPI = nx.from_scipy_sparse_matrix(A_KEGGPPI_sparr.astype(bool).astype(int))
        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]

        shortest_path_intermediaries_dict = {}
        shortest_path_edges_dict = {}
        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            multilink_counts = {}   

            path_ix1 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p1]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values
            path_ix2 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p2]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values

            shortest_path_edges_dict[(p1, p2)] = set()
            shortest_path_intermediaries_dict[(p1, p2)] = set()

            if sp_threshold == None:
                for pp1 in np.setdiff1d(path_ix1, path_ix2):
                    for pp2 in np.setdiff1d(path_ix2, path_ix1):
                        if nx.has_path(G_KEGGPPI, source=pp1, target=pp2):
                            temp_sp = nx.shortest_path(G_KEGGPPI, source=pp1, target=pp2)
                            intermediary_nodes = list(set(temp_sp) - (set(path_ix1) | set(path_ix2))) # intermediary nodes that are a part of neither pathway
                            if len(intermediary_nodes) > 0:
                                inter_sp = intermediary_sp(temp_sp, intermediary_nodes)
                                shortest_path_edges_dict[(p1, p2)].update(set([(i, j) for i, j in zip(inter_sp, inter_sp[1:])]))
                                shortest_path_intermediaries_dict[(p1, p2)].update(set(intermediary_nodes))
            else:
                for pp1 in np.setdiff1d(path_ix1, path_ix2):
                    for pp2 in np.setdiff1d(path_ix2, path_ix1):
                        if nx.has_path(G_KEGGPPI, source=pp1, target=pp2):
                            temp_sp = nx.shortest_path(G_KEGGPPI, source=pp1, target=pp2)
                            intermediary_nodes = list(set(temp_sp) - (set(path_ix1) | set(path_ix2))) # intermediary nodes that are a part of neither pathway
                            if (len(intermediary_nodes) > 0) & (len(intermediary_nodes) <= sp_threshold):
                                inter_sp = intermediary_sp(temp_sp, intermediary_nodes)
                                shortest_path_edges_dict[(p1, p2)].update(set([(i, j) for i, j in zip(inter_sp, inter_sp[1:])]))
                                shortest_path_intermediaries_dict[(p1, p2)].update(set(intermediary_nodes))   

        with open(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_edges_dict.pickle', 'wb') as fp:
            pickle.dump(shortest_path_edges_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
        with open(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_intermediaries_dict.pickle', 'wb') as fp:
            pickle.dump(shortest_path_intermediaries_dict, fp, protocol=pickle.HIGHEST_PROTOCOL) 
            
    else:
        
        print('Reading the shortest paths with threshold (%s) between all pairs of KEGG signaling pathways...' % str(sp_threshold))
        with open(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_edges_dict.pickle', 'rb') as fp:
            shortest_path_edges_dict = pickle.load(fp)
        with open(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_intermediaries_dict.pickle', 'rb') as fp:
            shortest_path_intermediaries_dict = pickle.load(fp)         
                    
    return (shortest_path_edges_dict, shortest_path_intermediaries_dict)

### multilink counts of the shortest path edges between pairs of KEGG signaling pathways
# note that the below is for the smaller subset of pathways that are in XTalkDB for cross-validation purposes
def shortest_paths_multilink_counts(common_crosstalk_paths, shortest_path_edges_dict, sp_threshold,
                                    A_GRN_sparr, A_KEGGPPI_sparr, all_motif_types, all_motif_types_list, KEGG_interaction_types_dict, 
                                    A_KEGG_e_sparr_dict, proj_path, input_GRN):

    if not os.path.exists(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts.pickle'):  
        
        print('Calculating multilink counts of shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways...' 
              % str(sp_threshold))
        shortest_paths_multilink_counts_df_dict = {}

        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]

        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            multilink_counts = {}   

            # fancy index the sparse arrays to get the shortest path edges between the pair of pathways
            A_2_arr = np.asarray(A_GRN_sparr[np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                             np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]
            A_1_arr = np.asarray(A_KEGGPPI_sparr[np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                                 np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]

            # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs [only for 0X]
            A_mux = np.dstack([A_1_arr, A_2_arr])
            e_multilink_counts = count_motifs_1D(A_mux, all_motif_types_list)

            for motif in e_multilink_counts.keys():
                if motif in ['0-1', '00', '01']:
                    multilink_counts[motif] = e_multilink_counts[motif]    


            for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):       

                # slice KEGG(e)+PPI to get the pathway nodes only
                A_1_arr = np.asarray(A_KEGG_e_sparr_dict[e][np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                                            np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                A_mux = np.dstack([A_1_arr, A_2_arr])
                e_multilink_counts = count_motifs_1D(A_mux, all_motif_types_list)

                # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                # (since we don't want to multiple-count these for every interaction type e)
                for motif in e_multilink_counts.keys():
                    if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                        multilink_counts[motif] = e_multilink_counts[motif]


            multilink_counts_df = pd.DataFrame(0, index=all_motif_types, columns=['motif_counts'])
            for key in multilink_counts.keys():
                if key in all_motif_types:
                    multilink_counts_df.at[key, 'motif_counts'] = multilink_counts[key]


            shortest_paths_multilink_counts_df_dict[(p1, p2)] = multilink_counts_df  

        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts.pickle', 'wb') as fp:
            pickle.dump(shortest_paths_multilink_counts_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        print('Reading multilink counts of shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways...' 
              % str(sp_threshold))
        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts.pickle', 'rb') as fp:
            shortest_paths_multilink_counts_df_dict = pickle.load(fp)        
            
    return shortest_paths_multilink_counts_df_dict

def shortest_paths_multilink_counts_rand(A_GRN_sparr_rand_dict, common_crosstalk_paths, shortest_path_edges_dict, sp_threshold,
                                        all_motif_types, A_KEGGPPI_sparr_rand_dict, all_motif_types_list, KEGG_interaction_types_dict, 
                                        A_KEGG_e_sparr_rand_dict, proj_path, input_GRN, rand_start = 0, N_rand = 100, N_swap = 10):

    if not os.path.exists(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_rand_%srands.pickle' % N_rand):     
    
        print('Calculating multilink counts of shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways for randomized multilayer network ensebles...' 
              % str(sp_threshold))
        shortest_paths_multilink_counts_rand_df_dict = {}
        all_pairs = [(i, j) for i, j in product(list(np.arange(len(A_GRN_sparr_rand_dict))), list(np.arange(len(A_GRN_sparr_rand_dict))))]
        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]

        for nn, (p1, p2) in enumerate(tqdm(all_perms, position=0, leave=True)): 

            multilink_counts_rand_df = pd.DataFrame(0, index=all_motif_types, columns=np.arange(N_rand))

            for nrand in np.arange(N_rand):    

                rand_ix = all_pairs[rand_start + (nn*N_rand+nrand)]

                multilink_counts_rand = {}   

                # slice GRN to get the pathway nodes only
                A_2_arr_rand = np.asarray(A_GRN_sparr_rand_dict[rand_ix[1]][np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                                                            np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]  
                A_1_arr_rand = np.asarray(A_KEGGPPI_sparr_rand_dict[rand_ix[0]][np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                                                                np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs [only for 0X]
                A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                e_multilink_counts_rand = count_motifs_1D(A_mux_rand, all_motif_types_list)

                for motif in e_multilink_counts_rand.keys():
                    if motif in ['0-1', '00', '01']:
                        multilink_counts_rand[motif] = e_multilink_counts_rand[motif]    


                for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):       

                    # slice KEGG(e)+PPI to get the pathway nodes only
                    A_1_arr_rand = np.asarray(A_KEGG_e_sparr_rand_dict[rand_ix[0]][e][np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0],
                                                                                      np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]

                    # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                    A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                    e_multilink_counts_rand = count_motifs_1D(A_mux_rand, all_motif_types_list)

                    # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                    # (since we don't want to multiple-count these for every interaction type e)
                    for motif in e_multilink_counts_rand.keys():
                        if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                            multilink_counts_rand[motif] = e_multilink_counts_rand[motif]


                for key in multilink_counts_rand.keys():
                    if key in all_motif_types:
                        multilink_counts_rand_df.at[key, nrand] = multilink_counts_rand[key]

            shortest_paths_multilink_counts_rand_df_dict[(p1, p2)] = multilink_counts_rand_df     

        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_rand_%srands.pickle' % N_rand, 'wb') as fp:
            pickle.dump(shortest_paths_multilink_counts_rand_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        print('Reading multilink counts of shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways for randomized multilayer network ensebles...' 
              % str(sp_threshold))        
        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_rand_%srands.pickle' % N_rand, 'rb') as fp:
            shortest_paths_multilink_counts_rand_df_dict = pickle.load(fp)        
            
    return shortest_paths_multilink_counts_rand_df_dict

def shortest_paths_multilink_zscores_pvals(common_crosstalk_paths, shortest_paths_multilink_counts_rand_df_dict, sp_threshold,
                                          shortest_paths_multilink_counts_df_dict, all_motif_types, proj_path, input_GRN, N_rand=100):

    if not (os.path.exists(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_stats.pickle') &
            os.path.exists(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_rand_zscores.pickle')):         
    
        print('Calculating multilink statistics for shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways...' 
              % str(sp_threshold))
        shortest_paths_multilink_counts_stats_df_dict = {}
        shortest_paths_rand_zscores_df_dict = {}

        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]

        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            shortest_paths_multilink_counts_stats_df = pd.DataFrame(index=all_motif_types)
            shortest_paths_multilink_counts_stats_df['rand_mean'] = shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].mean(axis=1)
            shortest_paths_multilink_counts_stats_df['rand_std'] = shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].std(axis=1)
            shortest_paths_multilink_counts_stats_df['actual_motif_counts'] = shortest_paths_multilink_counts_df_dict[(p1, p2)]['motif_counts']
            nostd_ix = shortest_paths_multilink_counts_stats_df[shortest_paths_multilink_counts_stats_df['rand_std'] != 0].index
            shortest_paths_multilink_counts_stats_df['z-score'] = ((shortest_paths_multilink_counts_stats_df.loc[nostd_ix]['actual_motif_counts'] - 
                                                                  shortest_paths_multilink_counts_stats_df.loc[nostd_ix]['rand_mean']) / 
                                                                  shortest_paths_multilink_counts_stats_df.loc[nostd_ix]['rand_std'])

            pvals_df = pd.DataFrame(index=all_motif_types, columns=['p-value'])
            for m in shortest_paths_multilink_counts_stats_df.index:
                if shortest_paths_multilink_counts_stats_df['z-score'].loc[m] > 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].loc[m] >= 
                                                     shortest_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif shortest_paths_multilink_counts_stats_df['z-score'].loc[m] < 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].loc[m] <= 
                                                     shortest_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif pd.isnull(shortest_paths_multilink_counts_stats_df['z-score'].loc[m]) & (shortest_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m] >
                                                                                            shortest_paths_multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0
                elif pd.isnull(shortest_paths_multilink_counts_stats_df['z-score'].loc[m]) & (shortest_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m] <
                                                                                            shortest_paths_multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0

            shortest_paths_multilink_counts_stats_df['p-value'] = pvals_df['p-value'] 

            shortest_paths_rand_zscores_df = shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].subtract(shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].mean(axis=1), axis=0).div(shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].std(axis=1), axis=0)

            shortest_paths_multilink_counts_stats_df_dict[(p1, p2)] = shortest_paths_multilink_counts_stats_df
            shortest_paths_rand_zscores_df_dict[(p1, p2)] = shortest_paths_rand_zscores_df  

        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_stats.pickle', 'wb') as fp:
            pickle.dump(shortest_paths_multilink_counts_stats_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)  
        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_rand_zscores.pickle', 'wb') as fp:
            pickle.dump(shortest_paths_rand_zscores_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:

        print('Reading multilink statistics for shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways...' 
              % str(sp_threshold))
        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_stats.pickle', 'rb') as fp:
            shortest_paths_multilink_counts_stats_df_dict =  pickle.load(fp)  
        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_rand_zscores.pickle', 'rb') as fp:
            shortest_paths_rand_zscores_df_dict = pickle.load(fp)
    
    return (shortest_paths_multilink_counts_stats_df_dict, shortest_paths_rand_zscores_df_dict)

### multilink counts of the edges between pairs of KEGG signaling pathways
# note that the below is for the smaller subset of pathways that are in XTalkDB for cross-validation purposes
# The output file has the suffix "discovery" to indicate that this is all the KEGG signaling pathways that we have (61 in total) 
# and is not limited to the 25 benchmark pathways.
def between_paths_multilink_counts_discovery(discovery_crosstalk_paths, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, 
                                  A_GRN_sparr, A_KEGGPPI_sparr, all_motif_types_list, all_motif_types, KEGG_interaction_types_dict, 
                                  A_KEGG_e_sparr_dict, proj_path, input_GRN):

    if not os.path.exists(proj_path + input_GRN + '_between_paths_multilink_counts_discovery.pickle'):    
        
        print('Calculating multilink counts of direct edges between all pairs of KEGG signaling pathways...')    
        between_paths_multilink_counts_df_dict = {}

        all_perms = [(i, j) for i, j in permutations(discovery_crosstalk_paths, 2)]

        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            multilink_counts = {}   

            path_ix1 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p1]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values
            path_ix2 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p2]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values

            # slice GRN to get the pathway nodes only
            A_2_arr = A_GRN_sparr[tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()
            A_1_arr = A_KEGGPPI_sparr[tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()

            # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs [only for 0X]
            A_mux = np.dstack([A_1_arr, A_2_arr])
            e_multilink_counts = count_motifs(A_mux, all_motif_types_list)

            for motif in e_multilink_counts.keys():
                if motif in ['0-1', '00', '01']:
                    multilink_counts[motif] = e_multilink_counts[motif]    


            for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):       

                # slice KEGG(e)+PPI to get the pathway nodes only
                A_1_arr = A_KEGG_e_sparr_dict[e][tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                A_mux = np.dstack([A_1_arr, A_2_arr])
                e_multilink_counts = count_motifs(A_mux, all_motif_types_list)

                # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                # (since we don't want to multiple-count these for every interaction type e)
                for motif in e_multilink_counts.keys():
                    if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                        multilink_counts[motif] = e_multilink_counts[motif]


            multilink_counts_df = pd.DataFrame(0, index=all_motif_types, columns=['motif_counts'])
            for key in multilink_counts.keys():
                if key in all_motif_types:
                    multilink_counts_df.at[key, 'motif_counts'] = multilink_counts[key]


            between_paths_multilink_counts_df_dict[(p1, p2)] = multilink_counts_df  

        with open(proj_path + input_GRN + '_between_paths_multilink_counts_discovery.pickle', 'wb') as fp:
            pickle.dump(between_paths_multilink_counts_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        print('Reading multilink counts of direct edges between all pairs of KEGG signaling pathways...')  
        with open(proj_path + input_GRN + '_between_paths_multilink_counts_discovery.pickle', 'rb') as fp:
            between_paths_multilink_counts_df_dict = pickle.load(fp)        
            
    return between_paths_multilink_counts_df_dict

# The version where randomized instances are selected uniformly at random, as opposed to in order for reproducibility. The reason we need is that
# we need more than 250k (500 x 500) unique layer combinations in the full pathways case (3660 x 100). The output file has the suffix "discovery"
# to indicate that this is all the KEGG signaling pathways that we have (61 in total) and is not limited to the 25 benchmark pathways.
def between_paths_multilink_counts_rand_discovery(A_GRN_sparr_rand_dict, discovery_crosstalk_paths, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, 
                                        all_motif_types,A_KEGGPPI_sparr_rand_dict, all_motif_types_list, KEGG_interaction_types_dict, 
                                        A_KEGG_e_sparr_rand_dict, proj_path, input_GRN, N_rand = 100, N_swap = 10):

    if not os.path.exists(proj_path + input_GRN + '_between_paths_multilink_counts_rand_%srands_discovery.pickle' % N_rand):     
        
        print('Calculating multilink counts of direct edges between all pairs of KEGG signaling pathways for randomized multilayer network ensebles...')    
    
        between_paths_multilink_counts_rand_df_dict = {}
        all_pairs = [(i, j) for i, j in product(list(np.arange(len(A_GRN_sparr_rand_dict))), list(np.arange(len(A_GRN_sparr_rand_dict))))]
        all_perms = [(i, j) for i, j in permutations(discovery_crosstalk_paths, 2)]
        all_pairs_shuff = [(i, j) for i, j in product(np.random.choice(len(A_GRN_sparr_rand_dict), 600), 
                                                      np.random.choice(len(A_GRN_sparr_rand_dict), 600))]

        for nn, (p1, p2) in enumerate(tqdm(all_perms, position=0, leave=True)): 

            path_ix1 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p1]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values
            path_ix2 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p2]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values

            multilink_counts_rand_df = pd.DataFrame(0, index=all_motif_types, columns=np.arange(N_rand))

            for nrand in np.arange(N_rand):    

                rand_ix = all_pairs_shuff[nn*N_rand+nrand]

                multilink_counts_rand = {}   

                # slice GRN to get the pathway nodes only
                A_2_arr_rand = A_GRN_sparr_rand_dict[rand_ix[1]][tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()   
                A_1_arr_rand = A_KEGGPPI_sparr_rand_dict[rand_ix[0]][tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs [only for 0X]
                A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                e_multilink_counts_rand = count_motifs(A_mux_rand, all_motif_types_list)

                for motif in e_multilink_counts_rand.keys():
                    if motif in ['0-1', '00', '01']:
                        multilink_counts_rand[motif] = e_multilink_counts_rand[motif]    


                for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):       

                    # slice KEGG(e)+PPI to get the pathway nodes only
                    A_1_arr_rand = A_KEGG_e_sparr_rand_dict[rand_ix[0]][e][tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()

                    # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                    A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                    e_multilink_counts_rand = count_motifs(A_mux_rand, all_motif_types_list)

                    # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                    # (since we don't want to multiple-count these for every interaction type e)
                    for motif in e_multilink_counts_rand.keys():
                        if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                            multilink_counts_rand[motif] = e_multilink_counts_rand[motif]


                for key in multilink_counts_rand.keys():
                    if key in all_motif_types:
                        multilink_counts_rand_df.at[key, nrand] = multilink_counts_rand[key]

            between_paths_multilink_counts_rand_df_dict[(p1, p2)] = multilink_counts_rand_df     

        with open(proj_path + input_GRN + '_between_paths_multilink_counts_rand_%srands_discovery.pickle' % N_rand, 'wb') as fp:
            pickle.dump(between_paths_multilink_counts_rand_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        print('Reading multilink counts of direct edges between all pairs of KEGG signaling pathways for randomized multilayer network ensebles...')  
        with open(proj_path + input_GRN + '_between_paths_multilink_counts_rand_%srands_discovery.pickle' % N_rand, 'rb') as fp:
            between_paths_multilink_counts_rand_df_dict = pickle.load(fp)         
            
    return between_paths_multilink_counts_rand_df_dict

def between_paths_multilink_zscores_pvals_discovery(discovery_crosstalk_paths, between_paths_multilink_counts_rand_df_dict, 
                                          between_paths_multilink_counts_df_dict, all_motif_types, proj_path, input_GRN, N_rand=100):

    if not (os.path.exists(proj_path + input_GRN + '_between_paths_multilink_counts_stats_discovery.pickle') &
            os.path.exists(proj_path + input_GRN + '_between_paths_rand_zscores_discovery.pickle')):       
       
        print('Calculating multilink statistics for direct edges between all pairs of KEGG signaling pathways...')
        between_paths_multilink_counts_stats_df_dict = {}
        between_paths_rand_zscores_df_dict = {}

        all_perms = [(i, j) for i, j in permutations(discovery_crosstalk_paths, 2)]

        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            between_paths_multilink_counts_stats_df = pd.DataFrame(index=all_motif_types)
            between_paths_multilink_counts_stats_df['rand_mean'] = between_paths_multilink_counts_rand_df_dict[(p1, p2)].mean(axis=1)
            between_paths_multilink_counts_stats_df['rand_std'] = between_paths_multilink_counts_rand_df_dict[(p1, p2)].std(axis=1)
            between_paths_multilink_counts_stats_df['actual_motif_counts'] = between_paths_multilink_counts_df_dict[(p1, p2)]['motif_counts']
            nostd_ix = between_paths_multilink_counts_stats_df[between_paths_multilink_counts_stats_df['rand_std'] != 0].index
            between_paths_multilink_counts_stats_df['z-score'] = ((between_paths_multilink_counts_stats_df.loc[nostd_ix]['actual_motif_counts'] - 
                                                                  between_paths_multilink_counts_stats_df.loc[nostd_ix]['rand_mean']) / 
                                                                  between_paths_multilink_counts_stats_df.loc[nostd_ix]['rand_std'])

            pvals_df = pd.DataFrame(index=all_motif_types, columns=['p-value'])
            for m in between_paths_multilink_counts_stats_df.index:
                if between_paths_multilink_counts_stats_df['z-score'].loc[m] > 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(between_paths_multilink_counts_rand_df_dict[(p1, p2)].loc[m] >= 
                                                     between_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif between_paths_multilink_counts_stats_df['z-score'].loc[m] < 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(between_paths_multilink_counts_rand_df_dict[(p1, p2)].loc[m] <= 
                                                     between_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif pd.isnull(between_paths_multilink_counts_stats_df['z-score'].loc[m]) & (between_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m] >
                                                                                            between_paths_multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0
                elif pd.isnull(between_paths_multilink_counts_stats_df['z-score'].loc[m]) & (between_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m] <
                                                                                            between_paths_multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0

            between_paths_multilink_counts_stats_df['p-value'] = pvals_df['p-value'] 

            between_paths_rand_zscores_df = between_paths_multilink_counts_rand_df_dict[(p1, p2)].subtract(between_paths_multilink_counts_rand_df_dict[(p1, p2)].mean(axis=1), axis=0).div(between_paths_multilink_counts_rand_df_dict[(p1, p2)].std(axis=1), axis=0)

            between_paths_multilink_counts_stats_df_dict[(p1, p2)] = between_paths_multilink_counts_stats_df
            between_paths_rand_zscores_df_dict[(p1, p2)] = between_paths_rand_zscores_df  

        with open(proj_path + input_GRN + '_between_paths_multilink_counts_stats_discovery.pickle', 'wb') as fp:
            pickle.dump(between_paths_multilink_counts_stats_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)  
        with open(proj_path + input_GRN + '_between_paths_rand_zscores_discovery.pickle', 'wb') as fp:
            pickle.dump(between_paths_rand_zscores_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)  
            
    else:
        
        print('Reading multilink statistics for direct edges between all pairs of KEGG signaling pathways...')
        with open(proj_path + input_GRN + '_between_paths_multilink_counts_stats_discovery.pickle', 'rb') as fp:
            between_paths_multilink_counts_stats_df_dict = pickle.load(fp)  
        with open(proj_path + input_GRN + '_between_paths_rand_zscores_discovery.pickle', 'rb') as fp:
            between_paths_rand_zscores_df_dict = pickle.load(fp)          

    return (between_paths_multilink_counts_stats_df_dict, between_paths_rand_zscores_df_dict)

def get_shortest_paths_discovery(A_KEGGPPI_sparr, discovery_crosstalk_paths, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, proj_path, 
                       sp_threshold=None):

    if not (os.path.exists(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_edges_dict_discovery.pickle') &
            os.path.exists(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_intermediaries_dict_discovery.pickle')):     
    
        print('Calculating the shortest paths with threshold (%s) between all pairs of KEGG signaling pathways...' % sp_threshold)
        G_KEGGPPI = nx.from_scipy_sparse_matrix(A_KEGGPPI_sparr.astype(bool).astype(int))
        all_perms = [(i, j) for i, j in permutations(discovery_crosstalk_paths, 2)]

        shortest_path_intermediaries_dict = {}
        shortest_path_edges_dict = {}
        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            multilink_counts = {}   

            path_ix1 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p1]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values
            path_ix2 = pd.merge(KEGG_PPI_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p2]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values

            shortest_path_edges_dict[(p1, p2)] = set()
            shortest_path_intermediaries_dict[(p1, p2)] = set()

            if sp_threshold == None:
                for pp1 in np.setdiff1d(path_ix1, path_ix2):
                    for pp2 in np.setdiff1d(path_ix2, path_ix1):
                        if nx.has_path(G_KEGGPPI, source=pp1, target=pp2):
                            temp_sp = nx.shortest_path(G_KEGGPPI, source=pp1, target=pp2)
                            intermediary_nodes = list(set(temp_sp) - (set(path_ix1) | set(path_ix2))) # intermediary nodes that are a part of neither pathway
                            if len(intermediary_nodes) > 0:
                                inter_sp = intermediary_sp(temp_sp, intermediary_nodes)
                                shortest_path_edges_dict[(p1, p2)].update(set([(i, j) for i, j in zip(inter_sp, inter_sp[1:])]))
                                shortest_path_intermediaries_dict[(p1, p2)].update(set(intermediary_nodes))
            else:
                for pp1 in np.setdiff1d(path_ix1, path_ix2):
                    for pp2 in np.setdiff1d(path_ix2, path_ix1):
                        if nx.has_path(G_KEGGPPI, source=pp1, target=pp2):
                            temp_sp = nx.shortest_path(G_KEGGPPI, source=pp1, target=pp2)
                            intermediary_nodes = list(set(temp_sp) - (set(path_ix1) | set(path_ix2))) # intermediary nodes that are a part of neither pathway
                            if (len(intermediary_nodes) > 0) & (len(intermediary_nodes) <= sp_threshold):
                                inter_sp = intermediary_sp(temp_sp, intermediary_nodes)
                                shortest_path_edges_dict[(p1, p2)].update(set([(i, j) for i, j in zip(inter_sp, inter_sp[1:])]))
                                shortest_path_intermediaries_dict[(p1, p2)].update(set(intermediary_nodes))      

        with open(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_edges_dict_discovery.pickle', 'wb') as fp:
            pickle.dump(shortest_path_edges_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
        with open(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_intermediaries_dict_discovery.pickle', 'wb') as fp:
            pickle.dump(shortest_path_intermediaries_dict, fp, protocol=pickle.HIGHEST_PROTOCOL) 
            
    else:

        print('Reading the shortest paths with threshold (%s) between all pairs of KEGG signaling pathways...' % sp_threshold)
        with open(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_edges_dict_discovery.pickle', 'rb') as fp:
            shortest_path_edges_dict = pickle.load(fp)
        with open(proj_path + 'spthreshold' + str(sp_threshold) + '_shortest_path_intermediaries_dict_discovery.pickle', 'rb') as fp:
            shortest_path_intermediaries_dict = pickle.load(fp) 
                       
    return (shortest_path_edges_dict, shortest_path_intermediaries_dict)

### multilink counts of the shortest path edges between pairs of KEGG signaling pathways
# note that the below is for the smaller subset of pathways that are in XTalkDB for cross-validation purposes
def shortest_paths_multilink_counts_discovery(discovery_crosstalk_paths, shortest_path_edges_dict, sp_threshold,
                                    A_GRN_sparr, A_KEGGPPI_sparr, all_motif_types, all_motif_types_list, KEGG_interaction_types_dict, 
                                    A_KEGG_e_sparr_dict, proj_path, input_GRN):
    
    if not os.path.exists(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_discovery.pickle'):  
        
        print('Calculating multilink counts of shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways...' 
              % str(sp_threshold))   
        shortest_paths_multilink_counts_df_dict = {}

        all_perms = [(i, j) for i, j in permutations(discovery_crosstalk_paths, 2)]

        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            multilink_counts = {}   

            # fancy index the sparse arrays to get the shortest path edges between the pair of pathways
            A_2_arr = np.asarray(A_GRN_sparr[np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                             np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]
            A_1_arr = np.asarray(A_KEGGPPI_sparr[np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                                 np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]

            # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs [only for 0X]
            A_mux = np.dstack([A_1_arr, A_2_arr])
            e_multilink_counts = count_motifs_1D(A_mux, all_motif_types_list)

            for motif in e_multilink_counts.keys():
                if motif in ['0-1', '00', '01']:
                    multilink_counts[motif] = e_multilink_counts[motif]    


            for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):       

                # slice KEGG(e)+PPI to get the pathway nodes only
                A_1_arr = np.asarray(A_KEGG_e_sparr_dict[e][np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                                            np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                A_mux = np.dstack([A_1_arr, A_2_arr])
                e_multilink_counts = count_motifs_1D(A_mux, all_motif_types_list)

                # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                # (since we don't want to multiple-count these for every interaction type e)
                for motif in e_multilink_counts.keys():
                    if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                        multilink_counts[motif] = e_multilink_counts[motif]


            multilink_counts_df = pd.DataFrame(0, index=all_motif_types, columns=['motif_counts'])
            for key in multilink_counts.keys():
                if key in all_motif_types:
                    multilink_counts_df.at[key, 'motif_counts'] = multilink_counts[key]


            shortest_paths_multilink_counts_df_dict[(p1, p2)] = multilink_counts_df  

        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_discovery.pickle', 'wb') as fp:
            pickle.dump(shortest_paths_multilink_counts_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        print('Reading multilink counts of shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways...' 
              % str(sp_threshold))
        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_discovery.pickle', 'rb') as fp:
            shortest_paths_multilink_counts_df_dict = pickle.load(fp)       
            
    return shortest_paths_multilink_counts_df_dict

# The version where randomized instances are selected uniformly at random, as opposed to in order for reproducibility. The reason we need is that
# we need more than 250k (500 x 500) unique layer combinations in the full pathways case (3660 x 100). The output file has the suffix "discovery"
# to indicate that this is all the KEGG signaling pathways that we have (61 in total) and is not limited to the 25 benchmark pathways.
def shortest_paths_multilink_counts_rand_discovery(A_GRN_sparr_rand_dict, discovery_crosstalk_paths, shortest_path_edges_dict, sp_threshold,
                                        all_motif_types, A_KEGGPPI_sparr_rand_dict, all_motif_types_list, KEGG_interaction_types_dict, 
                                        A_KEGG_e_sparr_rand_dict, proj_path, input_GRN, N_rand = 100, N_swap = 10):
    
    if not os.path.exists(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_rand_%srands_discovery.pickle' % N_rand):     
    
        print('Calculating multilink counts of shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways for randomized multilayer network ensebles...' 
              % str(sp_threshold))  
        shortest_paths_multilink_counts_rand_df_dict = {}
        all_pairs = [(i, j) for i, j in product(list(np.arange(len(A_GRN_sparr_rand_dict))), list(np.arange(len(A_GRN_sparr_rand_dict))))]
        all_perms = [(i, j) for i, j in permutations(discovery_crosstalk_paths, 2)]
        all_pairs_shuff = [(i, j) for i, j in product(np.random.choice(len(A_GRN_sparr_rand_dict), 600), 
                                                      np.random.choice(len(A_GRN_sparr_rand_dict), 600))]
        
        for nn, (p1, p2) in enumerate(tqdm(all_perms, position=0, leave=True)): 

            multilink_counts_rand_df = pd.DataFrame(0, index=all_motif_types, columns=np.arange(N_rand))

            for nrand in np.arange(N_rand):    

                rand_ix = all_pairs_shuff[nn*N_rand+nrand]

                multilink_counts_rand = {}   

                # slice GRN to get the pathway nodes only
                A_2_arr_rand = np.asarray(A_GRN_sparr_rand_dict[rand_ix[1]][np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                                                            np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]  
                A_1_arr_rand = np.asarray(A_KEGGPPI_sparr_rand_dict[rand_ix[0]][np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                                                                np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]

                # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs [only for 0X]
                A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                e_multilink_counts_rand = count_motifs_1D(A_mux_rand, all_motif_types_list)

                for motif in e_multilink_counts_rand.keys():
                    if motif in ['0-1', '00', '01']:
                        multilink_counts_rand[motif] = e_multilink_counts_rand[motif]    


                for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):       

                    # slice KEGG(e)+PPI to get the pathway nodes only
                    A_1_arr_rand = np.asarray(A_KEGG_e_sparr_rand_dict[rand_ix[0]][e][np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0],
                                                                                      np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]

                    # generate the multiplex (based on the selected KEGG pathway's nodes) and count the multilink motifs
                    A_mux_rand = np.dstack([A_1_arr_rand, A_2_arr_rand])
                    e_multilink_counts_rand = count_motifs_1D(A_mux_rand, all_motif_types_list)

                    # disregard the motifs with no links in the signaling network (0X type motifs) since these will be counted separately below
                    # (since we don't want to multiple-count these for every interaction type e)
                    for motif in e_multilink_counts_rand.keys():
                        if (motif in all_motif_types) & (motif not in (['0-1', '00', '01'])):
                            multilink_counts_rand[motif] = e_multilink_counts_rand[motif]


                for key in multilink_counts_rand.keys():
                    if key in all_motif_types:
                        multilink_counts_rand_df.at[key, nrand] = multilink_counts_rand[key]

            shortest_paths_multilink_counts_rand_df_dict[(p1, p2)] = multilink_counts_rand_df     

        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_rand_%srands_discovery.pickle' % N_rand, 'wb') as fp:
            pickle.dump(shortest_paths_multilink_counts_rand_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
    else:
        
        print('Reading multilink counts of shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways for randomized multilayer network ensebles...' 
              % str(sp_threshold))
        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_rand_%srands_discovery.pickle' % N_rand, 'rb') as fp:
            shortest_paths_multilink_counts_rand_df_dict = pickle.load(fp)        
            
    return shortest_paths_multilink_counts_rand_df_dict

def shortest_paths_multilink_zscores_pvals_discovery(discovery_crosstalk_paths, shortest_paths_multilink_counts_rand_df_dict, sp_threshold,
                                          shortest_paths_multilink_counts_df_dict, all_motif_types, proj_path, input_GRN, N_rand=100):

    if not (os.path.exists(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_stats_discovery.pickle') &
            os.path.exists(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_rand_zscores_discovery.pickle')):         
    
        print('Calculating multilink statistics for shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways...' 
              % str(sp_threshold))
    
        shortest_paths_multilink_counts_stats_df_dict = {}
        shortest_paths_rand_zscores_df_dict = {}

        all_perms = [(i, j) for i, j in permutations(discovery_crosstalk_paths, 2)]

        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            shortest_paths_multilink_counts_stats_df = pd.DataFrame(index=all_motif_types)
            shortest_paths_multilink_counts_stats_df['rand_mean'] = shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].mean(axis=1)
            shortest_paths_multilink_counts_stats_df['rand_std'] = shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].std(axis=1)
            shortest_paths_multilink_counts_stats_df['actual_motif_counts'] = shortest_paths_multilink_counts_df_dict[(p1, p2)]['motif_counts']
            nostd_ix = shortest_paths_multilink_counts_stats_df[shortest_paths_multilink_counts_stats_df['rand_std'] != 0].index
            shortest_paths_multilink_counts_stats_df['z-score'] = ((shortest_paths_multilink_counts_stats_df.loc[nostd_ix]['actual_motif_counts'] - 
                                                                  shortest_paths_multilink_counts_stats_df.loc[nostd_ix]['rand_mean']) / 
                                                                  shortest_paths_multilink_counts_stats_df.loc[nostd_ix]['rand_std'])

            pvals_df = pd.DataFrame(index=all_motif_types, columns=['p-value'])
            for m in shortest_paths_multilink_counts_stats_df.index:
                if shortest_paths_multilink_counts_stats_df['z-score'].loc[m] > 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].loc[m] >= 
                                                     shortest_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif shortest_paths_multilink_counts_stats_df['z-score'].loc[m] < 0:
                    pvals_df.at[m, 'p-value'] = 1.0*(shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].loc[m] <= 
                                                     shortest_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m]).sum()/N_rand
                elif pd.isnull(shortest_paths_multilink_counts_stats_df['z-score'].loc[m]) & (shortest_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m] >
                                                                                            shortest_paths_multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0
                elif pd.isnull(shortest_paths_multilink_counts_stats_df['z-score'].loc[m]) & (shortest_paths_multilink_counts_stats_df['actual_motif_counts'].loc[m] <
                                                                                            shortest_paths_multilink_counts_stats_df['rand_mean'].loc[m]):
                    pvals_df.at[m, 'p-value'] = 0.0

            shortest_paths_multilink_counts_stats_df['p-value'] = pvals_df['p-value'] 

            shortest_paths_rand_zscores_df = shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].subtract(shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].mean(axis=1), axis=0).div(shortest_paths_multilink_counts_rand_df_dict[(p1, p2)].std(axis=1), axis=0)

            shortest_paths_multilink_counts_stats_df_dict[(p1, p2)] = shortest_paths_multilink_counts_stats_df
            shortest_paths_rand_zscores_df_dict[(p1, p2)] = shortest_paths_rand_zscores_df  

        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_stats_discovery.pickle', 'wb') as fp:
            pickle.dump(shortest_paths_multilink_counts_stats_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)  
        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_rand_zscores_discovery.pickle', 'wb') as fp:
            pickle.dump(shortest_paths_rand_zscores_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)  
            
    else:
        print('Reading multilink statistics for shortest path edges with threshold (%s) between all pairs of KEGG signaling pathways...' 
              % str(sp_threshold))
        
        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_multilink_counts_stats_discovery.pickle', 'rb') as fp:
            shortest_paths_multilink_counts_stats_df_dict = pickle.load(fp)  
        with open(proj_path + input_GRN + '_sp' + str(sp_threshold) + '_shortest_paths_rand_zscores_discovery.pickle', 'rb') as fp:
            shortest_paths_rand_zscores_df_dict = pickle.load(fp)          

    return (shortest_paths_multilink_counts_stats_df_dict, shortest_paths_rand_zscores_df_dict)

# Simple way to construct the adjacency matrix for only the KEGG signaling network for the sole purpose of implementing XTalk, 
# without having to worry about individual interaction types, etc. 
def process_KEGG(KEGG_all_edges_entrez, proj_path, input_filenames_dict):
    
    print('Generating adjacency matrix for the KEGG signaling network...')
    ### import the KEGG signaling network, remove the GErel (gene regulatory) type
    KEGG_all_nodes_df = pd.read_csv(proj_path + input_filenames_dict['KEGG_all_nodes_df'])
    KEGG_allnodes_entrez = set([int(x.split('hsa:')[1]) for x in 
                                list(KEGG_all_nodes_df[KEGG_all_nodes_df['Node_type']=='gene']['KEGG_name(s)_expanded'].values)])

    KEGG_allnodes_entrez_df = pd.DataFrame(sorted(list(KEGG_allnodes_entrez)))
    KEGG_allnodes_entrez_df['ix'] = pd.DataFrame(sorted(list(KEGG_allnodes_entrez))).index

    A_KEGG = pd.DataFrame(0, index=sorted(list(KEGG_allnodes_entrez)), columns=sorted(list(KEGG_allnodes_entrez)))
    for s, t in tqdm(KEGG_all_edges_entrez[['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y']].drop_duplicates().values, 
                     position=0, leave=True):
        if s != t:
            A_KEGG.at[s, t] = 1

    A_KEGG_sparr = sparse.csr_matrix(A_KEGG.to_numpy()) 
    
    return (KEGG_allnodes_entrez_df, A_KEGG_sparr)

def node_overlap_sig(KEGG_allnodes_entrez_df, KEGG_PPI_all_edges_entrez, KEGG_all_paths, common_crosstalk_paths, proj_path):

    if not os.path.exists(proj_path + 'node_edge_overlap_pairs_df.csv'):
        
        print('Calculating node and edge overlap significance for all pairs of KEGG signaling pathways...')
        KEGG_allnodes_entrez = set(KEGG_allnodes_entrez_df[0])

        KEGG_path_edges_df_dict = {}
        for p in KEGG_all_paths:
            KEGG_path_edges_df_dict[p] = KEGG_PPI_all_edges_entrez[KEGG_PPI_all_edges_entrez['Path_label']==p]
            # remove the few edges with edge_subtype=NaN:
            KEGG_path_edges_df_dict[p] = KEGG_path_edges_df_dict[p][~pd.isnull(KEGG_path_edges_df_dict[p]['Edge_subtype'])]

        node_overlap_counts_df = pd.DataFrame(index=KEGG_all_paths, columns=KEGG_all_paths)
        node_overlap_J_df = pd.DataFrame(index=KEGG_all_paths, columns=KEGG_all_paths)
        node_overlap_hyperP_df = pd.DataFrame(index=KEGG_all_paths, columns=KEGG_all_paths)

        edge_overlap_counts_df = pd.DataFrame(index=KEGG_all_paths, columns=KEGG_all_paths)
        edge_overlap_J_df = pd.DataFrame(index=KEGG_all_paths, columns=KEGG_all_paths)
        edge_overlap_hyperP_df = pd.DataFrame(index=KEGG_all_paths, columns=KEGG_all_paths)

        #num_all_KEGG_edges = len(KEGG_all_edges_entrez[['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y']].drop_duplicates())
        num_all_KEGG_edges = (len(KEGG_allnodes_entrez) * (len(KEGG_allnodes_entrez) - 1)) / 2.0

        for p1 in tqdm(KEGG_all_paths, position=0, leave=True):
            for p2 in KEGG_all_paths:
                node_overlap_counts_df.at[p1, p2] = len(KEGG_path_nodes_entrez_dict[p1] & KEGG_path_nodes_entrez_dict[p2])
                node_overlap_J_df.at[p1, p2] = len(KEGG_path_nodes_entrez_dict[p1] & KEGG_path_nodes_entrez_dict[p2]) /\
                                                    len(KEGG_path_nodes_entrez_dict[p1] | KEGG_path_nodes_entrez_dict[p2])
                A = len(KEGG_path_nodes_entrez_dict[p1] & KEGG_path_nodes_entrez_dict[p2])
                B = len(KEGG_path_nodes_entrez_dict[p1] - KEGG_path_nodes_entrez_dict[p2])
                C = len(KEGG_path_nodes_entrez_dict[p2] - KEGG_path_nodes_entrez_dict[p1])
                D = len(KEGG_allnodes_entrez) - len(KEGG_path_nodes_entrez_dict[p1] | KEGG_path_nodes_entrez_dict[p2])        
                node_overlap_hyperP_df.at[p1, p2] = st.fisher_exact([[A, B], [C, D]])[1]
                if st.fisher_exact([[A, B], [C, D]])[1] != 0.0: # cap p-values at 1e-323
                    node_overlap_hyperP_df.at[p1, p2] = st.fisher_exact([[A, B], [C, D]])[1]
                else:
                    node_overlap_hyperP_df.at[p1, p2] = 1e-323


                temp_edge_inner_df = pd.merge(KEGG_path_edges_df_dict[p1], KEGG_path_edges_df_dict[p2],
                                                on=['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y'], 
                                              how='inner').drop_duplicates(['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y'])
                temp_edge_outer_df = pd.merge(KEGG_path_edges_df_dict[p1], KEGG_path_edges_df_dict[p2],
                                                on=['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y'], 
                                              how='outer').drop_duplicates(['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y'])

                edge_overlap_counts_df.at[p1, p2] = len(temp_edge_inner_df)
                edge_overlap_J_df.at[p1, p2] = len(temp_edge_inner_df) / len(temp_edge_outer_df)

                A = len(temp_edge_inner_df)
                B = len(KEGG_path_edges_df_dict[p1]) - len(temp_edge_inner_df)
                C = len(KEGG_path_edges_df_dict[p2]) - len(temp_edge_inner_df)
                D = num_all_KEGG_edges - (A + B + C)
                if st.fisher_exact([[A, B], [C, D]])[1] != 0.0: # cap p-values at 1e-323
                    edge_overlap_hyperP_df.at[p1, p2] = st.fisher_exact([[A, B], [C, D]])[1]
                else:
                    edge_overlap_hyperP_df.at[p1, p2] = 1e-323

        node_overlap_hyperFDR_df = pd.DataFrame(np.reshape(multi.multipletests(node_overlap_hyperP_df.values.flatten(), method='fdr_bh')[1], 
                                                           (len(node_overlap_hyperP_df.index), len(node_overlap_hyperP_df.columns))), 
                                                index=node_overlap_hyperP_df.index, columns=node_overlap_hyperP_df.columns)

        edge_overlap_hyperFDR_df = pd.DataFrame(np.reshape(multi.multipletests(edge_overlap_hyperP_df.values.flatten(), method='fdr_bh')[1], 
                                                           (len(edge_overlap_hyperP_df.index), len(edge_overlap_hyperP_df.columns))), 
                                                index=edge_overlap_hyperP_df.index, columns=edge_overlap_hyperP_df.columns)


        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]
        node_edge_overlap_pairs_df = pd.DataFrame(index=['%s-%s' % (i, j) for i, j in all_perms], columns=['Node overlap hypergeometric FDR', 
                                                                                                           'Edge overlap hypergeometric FDR'])

        for p1, p2 in all_perms:
            node_edge_overlap_pairs_df.at['%s-%s' % (p1, p2), 'Node overlap hypergeometric FDR'] = node_overlap_hyperFDR_df.at[p1, p2]
            node_edge_overlap_pairs_df.at['%s-%s' % (p1, p2), 'Edge overlap hypergeometric FDR'] = edge_overlap_hyperFDR_df.at[p1, p2] 

        
        node_edge_overlap_pairs_df.to_csv(proj_path + 'node_edge_overlap_pairs_df.csv')
    
    else:
        print('Reading node and edge overlap significance for all pairs of KEGG signaling pathways...')
        node_edge_overlap_pairs_df = pd.read_csv(proj_path + 'node_edge_overlap_pairs_df.csv', index_col=0)
        
    return node_edge_overlap_pairs_df

def KEGG_between_edges(common_crosstalk_paths, KEGG_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, A_KEGG_sparr, proj_path):

    if not os.path.exists(proj_path + 'KEGG_between_edges.csv'):
        
        print('Calculating the number of direct edges between all pairs of KEGG pathways in the KEGG signaling network...')
        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]
        KEGG_between_edges_df = pd.DataFrame(index=['%s-%s' % (i, j) for i, j in all_perms], columns=['Edge overlap'])

        for p1, p2 in tqdm(all_perms, position=0, leave=True): 

            path_ix1 = pd.merge(KEGG_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p1]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values
            path_ix2 = pd.merge(KEGG_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p2]))), how='inner', 
                           left_on=0, right_on=0)['ix'].values

            A_1_arr = A_KEGG_sparr[tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()

            KEGG_between_edges_df.at['%s-%s' % (p1, p2), 'Edge overlap'] = A_1_arr.sum()
            
            KEGG_between_edges_df.to_csv(proj_path + 'KEGG_between_edges.csv')
            
    else:
        
        print('Reading the number of direct edges between all pairs of KEGG pathways in the KEGG signaling network...')
        KEGG_between_edges_df = pd.read_csv(proj_path + 'KEGG_between_edges.csv', index_col=0)
        
    return KEGG_between_edges_df

def KEGG_between_edges_rand(common_crosstalk_paths, KEGG_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, A_KEGG_sparr_rand_dict, proj_path):
    
    if not os.path.exists(proj_path + 'KEGG_between_edges_rand.csv'):
    
        print('Calculating the number of direct edges between all pairs of KEGG pathways in the KEGG signaling network for randomized network ensebles...')
        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]
        KEGG_between_edges_rand_df = pd.DataFrame(index=['%s-%s' % (i, j) for i, j in all_perms], 
                                                  columns=np.arange(len(A_KEGG_sparr_rand_dict.keys())))

        for n in tqdm(np.arange(len(A_KEGG_sparr_rand_dict)), position=0, leave=True):


            for p1, p2 in all_perms: 

                path_ix1 = pd.merge(KEGG_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p1]))), how='inner', 
                               left_on=0, right_on=0)['ix'].values
                path_ix2 = pd.merge(KEGG_allnodes_entrez_df, pd.DataFrame(sorted(list(KEGG_path_nodes_entrez_dict[p2]))), how='inner', 
                               left_on=0, right_on=0)['ix'].values

                A_1_arr_rand = A_KEGG_sparr_rand_dict[n][tuple([np.setdiff1d(path_ix1, path_ix2)])][:, np.setdiff1d(path_ix2, path_ix1)].toarray()

                KEGG_between_edges_rand_df.at['%s-%s' % (p1, p2), n] = A_1_arr_rand.sum()
                
        KEGG_between_edges_rand_df.to_csv(proj_path + 'KEGG_between_edges_rand.csv')
    
    else:
        
        print('Reading the number of direct edges between all pairs of KEGG pathways in the KEGG signaling network for randomized network ensebles...')
        KEGG_between_edges_rand_df = pd.read_csv(proj_path + 'KEGG_between_edges_rand.csv', index_col=0)
        
    return KEGG_between_edges_rand_df

@jit(nopython=True)
def edge_swap(adj, ntry):
    
    adj_rand = np.copy(adj)
    nrew = 0
    ix1 = np.where(adj_rand>0)[0]
    ix2 = np.where(adj_rand>0)[1]

    # swap ntry * # of edges
    for i in np.arange(ntry*len(ix1)):
        
        # choose two edges at random
        swap_ix1 = int(random.random() * len(ix1)) #random.sample(set(np.arange(len(ix1))), 1)[0]
        swap_ix2 = int(random.random() * len(ix1)) #random.sample(set(np.arange(len(ix1))), 1)[0]
        A = ix1[swap_ix1]
        B = ix2[swap_ix1]
        C = ix1[swap_ix2]
        D = ix2[swap_ix2]

        if (len(set([A, B, C, D])) == 4):
            if ((adj_rand[A, D] == 0) & (adj_rand[C, B] == 0)):

                AB_val = adj_rand[A, B]
                CD_val = adj_rand[C, D]
                
                adj_rand[A, B] = 0
                adj_rand[C, D] = 0                
                adj_rand[A, D] = AB_val
                adj_rand[C, B] = CD_val   
                nrew = nrew + 1
                ix1[swap_ix1] = A
                ix2[swap_ix1] = D
                ix1[swap_ix2] = C
                ix2[swap_ix2] = B                
            
    return (adj_rand, nrew)      

def randomize_KEGG(A_KEGG_sparr, proj_path, N_runs=500, N_swap=10, return_output=False):

    if not os.path.exists(proj_path + 'A_KEGG_sparr_rand_dict_%sruns.pickle' % N_runs):    
    
        print('Computing ensemble of randomized KEGG sparse  matrices...')    
    
        A_KEGG_sparr_rand_dict = {}

        for nrand in tqdm(np.arange(N_runs), position=0, leave=True):
            A_KEGG_arr_rand, A_KEGG_nrew_rand = edge_swap(A_KEGG_sparr.toarray(), N_swap)
            A_KEGG_sparr_rand_dict[nrand] = sparse.csr_matrix(A_KEGG_arr_rand)
 
        with open(proj_path + 'A_KEGG_sparr_rand_dict_%sruns.pickle' % N_runs, 'wb') as fp:
            pickle.dump(A_KEGG_sparr_rand_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
        
        if return_output==True:        
            return A_KEGG_sparr_rand_dict
        
    else:
        
        print('Reading ensemble of randomized KEGG sparse  matrices...')    
        
        with open(proj_path + 'A_KEGG_sparr_rand_dict_%sruns.pickle' % N_runs, 'rb') as fp:
            A_KEGG_sparr_rand_dict = pickle.load(fp)
        
        if return_output==True:        
            return A_KEGG_sparr_rand_dict        

def process_receptor_TF_data(proj_path, input_filenames_dict):
    
    HUGO_symb_entrez_uniprot = pd.read_csv(proj_path + input_filenames_dict['HUGO_symb_entrez_uniprot'], sep='\t')
    
    # Get receptors from [Almén, Markus Sällman, et al. "Mapping the human membrane proteome: a majority of the human membrane proteins can be 
    # classified according to function and evolutionary origin." BMC biology 7.1 (2009): 1-14.].
    # Parse gene symbols from IPI Descriptions, match with recent annotations from HUGO Gene Names.
    # Match remaining old symbols with Previous Symbols from HUGO Genes.
    # In the end, only 1211 genes out of 1350 are matched to their Entrez IDs via their most recent gene symbols.
    
    print('Obtaining cell surface receptor information from Almen et al. BMC biology 7.1 (2009): 1-14. ...')
    almen_etal = pd.read_csv(proj_path + input_filenames_dict['almen_etal'])
    almen_etal_rec = almen_etal[almen_etal['Main Class'] == 'Receptors']
    print(len(almen_etal_rec['IPI Accession'].unique()))
    almen_etal_rec = almen_etal_rec.assign(gene_symbol=almen_etal_rec['IPI Description'].str.split(' ').str[0])
    print(len(almen_etal_rec['gene_symbol'].unique()))

    unmapped_symb = set(almen_etal_rec['gene_symbol'].unique()) - set(HUGO_symb_entrez_uniprot['Approved symbol']) - set(['-'])
    print(len(unmapped_symb))

    HUGO_prev_symb = HUGO_symb_entrez_uniprot[~pd.isnull(HUGO_symb_entrez_uniprot['Previous symbols'])]

    old_new_symb_dict = {}
    for i in unmapped_symb:
        if len( HUGO_prev_symb[HUGO_prev_symb['Previous symbols'].str.contains(i)]) > 0:
            old_new_symb_dict[i] = HUGO_prev_symb[HUGO_prev_symb['Previous symbols'].str.contains(i)]['Approved symbol'].values[0]

    print(len(old_new_symb_dict))

    for i in old_new_symb_dict.keys():
        old_ix = almen_etal_rec[almen_etal_rec['gene_symbol'] == i].index[0]
        almen_etal_rec.at[old_ix, 'gene_symbol'] = old_new_symb_dict[i]

    unmapped_symb = set(almen_etal_rec['gene_symbol'].unique()) - set(HUGO_symb_entrez_uniprot['Approved symbol']) - set(['-'])
    print(len(unmapped_symb))

    almen_etal_rec_entrez = pd.merge(almen_etal_rec, HUGO_symb_entrez_uniprot, left_on='gene_symbol', right_on='Approved symbol')
    almen_etal_rec_entrez = almen_etal_rec_entrez[~pd.isnull(almen_etal_rec_entrez['NCBI Gene ID(supplied by NCBI)'])]
    almen_etal_rec_entrez = almen_etal_rec_entrez.astype({'NCBI Gene ID(supplied by NCBI)': int})

    almen_etal_rec_entrez_set = set(almen_etal_rec_entrez['NCBI Gene ID(supplied by NCBI)'].unique()) - set(['-'])
    print(len(almen_etal_rec_entrez_set))
    
    # Get TFs from [Lambert, Samuel A., et al. "The human transcription factors." Cell 172.4 (2018): 650-665.]
    # Convert to entrez ids
    print('Obtaining transcription factor information from Lambert et al. Cell 172.4 (2018): 650-665. ...')
    lambert_etal = pd.read_csv(proj_path + input_filenames_dict['lambert_etal'])
    lambert_etal_TFs = lambert_etal[lambert_etal['Is TF?']=='Yes']
    lambert_etal_TFs_entrez = pd.merge(lambert_etal_TFs, HUGO_symb_entrez_uniprot, left_on='Unnamed: 1', right_on='Approved symbol')
    lambert_etal_TFs_entrez = lambert_etal_TFs_entrez.astype({'NCBI Gene ID(supplied by NCBI)': int})

    lambert_etal_TFs_entrez_set = set(lambert_etal_TFs_entrez['NCBI Gene ID(supplied by NCBI)'])
    print(len(lambert_etal_TFs_entrez_set))  
    
    return (almen_etal_rec_entrez_set, lambert_etal_TFs_entrez_set)

def build_diG_KEGG_forXTalk(A_KEGG_sparr, KEGG_allnodes_entrez_df):
    
    KEGG_edgelist_df = pd.DataFrame(columns=['node1_ix', 'node2_ix'])
    KEGG_edgelist_df['node1_ix'], KEGG_edgelist_df['node2_ix'] = np.nonzero(A_KEGG_sparr.toarray())

    KEGG_edgelist_entrez_df = pd.merge(pd.merge(KEGG_edgelist_df, KEGG_allnodes_entrez_df, left_on='node1_ix', right_on='ix', how='left'), 
                                       KEGG_allnodes_entrez_df, left_on='node2_ix', right_on='ix', how='left')[['0_x', '0_y']]

    diG_KEGG_forXTalk = nx.DiGraph()
    diG_KEGG_forXTalk.add_edges_from(KEGG_edgelist_entrez_df.to_numpy())

    return diG_KEGG_forXTalk

def XTalk_chi(common_crosstalk_paths, KEGG_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, A_KEGG_sparr, almen_etal_rec_entrez_set, 
              lambert_etal_TFs_entrez_set, proj_path, K):

    if not os.path.exists(proj_path + 'XTalk_chi_df_%sK.csv' % K):
        
        print('Calculating chi statistic (XTalk) with K = %s for all pairs of KEGG pathways... ' % K)
        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]
        XTalk_chi_df = pd.DataFrame(index=['%s-%s' % (i, j) for i, j in all_perms], columns=['XTalk_Chi'])

        diG_KEGG_forXTalk = build_diG_KEGG_forXTalk(A_KEGG_sparr, KEGG_allnodes_entrez_df)

        for A, B in tqdm(all_perms, position=0, leave=True): 
            rec_A = list(set(KEGG_path_nodes_entrez_dict[A]) & almen_etal_rec_entrez_set)
            TF_B = list(set(KEGG_path_nodes_entrez_dict[B]) & lambert_etal_TFs_entrez_set)
            if (len(rec_A) > 0) & (len(TF_B) > 0):
                d_rec_A = {}
                for t in TF_B:
                    d_rec_A[t] = []
                    for r in rec_A:
                        if nx.has_path(diG_KEGG_forXTalk, source=r, target=t):        
                            all_paths_gen = nx.shortest_simple_paths(diG_KEGG_forXTalk, source=r, target=t)
                            for path in islice(all_paths_gen, K):
                                d_rec_A[t].append(len(path) - 1)
                    
                TF_B_sum = []
                for t in d_rec_A.keys():
                    d_rec_A[t] = sorted(d_rec_A[t])
                    TF_B_sum.append(sum(d_rec_A[t][0:K]))

                if sum(TF_B_sum) > 0:
                    XTalk_chi_df.at['%s-%s' % (A, B), 'XTalk_Chi'] = (1/(K*len(TF_B))) * sum(TF_B_sum)            
                elif sum(TF_B_sum) == 0: # this sum being zero means there is no path between any receptor in A and TF t in B
                    XTalk_chi_df.at['%s-%s' % (A, B), 'XTalk_Chi'] = np.nan
            else:
                XTalk_chi_df.at['%s-%s' % (A, B), 'XTalk_Chi'] = np.nan

        #XTalk_chi_df = XTalk_chi_df[~pd.isnull(XTalk_chi_df['XTalk_Chi'])].sort_values('XTalk_Chi')

        XTalk_chi_df.to_csv(proj_path + 'XTalk_chi_df_%sK.csv' % K)
        
    else:
        
        print('Reading chi statistic (XTalk) with K = %s for all pairs of KEGG pathways... ' % K)
        XTalk_chi_df = pd.read_csv(proj_path + 'XTalk_chi_df_%sK.csv' % K, index_col=0)
           
    return XTalk_chi_df

def XTalk_chi_rand(common_crosstalk_paths, KEGG_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, A_KEGG_sparr_rand_dict, almen_etal_rec_entrez_set, 
                   lambert_etal_TFs_entrez_set, proj_path, K, Nrand):

    if not os.path.exists(proj_path + 'XTalk_chi_rand_df_%sK_%sNrand.csv' % (K, Nrand)):
    
        print('Calculating chi statistic (XTalk) with K = %s for all pairs of KEGG pathways for randomized network ensebles... ' % K)
        all_perms = [(i, j) for i, j in permutations(common_crosstalk_paths, 2)]
        XTalk_chi_rand_df = pd.DataFrame(index=['%s-%s' % (i, j) for i, j in all_perms], columns=np.arange(Nrand))

        for n in tqdm(np.arange(Nrand), position=0, leave=True):

            diG_KEGG_forXTalk_rand = build_diG_KEGG_forXTalk(A_KEGG_sparr_rand_dict[n], KEGG_allnodes_entrez_df)

            for A, B in all_perms: 
                rec_A = list(set(KEGG_path_nodes_entrez_dict[A]) & almen_etal_rec_entrez_set)
                TF_B = list(set(KEGG_path_nodes_entrez_dict[B]) & lambert_etal_TFs_entrez_set)
                if (len(rec_A) > 0) & (len(TF_B) > 0):
                    d_rec_A = {}
                    for t in TF_B:
                        d_rec_A[t] = []
                        for r in rec_A:
                            if nx.has_path(diG_KEGG_forXTalk_rand, source=r, target=t):      
                                all_paths_gen = nx.shortest_simple_paths(diG_KEGG_forXTalk_rand, source=r, target=t)
                                for path in islice(all_paths_gen, K):
                                    d_rec_A[t].append(len(path) - 1)
                                    
                    TF_B_sum = []
                    for t in d_rec_A.keys():
                        d_rec_A[t] = sorted(d_rec_A[t])
                        TF_B_sum.append(sum(d_rec_A[t][0:K]))

                    if sum(TF_B_sum) > 0:
                        XTalk_chi_rand_df.at['%s-%s' % (A, B), n] = (1/(K*len(TF_B))) * sum(TF_B_sum)            
                    elif sum(TF_B_sum) == 0: # this sum being zero means there is no path between any receptor in A and TF t in B
                        XTalk_chi_rand_df.at['%s-%s' % (A, B), n] = np.nan
                else:
                    XTalk_chi_rand_df.at['%s-%s' % (A, B), n] = np.nan

        XTalk_chi_rand_df.to_csv(proj_path + 'XTalk_chi_rand_df_%sK_%sNrand.csv' % (K, Nrand))
    
    else:
        
        print('Reading chi statistic (XTalk) with K = %s for all pairs of KEGG pathways for randomized network ensebles... ' % K)
        XTalk_chi_rand_df = pd.read_csv(proj_path + 'XTalk_chi_rand_df_%sK_%sNrand.csv' % (K, Nrand), index_col=0)
    
    return XTalk_chi_rand_df

def multilink_positives(multilink_counts_stats_df_dict):
    
    grn_edge_motifs = np.array([[str(x)+'1' for x in list(range(12))], [str(x)+'-1' for x in list(range(12))]]).ravel()

    all_perms = [(i, j) for i, j in multilink_counts_stats_df_dict.keys()]
    multilink_pos = {i: [] for i in grn_edge_motifs}

    for (p1, p2) in all_perms: 

        grn_edge_motifs_zscores = multilink_counts_stats_df_dict[(p1, p2)].loc[grn_edge_motifs]
        grn_edge_motifs_zscores_sig = grn_edge_motifs_zscores[(grn_edge_motifs_zscores['p-value'] <= 0.05) & 
                                                              (grn_edge_motifs_zscores['actual_motif_counts'] - grn_edge_motifs_zscores['rand_mean'] > 0)]
        for m in grn_edge_motifs_zscores_sig.index.to_numpy():
            multilink_pos[m].append('%s-%s' % (p1, p2)) 

    for m in grn_edge_motifs:
        if len(multilink_pos[m]) == 0:
            multilink_pos.pop(m)

    multilink_pos_all = []
    for m in multilink_pos.keys():
        multilink_pos_all.extend(multilink_pos[m])
    multilink_pos_all = set(multilink_pos_all)

    return (multilink_pos, multilink_pos_all)

def XTalk_positives(XTalk_chi_df, XTalk_chi_rand_df):
    XTalk_chi_stats_df = pd.DataFrame(index=XTalk_chi_df.index, columns=['z-score', 'emp. pval'])
    XTalk_chi_stats_df['z-score'] = (XTalk_chi_df['XTalk_Chi'] - XTalk_chi_rand_df.mean(axis=1)) / XTalk_chi_rand_df.std(axis=1)
    XTalk_chi_stats_df['emp. pval'] = XTalk_chi_rand_df.lt(XTalk_chi_df['XTalk_Chi'], axis=0).sum(axis=1) / XTalk_chi_rand_df.shape[1]
    XTalk_pos_all = set(XTalk_chi_stats_df[(XTalk_chi_stats_df['emp. pval'] <= 0.05) & (~pd.isnull(XTalk_chi_stats_df['z-score']))].index)

    XTalk_nopath_pairs = XTalk_chi_df[pd.isnull(XTalk_chi_df['XTalk_Chi'])].index

    return (XTalk_chi_stats_df, XTalk_pos_all, XTalk_nopath_pairs)

def KEGG_between_edges_positives(KEGG_between_edges_df, KEGG_between_edges_rand_df):
    KEGG_between_edges_stats_df = pd.DataFrame(index=KEGG_between_edges_df.index, columns=['z-score', 'emp. pval'])
    KEGG_between_edges_stats_df['z-score'] = (KEGG_between_edges_df['Edge overlap'] - 
                                              KEGG_between_edges_rand_df.mean(axis=1)) / KEGG_between_edges_rand_df.std(axis=1)
    KEGG_between_edges_stats_df['emp. pval'] = KEGG_between_edges_rand_df.gt(KEGG_between_edges_df['Edge overlap'], axis=0).sum(axis=1) / \
                                                    KEGG_between_edges_rand_df.shape[1]
    KEGG_between_edges_pos_all = set(KEGG_between_edges_stats_df[(KEGG_between_edges_stats_df['emp. pval'] <= 0.05) & 
                                                                 (~pd.isnull(KEGG_between_edges_stats_df['z-score']))].index)

    #KEGG_between_edges_nopath_pairs = KEGG_between_edges_df[pd.isnull(KEGG_between_edges_df['Edge overlap'])].index

    return (KEGG_between_edges_stats_df, KEGG_between_edges_pos_all)#, KEGG_between_edges_nopath_pairs)

def calculate_auroc(TPR, FPR):
    
    dTPR = np.concatenate((np.ediff1d(TPR), [0]))
    dFPR = np.concatenate((np.ediff1d(FPR), [0]))
    
    return sum(TPR * dFPR) + sum(dTPR * dFPR)/2

def detected_auroc_auprc(detected_ranked_pathway_pairs, common_crosstalk_paths, XTalk_DB_positives_common_str):
    
    detected_roc_pr_stats = {}
    
    binsize = 1
    bins = np.arange(1, len(detected_ranked_pathway_pairs), binsize)      
    TPR = np.zeros(len(bins))
    FPR = np.zeros(len(bins))
    precision = np.zeros(len(bins))
    recall = np.zeros(len(bins))    

    all_perms_str = ['%s-%s'%(i, j) for i, j in permutations(common_crosstalk_paths, 2)]
    updated_universe = set(detected_ranked_pathway_pairs)
    updated_XTalk_DB_positives_common_str = set(XTalk_DB_positives_common_str) & set(detected_ranked_pathway_pairs)
    # calculate ROC and PR curves
    for i, n in enumerate(bins):

        TP = len(set(detected_ranked_pathway_pairs[0:n]) & updated_XTalk_DB_positives_common_str)
        FP = len(set(detected_ranked_pathway_pairs[0:n]) - updated_XTalk_DB_positives_common_str)
        FN = len(updated_XTalk_DB_positives_common_str - set(detected_ranked_pathway_pairs[0:n]))
        TN = len(updated_universe  - (set(detected_ranked_pathway_pairs[0:n]) | updated_XTalk_DB_positives_common_str))
        TPR[i] = TP / (TP + FN)
        FPR[i] = FP / (FP + TN)
        precision[i] = TP / (TP + FP)
        recall[i] = TP / (TP + FN)

    detected_roc_pr_stats['TPR'] = TPR
    detected_roc_pr_stats['FPR'] = FPR
    detected_roc_pr_stats['precision'] = precision
    detected_roc_pr_stats['recall'] = recall
    detected_roc_pr_stats['AUROC'] = calculate_auroc(TPR, FPR)
    detected_roc_pr_stats['AUPRC'] = calculate_auroc(precision, recall)
   
    return detected_roc_pr_stats

def auroc_auprc(ranked_pathway_pairs, ranked_pathway_pairs_df, len_sig, common_crosstalk_paths, XTalk_DB_positives_common_str, nonsig_Nrand = 1000):

    roc_pr_stats = {}
    
    binsize = 1
    bins = np.arange(1, len(ranked_pathway_pairs_df), binsize) 
    TPR = np.zeros((len(bins), nonsig_Nrand))
    FPR = np.zeros((len(bins), nonsig_Nrand))
    precision = np.zeros((len(bins), nonsig_Nrand))
    recall = np.zeros((len(bins), nonsig_Nrand))  
    AUROC = np.zeros(nonsig_Nrand)
    AUPRC = np.zeros(nonsig_Nrand)

    for nn in tqdm(np.arange(nonsig_Nrand), position=0, leave=True):

        ranked_pathway_pairs_shuf = np.concatenate((ranked_pathway_pairs[0:len_sig], 
                                                    random.sample(list(ranked_pathway_pairs[len_sig:]), len(ranked_pathway_pairs[len_sig:]))))

        all_perms_str = ['%s-%s'%(i, j) for i, j in permutations(common_crosstalk_paths, 2)]
        universe = set(ranked_pathway_pairs_shuf)

        # calculate ROC and PR curves
        for i, n in enumerate(bins):

            TP = len(set(ranked_pathway_pairs_shuf[0:n]) & set(XTalk_DB_positives_common_str))
            FP = len(set(ranked_pathway_pairs_shuf[0:n]) - set(XTalk_DB_positives_common_str))
            FN = len(set(XTalk_DB_positives_common_str) - set(ranked_pathway_pairs_shuf[0:n]))
            TN = len(universe  - (set(ranked_pathway_pairs_shuf[0:n]) | set(XTalk_DB_positives_common_str)))
            TPR[i, nn] = TP / (TP + FN)
            FPR[i, nn] = FP / (FP + TN)
            precision[i, nn] = TP / (TP + FP)
            recall[i, nn] = TP / (TP + FN)

        AUROC[nn] = calculate_auroc(TPR[:, nn], FPR[:, nn])
        AUPRC[nn] = calculate_auroc(precision[:, nn], recall[:, nn])

    roc_pr_stats['TPR'] = TPR
    roc_pr_stats['FPR'] = FPR
    roc_pr_stats['TPR_lower'] = np.mean(TPR, axis=1) - np.std(TPR, axis=1)
    roc_pr_stats['TPR_upper'] = np.mean(TPR, axis=1) + np.std(TPR, axis=1)
    roc_pr_stats['TPR_mean'] = np.mean(TPR, axis=1)
    roc_pr_stats['FPR_mean'] = np.mean(FPR, axis=1)
    roc_pr_stats['precision'] = precision
    roc_pr_stats['recall'] = recall
    roc_pr_stats['precision_lower'] = np.mean(precision, axis=1) - np.std(precision, axis=1)
    roc_pr_stats['precision_upper'] = np.mean(precision, axis=1) + np.std(precision, axis=1)
    roc_pr_stats['precision_mean'] = np.mean(precision, axis=1)
    roc_pr_stats['recall_mean'] = np.mean(recall, axis=1)
    roc_pr_stats['AUROC_mean'] = np.mean(AUROC)
    roc_pr_stats['AUROC_std'] = np.std(AUROC)
    roc_pr_stats['AUPRC_mean'] = np.mean(AUPRC)
    roc_pr_stats['AUPRC_std'] = np.std(AUPRC)   
    roc_pr_stats['len_sig'] = len_sig

    return roc_pr_stats

def get_ranked_pathway_pairs_benchmark(XTalk_DB_positives_common_str, common_crosstalk_paths, proj_path, method='MuXTalk_shortest',
                             multilink_params={'GRN': 'HumanGRN10e6', 'sp_threshold': 1}, 
                             XTalk_params={'K': 1, 'Nrand': 100}):
    
    if method == 'MuXTalk_between':
        
        with open(proj_path + multilink_params['GRN'] + '_between_paths_multilink_counts_stats.pickle', 'rb') as fp:
            between_paths_multilink_counts_stats_df_dict = pickle.load(fp)
        (between_multilink_pos, between_multilink_pos_all) = multilink_positives(between_paths_multilink_counts_stats_df_dict)

        all_perms = [(i, j) for i, j in between_paths_multilink_counts_stats_df_dict.keys()]
        num_sig_multilinks = pd.DataFrame(index=['%s-%s' % (i, j) for i, j in all_perms], columns=['Pathway A', 'Pathway B', '# sig', 
                                                                                                   'sig multilink types',  'actual motif counts',
                                                                                                   'z-scores', 'best z-score', 'p-values',
                                                                                                  'best p-value', 'MuXTalk_score', 'in gold'])

        for (p1, p2) in all_perms: 
            temp_df = between_paths_multilink_counts_stats_df_dict[(p1, p2)].loc[between_multilink_pos.keys()]
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'Pathway A'] = p1
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'Pathway B'] = p2
            num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] = len(temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                             (temp_df['p-value']<=0.05)])
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'sig multilink types'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                       (temp_df['p-value']<=0.05)].index.to_numpy()
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'actual motif counts'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                       (temp_df['p-value']<=0.05)]['actual_motif_counts'].to_numpy() 
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'z-scores'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                            (temp_df['p-value']<=0.05)]['z-score'].to_numpy()
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                (temp_df['p-value']<=0.05)]['z-score'].max()
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'p-values'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                (temp_df['p-value']<=0.05)]['p-value'].to_numpy()
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                (temp_df['p-value']<=0.05)]['p-value'].min()
                
            if (num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] > 0) & (~np.isnan(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'])):
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = ((num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] * 1000) +
                                                                              num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'] *
                                                                              -np.log10(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] + 0.001))
            elif (num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] > 0) & (np.isnan(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'])):
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = ((num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] * 1000) +
                                                                              -np.log10(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] + 0.001))          
            else:
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = 0.0

            if '%s-%s' % (p1, p2) in XTalk_DB_positives_common_str:
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'in gold'] = 'Yes'
            else:
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'in gold'] = 'No'

        ranked_pathway_pairs_df = num_sig_multilinks.sort_values(['MuXTalk_score'], ascending=False)
        ranked_pathway_pairs = ranked_pathway_pairs_df.index
        detected_ranked_pathway_pairs_df = ranked_pathway_pairs_df[ranked_pathway_pairs_df['# sig'] > 0]
        detected_ranked_pathway_pairs = detected_ranked_pathway_pairs_df.index   
        len_sig = len(ranked_pathway_pairs_df[ranked_pathway_pairs_df['# sig'] > 0])
        
        
    elif method == 'MuXTalk_shortest':
        
        with open(proj_path + multilink_params['GRN'] + 
                  '_sp%s_shortest_paths_multilink_counts_stats.pickle' % multilink_params['sp_threshold'], 'rb') as fp:
            shortest_paths_multilink_counts_stats_df_dict = pickle.load(fp)
        (shortest_multilink_pos, shortest_multilink_pos_all) = multilink_positives(shortest_paths_multilink_counts_stats_df_dict)

        all_perms = [(i, j) for i, j in shortest_paths_multilink_counts_stats_df_dict.keys()]
        num_sig_multilinks = pd.DataFrame(index=['%s-%s' % (i, j) for i, j in all_perms], columns=['Pathway A', 'Pathway B', '# sig', 
                                                                                                   'sig multilink types',  'actual motif counts',
                                                                                                   'z-scores', 'best z-score', 'p-values',
                                                                                                  'best p-value', 'MuXTalk_score', 'in gold'])

        for (p1, p2) in all_perms: 
            temp_df = shortest_paths_multilink_counts_stats_df_dict[(p1, p2)].loc[shortest_multilink_pos.keys()]
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'Pathway A'] = p1
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'Pathway B'] = p2
            num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] = len(temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                             (temp_df['p-value']<=0.05)])
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'sig multilink types'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                       (temp_df['p-value']<=0.05)].index.to_numpy()
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'actual motif counts'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                       (temp_df['p-value']<=0.05)]['actual_motif_counts'].to_numpy() 
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'z-scores'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                            (temp_df['p-value']<=0.05)]['z-score'].to_numpy()
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                (temp_df['p-value']<=0.05)]['z-score'].max()
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'p-values'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                (temp_df['p-value']<=0.05)]['p-value'].to_numpy()
            num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                (temp_df['p-value']<=0.05)]['p-value'].min()

            if (num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] > 0) & (~np.isnan(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'])):
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = ((num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] * 1000) +
                                                                              num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'] *
                                                                              -np.log10(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] + 0.001))
            elif (num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] > 0) & (np.isnan(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'])):
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = ((num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] * 1000) +
                                                                              -np.log10(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] + 0.001))          
            else:
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = 0.0

            if '%s-%s' % (p1, p2) in XTalk_DB_positives_common_str:
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'in gold'] = 'Yes'
            else:
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'in gold'] = 'No'

        ranked_pathway_pairs_df = num_sig_multilinks.sort_values(['MuXTalk_score'], ascending=False)
        ranked_pathway_pairs = ranked_pathway_pairs_df.index
        detected_ranked_pathway_pairs_df = ranked_pathway_pairs_df[ranked_pathway_pairs_df['# sig'] > 0]
        detected_ranked_pathway_pairs = detected_ranked_pathway_pairs_df.index   
        len_sig = len(ranked_pathway_pairs_df[ranked_pathway_pairs_df['# sig'] > 0])
        
    elif method == 'XTalk':
        
        XTalk_chi_df = pd.read_csv(proj_path + 'XTalk_chi_df_%sK.csv' % XTalk_params['K'], index_col=0)
        XTalk_chi_rand_df = pd.read_csv(proj_path + 'XTalk_chi_rand_df_%sK_%sNrand.csv' % (XTalk_params['K'], XTalk_params['Nrand']), index_col=0)
        (XTalk_chi_stats_df, XTalk_pos_all, XTalk_nopath_pairs) = XTalk_positives(XTalk_chi_df, XTalk_chi_rand_df)

        ranked_pathway_pairs_df = XTalk_chi_stats_df.sort_values(['emp. pval', 'z-score'], ascending=[True, True])
        ranked_pathway_pairs = ranked_pathway_pairs_df.index
        detected_ranked_pathway_pairs_df = ranked_pathway_pairs_df[~pd.isnull(ranked_pathway_pairs_df['z-score'])]
        detected_ranked_pathway_pairs = detected_ranked_pathway_pairs_df.index  
        len_sig = len(ranked_pathway_pairs_df[~pd.isnull(ranked_pathway_pairs_df['z-score'])])
        
    elif method == 'node_overlap':
        
        node_edge_overlap_pairs_df = pd.read_csv(proj_path + 'node_edge_overlap_pairs_df.csv', index_col=0)
        ranked_pathway_pairs_df = node_edge_overlap_pairs_df.sort_values('Node overlap hypergeometric FDR')
        ranked_pathway_pairs = ranked_pathway_pairs_df.index
        detected_ranked_pathway_pairs_df = ranked_pathway_pairs_df[ranked_pathway_pairs_df['Node overlap hypergeometric FDR'] != 1.0]
        detected_ranked_pathway_pairs = detected_ranked_pathway_pairs_df.index 
        len_sig = len(ranked_pathway_pairs_df[ranked_pathway_pairs_df['Node overlap hypergeometric FDR'] != 1.0])
        
    elif method == 'edge_overlap':
        
        node_edge_overlap_pairs_df = pd.read_csv(proj_path + 'node_edge_overlap_pairs_df.csv', index_col=0)
        ranked_pathway_pairs_df = node_edge_overlap_pairs_df.sort_values('Edge overlap hypergeometric FDR')
        ranked_pathway_pairs = ranked_pathway_pairs_df.index
        detected_ranked_pathway_pairs_df = ranked_pathway_pairs_df[ranked_pathway_pairs_df['Edge overlap hypergeometric FDR'] != 1.0]
        detected_ranked_pathway_pairs = detected_ranked_pathway_pairs_df.index
        len_sig = len(ranked_pathway_pairs_df[ranked_pathway_pairs_df['Edge overlap hypergeometric FDR'] != 1.0])
        
    elif method == 'direct_edges':  
        
        KEGG_between_edges_df = pd.read_csv(proj_path + 'KEGG_between_edges.csv', index_col=0)
        KEGG_between_edges_rand_df = pd.read_csv(proj_path + 'KEGG_between_edges_rand.csv', index_col=0)
        (KEGG_between_edges_stats_df, KEGG_between_edges_pos_all) = KEGG_between_edges_positives(KEGG_between_edges_df, KEGG_between_edges_rand_df)

        ranked_pathway_pairs_df = KEGG_between_edges_stats_df.sort_values(['emp. pval', 'z-score'], ascending=[True, False])
        ranked_pathway_pairs = ranked_pathway_pairs_df.index
        detected_ranked_pathway_pairs_df = ranked_pathway_pairs_df
        detected_ranked_pathway_pairs = detected_ranked_pathway_pairs_df.index
        len_sig = len(ranked_pathway_pairs_df)

        

    return (detected_ranked_pathway_pairs_df, detected_ranked_pathway_pairs, ranked_pathway_pairs_df, ranked_pathway_pairs, len_sig)      


def get_ranked_pathway_pairs_discovery(XTalk_DB_positives_common_str, proj_path, method='MuXTalk_shortest', 
                                       multilink_params={'GRN': 'HumanGRN10e6', 'sp_threshold': 1}, parquet=False):
    
    if method == 'MuXTalk_between':
        
        if not os.path.exists(proj_path + multilink_params['GRN'] + '_between_detected_discovery.parquet'):
        
            with open(proj_path + multilink_params['GRN'] + 
                      '_between_paths_multilink_counts_stats_discovery.pickle', 'rb') as fp:
                between_paths_multilink_counts_stats_df_dict = pickle.load(fp)
            (between_multilink_pos, between_multilink_pos_all) = multilink_positives(between_paths_multilink_counts_stats_df_dict)

            all_perms = [(i, j) for i, j in between_paths_multilink_counts_stats_df_dict.keys()]
            num_sig_multilinks = pd.DataFrame(index=['%s-%s' % (i, j) for i, j in all_perms], columns=['Pathway A', 'Pathway B', '# sig', 
                                                                                                       'sig multilink types',  'actual motif counts',
                                                                                                       'z-scores', 'best z-score', 'p-values',
                                                                                                      'best p-value', 'MuXTalk_score', 'in gold'])

            for (p1, p2) in tqdm(all_perms, position=0, leave=True): 
                temp_df = between_paths_multilink_counts_stats_df_dict[(p1, p2)].loc[between_multilink_pos.keys()]
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'Pathway A'] = p1
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'Pathway B'] = p2
                num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] = len(temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                 (temp_df['p-value']<=0.05)])
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'sig multilink types'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                           (temp_df['p-value']<=0.05)].index.to_numpy()
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'actual motif counts'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                           (temp_df['p-value']<=0.05)]['actual_motif_counts'].to_numpy() 
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'z-scores'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                (temp_df['p-value']<=0.05)]['z-score'].to_numpy()
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                    (temp_df['p-value']<=0.05)]['z-score'].max()
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'p-values'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                    (temp_df['p-value']<=0.05)]['p-value'].to_numpy()
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                    (temp_df['p-value']<=0.05)]['p-value'].min()

                if (num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] > 0) & (~np.isnan(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'])):
                    num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = ((num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] * 1000) +
                                                                                  num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'] *
                                                                                  -np.log10(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] + 0.001))
                elif (num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] > 0) & (np.isnan(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'])):
                    num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = ((num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] * 1000) +
                                                                                  -np.log10(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] + 0.001))          
                else:
                    num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = 0.0

                if '%s-%s' % (p1, p2) in XTalk_DB_positives_common_str:
                    num_sig_multilinks.at['%s-%s' % (p1, p2), 'in gold'] = 'Yes'
                else:
                    num_sig_multilinks.at['%s-%s' % (p1, p2), 'in gold'] = 'No'

            ranked_pathway_pairs_df = num_sig_multilinks.sort_values(['MuXTalk_score'], ascending=False)
            ranked_pathway_pairs = ranked_pathway_pairs_df.index
            detected_ranked_pathway_pairs_df = ranked_pathway_pairs_df[ranked_pathway_pairs_df['# sig'] > 0]
            detected_ranked_pathway_pairs = detected_ranked_pathway_pairs_df.index   
            len_sig = len(ranked_pathway_pairs_df[ranked_pathway_pairs_df['# sig'] > 0])
            
            if parquet==True:
                ranked_pathway_pairs_df.to_parquet(proj_path + multilink_params['GRN'] + '_between_discovery.parquet')
                detected_ranked_pathway_pairs_df.to_parquet(proj_path + multilink_params['GRN'] + '_between_detected_discovery.parquet')
            else:
                ranked_pathway_pairs_df.to_csv(proj_path + multilink_params['GRN'] + '_between_discovery.csv')
                detected_ranked_pathway_pairs_df.to_csv(proj_path + multilink_params['GRN'] + '_between_detected_discovery.csv')                
        
        else:
            
            if parquet==True:
                ranked_pathway_pairs_df = pd.read_parquet(proj_path + multilink_params['GRN'] + '_between_discovery.parquet')
                detected_ranked_pathway_pairs_df = pd.read_parquet(proj_path + multilink_params['GRN'] + '_between_detected_discovery.parquet')
            else:
                ranked_pathway_pairs_df = pd.read_csv(proj_path + multilink_params['GRN'] + '_between_discovery.csv', index_col=0)
                detected_ranked_pathway_pairs_df = pd.read_csv(proj_path + multilink_params['GRN'] + '_between_detected_discovery.csv', index_col=0)  
            
            ranked_pathway_pairs = ranked_pathway_pairs_df.index
            detected_ranked_pathway_pairs = detected_ranked_pathway_pairs_df.index   
            len_sig = len(ranked_pathway_pairs_df[ranked_pathway_pairs_df['# sig'] > 0])
        
    elif method == 'MuXTalk_shortest':
        
        if not os.path.exists(proj_path + multilink_params['GRN'] + '_shortest_sp%s_detected_discovery.parquet' % multilink_params['sp_threshold']):
            
            with open(proj_path + multilink_params['GRN'] + 
                      '_sp%s_shortest_paths_multilink_counts_stats_discovery.pickle' % multilink_params['sp_threshold'], 'rb') as fp:
                shortest_paths_multilink_counts_stats_df_dict = pickle.load(fp)
            (shortest_multilink_pos, shortest_multilink_pos_all) = multilink_positives(shortest_paths_multilink_counts_stats_df_dict)

            all_perms = [(i, j) for i, j in shortest_paths_multilink_counts_stats_df_dict.keys()]
            num_sig_multilinks = pd.DataFrame(index=['%s-%s' % (i, j) for i, j in all_perms], columns=['Pathway A', 'Pathway B', '# sig', 
                                                                                                       'sig multilink types',  'actual motif counts',
                                                                                                       'z-scores', 'best z-score', 'p-values',
                                                                                                      'best p-value', 'MuXTalk_score', 'in gold'])

            for (p1, p2) in tqdm(all_perms, position=0, leave=True): 
                temp_df = shortest_paths_multilink_counts_stats_df_dict[(p1, p2)].loc[shortest_multilink_pos.keys()]
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'Pathway A'] = p1
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'Pathway B'] = p2
                num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] = len(temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                 (temp_df['p-value']<=0.05)])
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'sig multilink types'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                           (temp_df['p-value']<=0.05)].index.to_numpy()
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'actual motif counts'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                           (temp_df['p-value']<=0.05)]['actual_motif_counts'].to_numpy() 
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'z-scores'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                (temp_df['p-value']<=0.05)]['z-score'].to_numpy()
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                    (temp_df['p-value']<=0.05)]['z-score'].max()
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'p-values'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                    (temp_df['p-value']<=0.05)]['p-value'].to_numpy()
                num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] = temp_df[(temp_df['actual_motif_counts'] - temp_df['rand_mean'] > 0) & 
                                                                                    (temp_df['p-value']<=0.05)]['p-value'].min()

                if (num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] > 0) & (~np.isnan(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'])):
                    num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = ((num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] * 1000) +
                                                                                  num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'] *
                                                                                  -np.log10(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] + 0.001))
                elif (num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] > 0) & (np.isnan(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best z-score'])):
                    num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = ((num_sig_multilinks.at['%s-%s' % (p1, p2), '# sig'] * 1000) +
                                                                                  -np.log10(num_sig_multilinks.at['%s-%s' % (p1, p2), 'best p-value'] + 0.001))          
                else:
                    num_sig_multilinks.at['%s-%s' % (p1, p2), 'MuXTalk_score'] = 0.0

                if '%s-%s' % (p1, p2) in XTalk_DB_positives_common_str:
                    num_sig_multilinks.at['%s-%s' % (p1, p2), 'in gold'] = 'Yes'
                else:
                    num_sig_multilinks.at['%s-%s' % (p1, p2), 'in gold'] = 'No'

            ranked_pathway_pairs_df = num_sig_multilinks.sort_values(['MuXTalk_score'], ascending=False) 
            ranked_pathway_pairs = ranked_pathway_pairs_df.index
            detected_ranked_pathway_pairs_df = ranked_pathway_pairs_df[ranked_pathway_pairs_df['# sig'] > 0]
            detected_ranked_pathway_pairs = detected_ranked_pathway_pairs_df.index 
            len_sig = len(ranked_pathway_pairs_df[ranked_pathway_pairs_df['# sig'] > 0])
            
            if parquet==True:
                ranked_pathway_pairs_df.to_parquet(proj_path + multilink_params['GRN'] + '_shortest_sp%s_discovery.parquet' % multilink_params['sp_threshold'])
                detected_ranked_pathway_pairs_df.to_parquet(proj_path + multilink_params['GRN'] + '_shortest_sp%s_detected_discovery.parquet' % multilink_params['sp_threshold'])
            else:
                ranked_pathway_pairs_df.to_csv(proj_path + multilink_params['GRN'] + '_shortest_sp%s_discovery.csv' % multilink_params['sp_threshold'])
                detected_ranked_pathway_pairs_df.to_csv(proj_path + multilink_params['GRN'] + '_shortest_sp%s_detected_discovery.csv' % multilink_params['sp_threshold'])            
            
        else:
            
            if parquet==True:
                ranked_pathway_pairs_df = pd.read_parquet(proj_path + multilink_params['GRN'] + '_shortest_sp%s_discovery.parquet' % multilink_params['sp_threshold'])
                detected_ranked_pathway_pairs_df = pd.read_parquet(proj_path + multilink_params['GRN'] + '_shortest_sp%s_detected_discovery.parquet' % multilink_params['sp_threshold'])
            else:
                ranked_pathway_pairs_df = pd.read_csv(proj_path + multilink_params['GRN'] + '_shortest_sp%s_discovery.csv' % multilink_params['sp_threshold'], index_col=0)
                detected_ranked_pathway_pairs_df = pd.read_csv(proj_path + multilink_params['GRN'] + '_shortest_sp%s_detected_discovery.csv' % multilink_params['sp_threshold'], index_col=0)             
            
            ranked_pathway_pairs = ranked_pathway_pairs_df.index                       
            detected_ranked_pathway_pairs = detected_ranked_pathway_pairs_df.index
            len_sig = len(ranked_pathway_pairs_df[ranked_pathway_pairs_df['# sig'] > 0])
            
    return (detected_ranked_pathway_pairs_df, detected_ranked_pathway_pairs, ranked_pathway_pairs_df, ranked_pathway_pairs, len_sig)      

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

def return_shortest_multilink_edges(p1, p2, multilink_type, KEGG_PPI_allnodes_entrez_df, KEGG_path_nodes_entrez_dict, A_GRN_sparr, 
                                  A_KEGGPPI_sparr, KEGG_interaction_types_dict, A_KEGG_e_sparr_dict, shortest_path_edges_dict):

    multilink_counts = {}   

    A_2_arr = np.asarray(A_GRN_sparr[np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                     np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]

    if '-' in multilink_type:
        sig_ix = int(multilink_type[:len(multilink_type)-2])
        grn_ix = int(multilink_type[-2:])
    else:
        sig_ix = int(multilink_type[:len(multilink_type)-1])
        grn_ix = int(multilink_type[-1:])  

    if multilink_type in ['0-1', '00', '01']:       
        A_1_arr = np.asarray(A_KEGGPPI_sparr[np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                             np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]
        A_mux = np.dstack([A_1_arr, A_2_arr])

    else:
        e = list(KEGG_interaction_types_dict.keys())[list(KEGG_interaction_types_dict.values()).index(sig_ix)]
        A_1_arr = np.asarray(A_KEGG_e_sparr_dict[e][np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 0], 
                                                    np.array(list(shortest_path_edges_dict[(p1, p2)]))[:, 1]])[0]
        A_mux = np.dstack([A_1_arr, A_2_arr])


    sig_edges_ix = np.array(list(shortest_path_edges_dict[(p1, p2)]))[np.argwhere((A_mux[0] == np.array([sig_ix,  grn_ix])).all(-1)).flatten()]    
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

def KEGG_netstats(KEGG_path_nodes_entrez_dict, KEGG_all_edges_entrez):
    
    KEGG_netstats_df = pd.DataFrame(index=KEGG_path_nodes_entrez_dict.keys(), columns=['n_nodes', 'n_edges', 'density'])
    for i in KEGG_path_nodes_entrez_dict.keys():
        KEGG_netstats_df.at[i, 'n_nodes'] = len(KEGG_path_nodes_entrez_dict[i])
        KEGG_netstats_df.at[i, 'n_edges'] = len(KEGG_all_edges_entrez[KEGG_all_edges_entrez['Path_label']==i][['NCBI Gene ID(supplied by NCBI)_x',
                                                                                                            'NCBI Gene ID(supplied by NCBI)_y']]
                                             .drop_duplicates())
        KEGG_netstats_df.at[i, 'density'] = KEGG_netstats_df.at[i, 'n_edges'] / ((KEGG_netstats_df.at[i, 'n_nodes'] * 
                                                                                  (KEGG_netstats_df.at[i, 'n_nodes']-1)))

    return KEGG_netstats_df

def run_MuXTalk(proj_path, input_filenames_dict, input_GRN, MuXTalk_method='MuXTalk_shortest', N_runs=500, N_swap=10, sp_threshold=1, 
                parquet=False):
    
    (KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, PPI_all_edges_entrez, KEGG_PPI_all_edges_entrez, KEGG_PPI_allnodes_entrez_df, 
     KEGG_interaction_types_dict, KEGG_all_edges_entrez, KEGG_all_paths, KEGG_path_nodes_entrez_dict, 
     all_motif_types, all_motif_types_list) = process_data(proj_path, input_GRN, input_filenames_dict)

    (A_GRN_sparr, A_PPI_sparr, A_KEGGPPI_sparr, A_KEGG_e_sparr_dict) = sparse_layers(KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, 
                                                                                     PPI_all_edges_entrez, 
                                                                                     KEGG_all_edges_entrez, KEGG_PPI_all_edges_entrez, 
                                                                                     KEGG_interaction_types_dict)
    
    (common_crosstalk_paths, XTalk_DB_positives_common, XTalk_DB_positives_common_str) = XTalk_DB(KEGG_all_paths, proj_path, input_filenames_dict)
    
    KEGG_all_paths_sansIL17 = set(KEGG_all_paths) - set(['IL-17 signaling pathway'])
    
    A_GRN_sparr_rand_dict = randomize_GRN(A_GRN_sparr, proj_path, input_GRN, N_runs=N_runs, N_swap=N_swap, return_output=True)

    A_KEGG_e_sparr_rand_dict = randomize_KEGG_e(A_KEGG_e_sparr_dict, KEGG_interaction_types_dict, proj_path, input_GRN, N_runs=N_runs, 
                                                N_swap=N_swap, return_output=True)

    A_KEGG_e_sparr_rand_dict_dir = proj_path + 'A_KEGG_e_sparr_rand_dict_500runs.pickle'
    A_KEGGPPI_sparr_rand_dict = randomize_KEGGPPI(KEGG_PPI_allnodes_entrez, KEGG_interaction_types_dict, A_KEGG_e_sparr_rand_dict_dir, proj_path, 
                                                  input_GRN, return_output=True)
    
    if MuXTalk_method == 'MuXTalk_between':
        
        between_paths_multilink_counts_df_dict = between_paths_multilink_counts_discovery(KEGG_all_paths_sansIL17, KEGG_PPI_allnodes_entrez_df, 
                                                                                          KEGG_path_nodes_entrez_dict, A_GRN_sparr, A_KEGGPPI_sparr, 
                                                                                          all_motif_types_list, all_motif_types, 
                                                                                          KEGG_interaction_types_dict, A_KEGG_e_sparr_dict, 
                                                                                          proj_path, input_GRN)

        between_paths_multilink_counts_rand_df_dict = between_paths_multilink_counts_rand_discovery(A_GRN_sparr_rand_dict, KEGG_all_paths_sansIL17,
                                                                                                    KEGG_PPI_allnodes_entrez_df, 
                                                                                                    KEGG_path_nodes_entrez_dict, all_motif_types, 
                                                                                                    A_KEGGPPI_sparr_rand_dict, all_motif_types_list, 
                                                                                                    KEGG_interaction_types_dict, 
                                                                                                    A_KEGG_e_sparr_rand_dict, 
                                                                                                    proj_path, input_GRN, N_rand = 100, N_swap = 10)

        (between_paths_multilink_counts_stats_df_dict, 
         between_paths_rand_zscores_df_dict) = between_paths_multilink_zscores_pvals_discovery(KEGG_all_paths_sansIL17, 
                                                                                               between_paths_multilink_counts_rand_df_dict, 
                                                                                               between_paths_multilink_counts_df_dict, 
                                                                                               all_motif_types, proj_path, input_GRN, N_rand=100)

        (between_detected_ranked_pathway_pairs_df, between_detected_ranked_pathway_pairs, 
         between_ranked_pathway_pairs_df, between_ranked_pathway_pairs, 
         between_len_sig) = get_ranked_pathway_pairs_discovery(XTalk_DB_positives_common_str, proj_path, method=MuXTalk_method, 
                                                               multilink_params={'GRN': input_GRN, 'sp_threshold': sp_threshold}, parquet=parquet)  
        
    
    elif MuXTalk_method == 'MuXTalk_shortest':
        
        shortest_path_edges_dict, shortest_path_intermediaries_dict = get_shortest_paths_discovery(A_KEGGPPI_sparr, KEGG_all_paths_sansIL17, 
                                                                                                   KEGG_PPI_allnodes_entrez_df, 
                                                                                                   KEGG_path_nodes_entrez_dict, proj_path, 
                                                                                                   sp_threshold=sp_threshold)

        shortest_paths_multilink_counts_df_dict = shortest_paths_multilink_counts_discovery(KEGG_all_paths_sansIL17, shortest_path_edges_dict, 
                                                                                            sp_threshold, A_GRN_sparr, A_KEGGPPI_sparr, 
                                                                                            all_motif_types, all_motif_types_list, 
                                                                                            KEGG_interaction_types_dict, A_KEGG_e_sparr_dict, 
                                                                                            proj_path, input_GRN)

        shortest_paths_multilink_counts_rand_df_dict = shortest_paths_multilink_counts_rand_discovery(A_GRN_sparr_rand_dict, KEGG_all_paths_sansIL17, 
                                                                                                      shortest_path_edges_dict, sp_threshold, 
                                                                                                      all_motif_types, A_KEGGPPI_sparr_rand_dict, 
                                                                                                      all_motif_types_list, 
                                                                                                      KEGG_interaction_types_dict, 
                                                                                                      A_KEGG_e_sparr_rand_dict, proj_path, input_GRN, 
                                                                                                      N_rand = 100, N_swap = 10)

        (shortest_paths_multilink_counts_stats_df_dict, 
         shortest_paths_rand_zscores_df_dict) = shortest_paths_multilink_zscores_pvals_discovery(KEGG_all_paths_sansIL17, 
                                                                                                 shortest_paths_multilink_counts_rand_df_dict, 
                                                                                                 sp_threshold, 
                                                                                                 shortest_paths_multilink_counts_df_dict, 
                                                                                                 all_motif_types, proj_path, input_GRN, N_rand=100)

        (shortest_detected_ranked_pathway_pairs_df, shortest_detected_ranked_pathway_pairs, 
         shortest_ranked_pathway_pairs_df, shortest_ranked_pathway_pairs, 
         shortest_len_sig) = get_ranked_pathway_pairs_discovery(XTalk_DB_positives_common_str, proj_path, method=MuXTalk_method, 
                                                                multilink_params={'GRN': input_GRN, 'sp_threshold': sp_threshold}, parquet=parquet)
        

def read_A_KEGG_e_sparr_rand_from_npz(proj_path, npz_folder_name, get_n, get_randomly=True):

    print('Reading ensemble of randomized edge type-specific KEGG sparse matrices from npz files...')
    
    sparr_rand_npz_file_names_dict = {'activation': 'Activation', 'binding/association': 'Binding_Association', 'compound': 'Compound',
                                      'dephosphorylation': 'Dephosphorylation', 'dissociation': 'Dissociation', 
                                      'indirect effect': 'Indirect_effect', 'inhibition': 'Inhibition', 'phosphorylation': 'Phosphorylation', 
                                      'state change': 'State_change', 'ubiquitination': 'Ubiquitination', 'ppi': 'PPI'}

    A_KEGG_e_sparr_rand_dict = {i: {} for i in np.arange(get_n)}
    
    
    if get_randomly:
        rand_ix = random.sample(range(len([name for name in os.listdir(proj_path + npz_folder_name + 'Activation/') if 
                                           name.startswith('A_KEGG_e_')])), get_n)
        for i in tqdm(np.arange(get_n), position=0, leave=True):
            for e in np.concatenate([sorted(list(set(sparr_rand_npz_file_names_dict.keys()) - set(['ppi']))), ['ppi']]):      
                A_KEGG_e_sparr_rand_dict[i][e] = sparse.load_npz(proj_path + npz_folder_name + '%s/A_KEGG_e_sparr_rand_%s_%s.npz' % 
                                                                 (sparr_rand_npz_file_names_dict[e], rand_ix[i], sparr_rand_npz_file_names_dict[e]))
                
    else:
        for i in tqdm(np.arange(get_n), position=0, leave=True):
            for e in np.concatenate([sorted(list(set(sparr_rand_npz_file_names_dict.keys()) - set(['ppi']))), ['ppi']]):      
                A_KEGG_e_sparr_rand_dict[i][e] = sparse.load_npz(proj_path + npz_folder_name + '%s/A_KEGG_e_sparr_rand_%s_%s.npz' % 
                                                                 (sparr_rand_npz_file_names_dict[e], i, sparr_rand_npz_file_names_dict[e]))   
                
    return A_KEGG_e_sparr_rand_dict

def read_A_GRN_sparr_rand_from_npz(proj_path, input_GRN, npz_folder_name, get_n, get_randomly=True):

    print('Reading ensemble of randomized GRN sparse matrices from npz files...')
    
    A_GRN_sparr_rand_dict = {}
    
    if get_randomly:
        rand_ix = random.sample(range(len([name for name in os.listdir(proj_path + input_GRN + '_' + npz_folder_name) if 
                                           name.startswith('%s_A_GRN_' % input_GRN)])), get_n)    
        for i in tqdm(np.arange(get_n), position=0, leave=True):
            A_GRN_sparr_rand_dict[i] = sparse.load_npz(proj_path + input_GRN + '_' + npz_folder_name + '%s_A_GRN_sparr_rand_%s.npz' % 
                                                       (input_GRN, rand_ix[i]))
    else:
        for i in tqdm(np.arange(get_n), position=0, leave=True):
            A_GRN_sparr_rand_dict[i] = sparse.load_npz(proj_path + input_GRN + '_' + npz_folder_name + '%s_A_GRN_sparr_rand_%s.npz' % 
                                                       (input_GRN, i))       
            
    return A_GRN_sparr_rand_dict

def read_A_KEGGPPI_sparr_rand_from_npz(proj_path, npz_folder_name, get_n, get_randomly=True):
    
    print('Reading ensemble of randomized KEGG+PPI combined sparse matrices from npz files...')
    
    A_KEGGPPI_sparr_rand_dict = {}
    
    if get_randomly:
        rand_ix = random.sample(range(len([name for name in os.listdir(proj_path + npz_folder_name) if name.startswith('A_KEGGPPI_')])), get_n) 
        for i in tqdm(np.arange(get_n), position=0, leave=True):
            A_KEGGPPI_sparr_rand_dict[i] = sparse.load_npz(proj_path + npz_folder_name + 'A_KEGGPPI_sparr_rand_%s.npz' % rand_ix[i])
    else:
        for i in tqdm(np.arange(get_n), position=0, leave=True):
            A_KEGGPPI_sparr_rand_dict[i] = sparse.load_npz(proj_path + npz_folder_name + 'A_KEGGPPI_sparr_rand_%s.npz' % i)
            
    return A_KEGGPPI_sparr_rand_dict

def read_A_KEGG_sparr_rand_from_npz(proj_path, npz_folder_name, get_n, get_randomly=True):
    
    print('Reading ensemble of randomized KEGG sparse matrices from npz files...')
    
    A_KEGG_sparr_rand_dict = {}
    
    if get_randomly:
        rand_ix = random.sample(range(len([name for name in os.listdir(proj_path + npz_folder_name) if name.startswith('A_KEGG_')])), get_n)    
        for i in tqdm(np.arange(get_n), position=0, leave=True):
            A_KEGG_sparr_rand_dict[i] = sparse.load_npz(proj_path + npz_folder_name + 'A_KEGG_sparr_rand_%s.npz' % rand_ix[i])
    else:
        for i in tqdm(np.arange(get_n), position=0, leave=True):
            A_KEGG_sparr_rand_dict[i] = sparse.load_npz(proj_path + npz_folder_name + 'A_KEGG_sparr_rand_%s.npz' % i)
            
    return A_KEGG_sparr_rand_dict

def randomize_GRN_npz(A_GRN_sparr, proj_path, input_GRN, npz_folder_name, N_runs=500, N_swap=10, get_n=100, get_randomly=True, return_output=False):
    
    if not os.path.exists(proj_path + input_GRN + '_' + npz_folder_name):
        
        os.makedirs(proj_path + input_GRN + '_' + npz_folder_name)
        
        print('Computing ensemble of randomized GRN sparse  matrices...')
        A_GRN_sparr_rand_dict = {}
        for nrand in tqdm(np.arange(N_runs), position=0, leave=True):
            A_GRN_arr_rand, A_GRN_nrew_rand = edge_swap_undir_sign(A_GRN_sparr.toarray(), N_swap)
            A_GRN_sparr_rand_dict[nrand] = sparse.csr_matrix(A_GRN_arr_rand)
            sparse.save_npz(proj_path + input_GRN + '_' + npz_folder_name + '%s_A_GRN_sparr_rand_%s.npz' % (input_GRN, nrand), 
                            A_GRN_sparr_rand_dict[nrand])          
        
        if return_output==True:        
            return A_GRN_sparr_rand_dict
            
    elif len(os.listdir(proj_path + input_GRN + '_' + npz_folder_name)) == 0:
        
        print('Computing ensemble of randomized GRN sparse  matrices...')
        A_GRN_sparr_rand_dict = {}
        for nrand in tqdm(np.arange(N_runs), position=0, leave=True):
            A_GRN_arr_rand, A_GRN_nrew_rand = edge_swap_undir_sign(A_GRN_sparr.toarray(), N_swap)
            A_GRN_sparr_rand_dict[nrand] = sparse.csr_matrix(A_GRN_arr_rand)
            sparse.save_npz(proj_path + input_GRN + '_' + npz_folder_name + '%s_A_GRN_sparr_rand_%s.npz' % (input_GRN, nrand), 
                            A_GRN_sparr_rand_dict[nrand])          
        
        if return_output==True:        
            return A_GRN_sparr_rand_dict    
            
    else:
        
        A_GRN_sparr_rand_dict = read_A_GRN_sparr_rand_from_npz(proj_path, input_GRN, npz_folder_name, get_n=get_n, get_randomly=get_randomly)
                
        return A_GRN_sparr_rand_dict

def randomize_KEGG_e_npz(A_KEGG_e_sparr_dict, KEGG_interaction_types_dict, proj_path, input_GRN, npz_folder_name, N_runs=500, N_swap=10, 
                         get_n=100, get_randomly=True, return_output=False):

    if not os.path.exists(proj_path + npz_folder_name):   
 
        sparr_rand_npz_file_names_dict = {'activation': 'Activation', 'binding/association': 'Binding_Association', 'compound': 'Compound',
                              'dephosphorylation': 'Dephosphorylation', 'dissociation': 'Dissociation', 'indirect effect': 'Indirect_effect',
                              'inhibition': 'Inhibition', 'phosphorylation': 'Phosphorylation', 'state change': 'State_change', 
                              'ubiquitination': 'Ubiquitination', 'ppi': 'PPI'}

        os.makedirs(proj_path + npz_folder_name)
        for e in sparr_rand_npz_file_names_dict.keys():
            os.makedirs(proj_path + npz_folder_name + sparr_rand_npz_file_names_dict[e])
           
        print('Computing ensemble of randomized edge type-specific KEGG sparse matrices...')
        A_KEGG_e_sparr_rand_dict = {i: {} for i in np.arange(N_runs)}

        for nrand in tqdm(np.arange(N_runs), position=0, leave=True):
            for e in np.concatenate([sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))), ['ppi']]):
                A_KEGG_e_arr_rand, A_KEGG_e_nrew_rand = edge_swap_undir_sign(A_KEGG_e_sparr_dict[e].toarray(), N_swap)
                A_KEGG_e_sparr_rand_dict[nrand][e] = sparse.csr_matrix(A_KEGG_e_arr_rand) 
                sparse.save_npz(proj_path + npz_folder_name + '%s/A_KEGG_e_sparr_rand_%s_%s.npz' % 
                                (sparr_rand_npz_file_names_dict[e], nrand, sparr_rand_npz_file_names_dict[e]), A_KEGG_e_sparr_rand_dict[nrand][e])

        if return_output==True:
            return A_KEGG_e_sparr_rand_dict
        
    else:
        
        A_KEGG_e_sparr_rand_dict = read_A_KEGG_e_sparr_rand_from_npz(proj_path, npz_folder_name, get_n=get_n, get_randomly=get_randomly)
         
        return A_KEGG_e_sparr_rand_dict        

def randomize_KEGGPPI_npz(KEGG_PPI_allnodes_entrez, KEGG_interaction_types_dict, proj_path, input_GRN, 
                      A_KEGG_e_npz_folder_name, A_KEGGPPI_npz_folder_name, get_n=100, get_randomly=True, return_output=False):
    
    if not os.path.exists(proj_path + A_KEGGPPI_npz_folder_name):  
        
        os.makedirs(proj_path + A_KEGGPPI_npz_folder_name)
        
        A_KEGG_e_sparr_rand_dict = read_A_KEGG_e_sparr_rand_from_npz(proj_path, A_KEGG_e_npz_folder_name, get_n=get_n, get_randomly=get_randomly)
    
        print('Computing ensemble of randomized KEGG+PPI combined sparse  matrices...')
        A_KEGGPPI_sparr_rand_dict = {}

        for nrand in tqdm(np.arange(len(A_KEGG_e_sparr_rand_dict)), position=0, leave=True):

            temp_sparse = sparse.csr_matrix((len(KEGG_PPI_allnodes_entrez), len(KEGG_PPI_allnodes_entrez)))

            # get only the KEGG interaction types from the A_KEGG_e adjacency matrices (which include both the KEGG interaction type and ppi (9))
            for e in sorted(list(set(KEGG_interaction_types_dict.keys()) - set(['ppi']))):
                temp_sparse = temp_sparse + (abs(A_KEGG_e_sparr_rand_dict[nrand][e])==KEGG_interaction_types_dict[e]).astype(int)

            # add at the end the ppi type
            temp_sparse = temp_sparse + (abs(A_KEGG_e_sparr_rand_dict[nrand]['ppi'])==KEGG_interaction_types_dict['ppi']).astype(int)      

            A_KEGGPPI_sparr_rand_dict[nrand] = (temp_sparse!=0).astype(int)
            sparse.save_npz(proj_path + A_KEGGPPI_npz_folder_name + 'A_KEGGPPI_sparr_rand_%s.npz' % nrand, A_KEGGPPI_sparr_rand_dict[nrand])            
      
        if return_output==True:
            return A_KEGGPPI_sparr_rand_dict
        
    else:
        
        A_KEGGPPI_sparr_rand_dict = read_A_KEGGPPI_sparr_rand_from_npz(proj_path, A_KEGGPPI_npz_folder_name, get_n=get_n, 
                                                                       get_randomly=get_randomly)

        return A_KEGGPPI_sparr_rand_dict     

def randomize_KEGG_npz(A_KEGG_sparr, proj_path, npz_folder_name, N_runs=500, N_swap=10, get_n=100, get_randomly=True, return_output=False):

    if not os.path.exists(proj_path + npz_folder_name):    

        os.makedirs(proj_path + npz_folder_name)
        
        print('Computing ensemble of randomized KEGG sparse matrices...')    
    
        A_KEGG_sparr_rand_dict = {}

        for nrand in tqdm(np.arange(N_runs), position=0, leave=True):
            A_KEGG_arr_rand, A_KEGG_nrew_rand = edge_swap(A_KEGG_sparr.toarray(), N_swap)
            A_KEGG_sparr_rand_dict[nrand] = sparse.csr_matrix(A_KEGG_arr_rand)
            sparse.save_npz(proj_path + npz_folder_name + 'A_KEGG_sparr_rand_%s.npz' % nrand, A_KEGG_sparr_rand_dict[nrand])            
                
        if return_output==True:        
            return A_KEGG_sparr_rand_dict
        
    else:
        
        A_KEGG_sparr_rand_dict = read_A_KEGG_sparr_rand_from_npz(proj_path, npz_folder_name, get_n=get_n, get_randomly=get_randomly)
    
        return A_KEGG_sparr_rand_dict        

def run_MuXTalk_npz(proj_path, input_filenames_dict, npz_filenames_dict, input_GRN, MuXTalk_method='MuXTalk_shortest', 
                    get_n=100, get_randomly=True, N_runs=500, N_swap=10, sp_threshold=1, parquet=False):
    
    sparr_rand_npz_file_names_dict = {'activation': 'Activation', 'binding/association': 'Binding_Association', 'compound': 'Compound',
                                  'dephosphorylation': 'Dephosphorylation', 'dissociation': 'Dissociation', 'indirect effect': 'Indirect_effect',
                                  'inhibition': 'Inhibition', 'phosphorylation': 'Phosphorylation', 'state change': 'State_change', 
                                  'ubiquitination': 'Ubiquitination', 'ppi': 'PPI'}
    
    (KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, PPI_all_edges_entrez, KEGG_PPI_all_edges_entrez, KEGG_PPI_allnodes_entrez_df, 
     KEGG_interaction_types_dict, KEGG_all_edges_entrez, KEGG_all_paths, KEGG_path_nodes_entrez_dict, 
     all_motif_types, all_motif_types_list) = process_data(proj_path, input_GRN, input_filenames_dict)

    (A_GRN_sparr, A_PPI_sparr, A_KEGGPPI_sparr, A_KEGG_e_sparr_dict) = sparse_layers(KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, 
                                                                                     PPI_all_edges_entrez, 
                                                                                     KEGG_all_edges_entrez, KEGG_PPI_all_edges_entrez, 
                                                                                     KEGG_interaction_types_dict)
    
    (common_crosstalk_paths, XTalk_DB_positives_common, XTalk_DB_positives_common_str) = XTalk_DB(KEGG_all_paths, proj_path, input_filenames_dict)
    
    KEGG_all_paths_sansIL17 = set(KEGG_all_paths) - set(['IL-17 signaling pathway'])
    
    A_GRN_sparr_rand_dict = randomize_GRN_npz(A_GRN_sparr, proj_path, input_GRN, npz_filenames_dict['GRN'], N_runs=500, N_swap=10, 
                                              get_n=get_n, get_randomly=get_randomly, return_output=False)
    A_KEGG_e_sparr_rand_dict = randomize_KEGG_e_npz(A_KEGG_e_sparr_dict, KEGG_interaction_types_dict, proj_path, input_GRN, 
                                                     npz_filenames_dict['KEGG_e'], N_runs=500, N_swap=10, get_n=get_n, get_randomly=get_randomly, 
                                                    return_output=False)
    A_KEGGPPI_sparr_rand_dict = randomize_KEGGPPI_npz(KEGG_PPI_allnodes_entrez, KEGG_interaction_types_dict, proj_path, 
                                                  input_GRN, npz_filenames_dict['KEGG_e'], npz_filenames_dict['KEGGPPI'], 
                                                  get_n=get_n, get_randomly=get_randomly, return_output=False)
    
    if MuXTalk_method == 'MuXTalk_between':
        
        between_paths_multilink_counts_df_dict = between_paths_multilink_counts_discovery(KEGG_all_paths_sansIL17, KEGG_PPI_allnodes_entrez_df, 
                                                                                          KEGG_path_nodes_entrez_dict, A_GRN_sparr, A_KEGGPPI_sparr, 
                                                                                          all_motif_types_list, all_motif_types, 
                                                                                          KEGG_interaction_types_dict, A_KEGG_e_sparr_dict, 
                                                                                          proj_path, input_GRN)

        between_paths_multilink_counts_rand_df_dict = between_paths_multilink_counts_rand_discovery(A_GRN_sparr_rand_dict, KEGG_all_paths_sansIL17,
                                                                                                    KEGG_PPI_allnodes_entrez_df, 
                                                                                                    KEGG_path_nodes_entrez_dict, all_motif_types, 
                                                                                                    A_KEGGPPI_sparr_rand_dict, all_motif_types_list, 
                                                                                                    KEGG_interaction_types_dict, 
                                                                                                    A_KEGG_e_sparr_rand_dict, 
                                                                                                    proj_path, input_GRN, N_rand = 100, N_swap = 10)

        (between_paths_multilink_counts_stats_df_dict, 
         between_paths_rand_zscores_df_dict) = between_paths_multilink_zscores_pvals_discovery(KEGG_all_paths_sansIL17, 
                                                                                               between_paths_multilink_counts_rand_df_dict, 
                                                                                               between_paths_multilink_counts_df_dict, 
                                                                                               all_motif_types, proj_path, input_GRN, N_rand=100)

        (between_detected_ranked_pathway_pairs_df, between_detected_ranked_pathway_pairs, 
         between_ranked_pathway_pairs_df, between_ranked_pathway_pairs, 
         between_len_sig) = get_ranked_pathway_pairs_discovery(XTalk_DB_positives_common_str, proj_path, method=MuXTalk_method, 
                                                               multilink_params={'GRN': input_GRN, 'sp_threshold': sp_threshold}, parquet=parquet)  
        
    
    elif MuXTalk_method == 'MuXTalk_shortest':
        
        shortest_path_edges_dict, shortest_path_intermediaries_dict = get_shortest_paths_discovery(A_KEGGPPI_sparr, KEGG_all_paths_sansIL17, 
                                                                                                   KEGG_PPI_allnodes_entrez_df, 
                                                                                                   KEGG_path_nodes_entrez_dict, proj_path, 
                                                                                                   sp_threshold=sp_threshold)

        shortest_paths_multilink_counts_df_dict = shortest_paths_multilink_counts_discovery(KEGG_all_paths_sansIL17, shortest_path_edges_dict, 
                                                                                            sp_threshold, A_GRN_sparr, A_KEGGPPI_sparr, 
                                                                                            all_motif_types, all_motif_types_list, 
                                                                                            KEGG_interaction_types_dict, A_KEGG_e_sparr_dict, 
                                                                                            proj_path, input_GRN)

        shortest_paths_multilink_counts_rand_df_dict = shortest_paths_multilink_counts_rand_discovery(A_GRN_sparr_rand_dict, KEGG_all_paths_sansIL17, 
                                                                                                      shortest_path_edges_dict, sp_threshold, 
                                                                                                      all_motif_types, A_KEGGPPI_sparr_rand_dict, 
                                                                                                      all_motif_types_list, 
                                                                                                      KEGG_interaction_types_dict, 
                                                                                                      A_KEGG_e_sparr_rand_dict, proj_path, input_GRN, 
                                                                                                      N_rand = 100, N_swap = 10)

        (shortest_paths_multilink_counts_stats_df_dict, 
         shortest_paths_rand_zscores_df_dict) = shortest_paths_multilink_zscores_pvals_discovery(KEGG_all_paths_sansIL17, 
                                                                                                 shortest_paths_multilink_counts_rand_df_dict, 
                                                                                                 sp_threshold, 
                                                                                                 shortest_paths_multilink_counts_df_dict, 
                                                                                                 all_motif_types, proj_path, input_GRN, N_rand=100)

        (shortest_detected_ranked_pathway_pairs_df, shortest_detected_ranked_pathway_pairs, 
         shortest_ranked_pathway_pairs_df, shortest_ranked_pathway_pairs, 
         shortest_len_sig) = get_ranked_pathway_pairs_discovery(XTalk_DB_positives_common_str, proj_path, method=MuXTalk_method, 
                                                                multilink_params={'GRN': input_GRN, 'sp_threshold': sp_threshold}, parquet=parquet)
        
        