import pandas as pd
import pickle
from tqdm import tqdm
import numpy as np
import os
import networkx as nx
import scipy.stats as st
import statsmodels.api as sm


def input_paths(proj_path):

    input_paths_dict = {'QUIZC_cids_inchi_smiles_path': proj_path + 'Pathopticon_intermediary_outputs/' + 'QUIZC_drug_CIDs_inchi_smiles.csv',
                        'tool_path': proj_path + 'Pathopticon_intermediary_outputs/' + 'QUIZC_activityStats_nooutliers_df_besttool.csv',
                        'benchmark_genesets_path': proj_path + 'Pathopticon_intermediary_outputs/' + 'benchmark_genesets.csv',
                        
                        'L1000_gene_info_path': proj_path + 'Pathopticon_external_data/L1000/' + 'GSE92742_Broad_LINCS_gene_info.txt',
                        'L1000_cell_info_path': proj_path + 'Pathopticon_external_data/L1000/' + 'GSE92742_Broad_LINCS_cell_info.txt',
                        'L1000_inst_info_path': proj_path + 'Pathopticon_external_data/L1000/' + 'GSE92742_Broad_LINCS_inst_info.txt',
                        'L1000_pert_info_path': proj_path + 'Pathopticon_external_data/L1000/' + 'GSE92742_Broad_LINCS_pert_info.txt',
                        'L1000_sig_info_path': proj_path + 'Pathopticon_external_data/L1000/' + 'GSE92742_Broad_LINCS_sig_info.txt',

                        'pos_edges_dict_path': proj_path + 'QUIZC75/' + 'cell_pos_edges_dict_75.pickle',
                        'neg_edges_dict_path': proj_path + 'QUIZC75/' + 'cell_neg_edges_dict_75.pickle',
                        'drugs_dict_path': proj_path + 'QUIZC75/' + 'cell_drugs_dict_75.pickle',

                        'cgp_dir': proj_path + 'Pathopticon_external_data/MSigDB/' + 'c2.cgp.v7.1.symbols.gmt',

                        'Enrichr_GEO_up_path': proj_path + 'Pathopticon_external_data/Enrichr/' + 'Disease_Perturbations_from_GEO_up.txt',
                        'Enrichr_GEO_dn_path': proj_path + 'Pathopticon_external_data/Enrichr/' + 'Disease_Perturbations_from_GEO_down.txt',

                        'TTD_drugs_path': proj_path + 'Pathopticon_external_data/TTD/' + 'P1-02-TTD_drug_download.txt',
                        'TTD_InChI2CID_path': proj_path + 'Pathopticon_external_data/TTD/' + 'TTD_drugs_InChI2CID.txt',
                        'TTD_drug_target_path': proj_path + 'Pathopticon_external_data/TTD/' + 'P1-07-Drug-TargetMapping.csv',
                        'TTD_target_path': proj_path + 'Pathopticon_external_data/TTD/' + 'P1-01-TTD_target_download.txt',
                        'MODZ_networks_path': proj_path + 'MODZ_networks/',
                        'CD_networks_path': proj_path + 'CD_networks/',
                        
                        'benchmark_path': proj_path + 'Pathopticon_benchmark_outputs/'}
    
    return input_paths_dict

def import_L1000_metadata(L1000_gene_info_path, L1000_cell_info_path, L1000_inst_info_path, L1000_pert_info_path, L1000_sig_info_path):

    L1000_gene_info = pd.read_csv(L1000_gene_info_path, sep='\t', low_memory=False)
    L1000_cell_info = pd.read_csv(L1000_cell_info_path, sep='\t', low_memory=False)
    L1000_inst_info = pd.read_csv(L1000_inst_info_path, sep='\t', low_memory=False)
    L1000_pert_info = pd.read_csv(L1000_pert_info_path, sep='\t', low_memory=False)
    L1000_sig_info = pd.read_csv(L1000_sig_info_path, sep='\t', low_memory=False)
    
    return L1000_gene_info, L1000_cell_info, L1000_inst_info, L1000_pert_info, L1000_sig_info

def process_QUIZC_output(pos_edges_dict_path, neg_edges_dict_path, drugs_dict_path, L1000_gene_info):

    with open(pos_edges_dict_path, 'rb') as handle:
        cell_pos_edges_dict = pickle.load(handle)
    with open(neg_edges_dict_path, 'rb') as handle:
        cell_neg_edges_dict = pickle.load(handle)
    with open(drugs_dict_path, 'rb') as handle:
        cell_drugs_dict = pickle.load(handle)     
    
    edgelist_df_dict = {}
    nodelist_df_dict = {}
    for c in tqdm(cell_drugs_dict.keys(), position=0, leave=True):

        pos_edges_df = pd.DataFrame(np.array([(i, j) for i, j in zip(cell_drugs_dict[c], cell_pos_edges_dict[c])])).explode(1)
        neg_edges_df = pd.DataFrame(np.array([(i, j) for i, j in zip(cell_drugs_dict[c], cell_neg_edges_dict[c])])).explode(1)

        pos_edges_df = pos_edges_df[~pd.isnull(pos_edges_df[1])]
        neg_edges_df = neg_edges_df[~pd.isnull(neg_edges_df[1])]

        pos_edges_df = pos_edges_df.reset_index(drop=True).rename(columns={0: 'Drug', 1: 'Target'})
        neg_edges_df = neg_edges_df.reset_index(drop=True).rename(columns={0: 'Drug', 1: 'Target'})

        pos_edges_df['Target'] = pos_edges_df['Target'].astype(int)
        neg_edges_df['Target'] = neg_edges_df['Target'].astype(int)

        pos_edges_df['Direction'] = 'Up'
        neg_edges_df['Direction'] = 'Down'

        pos_edges_df = pd.merge(pos_edges_df, L1000_gene_info, left_on='Target', 
                                right_on='pr_gene_id', how='left')[['Drug', 'pr_gene_symbol', 'Direction']]
        neg_edges_df = pd.merge(neg_edges_df, L1000_gene_info, left_on='Target', 
                                right_on='pr_gene_id', how='left')[['Drug', 'pr_gene_symbol', 'Direction']]

        all_edges_df = pd.concat([pos_edges_df, neg_edges_df]).rename(columns={'pr_gene_symbol': 'Target'})

        if len(all_edges_df) > 0:
            edgelist_df_dict[c] = all_edges_df 

            nodelist_df_dict_drug = pd.DataFrame(columns=['Id', 'Label', 'Type'])
            nodelist_df_dict_drug['Id'] = edgelist_df_dict[c]['Drug'].unique()
            nodelist_df_dict_drug['Label'] = edgelist_df_dict[c]['Drug'].unique()
            nodelist_df_dict_drug['Type'] = 'Drug'
            nodelist_df_dict_target = pd.DataFrame(columns=['Id', 'Label', 'Type'])
            nodelist_df_dict_target['Id'] = edgelist_df_dict[c]['Target'].unique()
            nodelist_df_dict_target['Label'] = edgelist_df_dict[c]['Target'].unique()
            nodelist_df_dict_target['Type'] = 'Gene'

            nodelist_df_dict[c] = pd.concat([nodelist_df_dict_drug, nodelist_df_dict_target])


    allcells = pd.DataFrame(sorted(list(nodelist_df_dict.keys())), columns=['Cell_type'])
    allgenes = pd.DataFrame(sorted(list(set(pd.concat(nodelist_df_dict)[pd.concat(nodelist_df_dict)['Type']=='Gene']['Id']))), columns=['Gene_symbol'])
    alldrugs = pd.DataFrame(sorted(list(set(pd.concat(nodelist_df_dict)[pd.concat(nodelist_df_dict)['Type']=='Drug']['Id']))), columns=['Pert_iname'])

    allnodes = pd.DataFrame(columns=['Node_name'])
    allnodes['Node_name'] = np.concatenate([allgenes['Gene_symbol'].values, alldrugs['Pert_iname'].values])
    allnodes['Node_ID'] = allnodes.index.values
    
    return edgelist_df_dict, nodelist_df_dict, allcells, allgenes, alldrugs, allnodes

def import_CD_MODZ_networks(path):
    nodelist_df_dict = {}
    edgelist_df_dict = {}
    for f in os.listdir(path):
        if f.endswith('_nodes.csv'):
            nodelist_df_dict[f.split('_')[1]] = pd.read_csv(path + f)
            nodelist_df_dict[f.split('_')[1]] = nodelist_df_dict[f.split('_')[1]].rename(columns={'ID': 'Id'})
        elif f.endswith('_edges_dir.csv'):
            edgelist_df_dict[f.split('_')[1]] = pd.read_csv(path + f)
            edgelist_df_dict[f.split('_')[1]] = edgelist_df_dict[f.split('_')[1]].rename(columns={'pert_iname': 'Drug', 'gene_name': 'Target'})[['Drug', 'Target', 'Direction']]

    allcells = pd.DataFrame(sorted(list(nodelist_df_dict.keys())), columns=['Cell_type'])
    allgenes = pd.DataFrame(sorted(list(set(pd.concat(nodelist_df_dict)[pd.concat(nodelist_df_dict)['Type']=='Gene']['Id']))), 
                                columns=['Gene_symbol'])
    allgenes = allgenes.drop(0) # remove the row with gene name "-666", this is NA in Broad Institute convention
    alldrugs = pd.DataFrame(sorted(list(set(pd.concat(nodelist_df_dict)[pd.concat(nodelist_df_dict)['Type']=='Drug']['Id']))), 
                                columns=['Pert_iname'])

    allnodes = pd.DataFrame(columns=['Node_name'])
    allnodes['Node_name'] = np.concatenate([allgenes['Gene_symbol'].values, alldrugs['Pert_iname'].values])
    allnodes['Node_ID'] = allnodes.index.values
    
    return edgelist_df_dict, nodelist_df_dict, allcells, allgenes, alldrugs, allnodes

def generate_diG_dict(nodelist_df_dict, edgelist_df_dict, allcells):

    diG_dict = {}
    for c in tqdm(allcells['Cell_type'], position=0, leave=True):
        diG_dict[c] = nx.DiGraph()
        for ix, node in enumerate(nodelist_df_dict[c]['Id'].values):
            if nodelist_df_dict[c].iloc[ix]['Type']=='Drug':
                diG_dict[c].add_node(node, node_type='Drug', node_color='limegreen', node_size=150)
            elif nodelist_df_dict[c].iloc[ix]['Type']=='Gene':
                diG_dict[c].add_node(node, node_type='Gene', node_color='darkslateblue', node_size=50) 

        for edge in edgelist_df_dict[c].values:       
            if edge[2]=='Up':
                diG_dict[c].add_edge(edge[0], edge[1], edge_type='Up', edge_color='red')
            elif edge[2]=='Down':
                diG_dict[c].add_edge(edge[0], edge[1], edge_type='Down', edge_color='deepskyblue')
                
    return diG_dict

def process_MSigDB_CGP(cgp_dir):
    
    cgp = {}
    f = open(cgp_dir)
    for line in f:
        line = line.rstrip()
        cgp[line.split('http://')[0].split('\t')[0]] = line.split('http://')[1].split('\t')[1:]
    f.close()
        
    # identify the indices of the duplicates (after removing _UP and DN) to determine the genesets that have both _UP and _DN
    cgp_ix = pd.DataFrame(index=cgp.keys())
    cgp_ix = cgp_ix[(cgp_ix.index.str.endswith('_UP')) | (cgp_ix.index.str.endswith('_DN'))]
    cgp_ix = cgp_ix[cgp_ix.index.str[:-3].duplicated(keep=False)]

    # CGP dictionary with only the genesets with UP and DN
    cgp_updn = {i: cgp[i] for i in cgp_ix.index.values}
        
    cgp_updn_allgenes = set([j for i in cgp_updn.values() for j in i])
    cgp_updn_labels = sorted(list(set(cgp_ix.index.str[:-3])))
   
    return cgp_updn, cgp_updn_labels, cgp_updn_allgenes

def import_TDD(TTD_drugs_path, TTD_InChI2CID_path, TTD_drug_target_path, TTD_target_path):
    
    TTD_drugs = pd.read_csv(TTD_drugs_path, sep='\t', skiprows=29)
    TTD_InChI2CID = pd.read_csv(TTD_InChI2CID_path, sep='\t', header=None)
    TTD_drugs_CID = pd.merge(TTD_drugs[TTD_drugs['DRUG__ID']=='DRUGINKE'], TTD_InChI2CID, left_on='D00AAN.1', right_on=0)
    TTD_drug_target = pd.read_csv(TTD_drug_target_path)

    f = open(TTD_target_path, 'r')
    TTD_targID_dict = {}
    for line in f:
        line = line.split('\t')
        if len(line) == 3:
            if line[1] == 'GENENAME':
                TTD_targID_dict[line[2].strip('\n')] = line[0]
    f.close()
    
    return TTD_drugs, TTD_InChI2CID, TTD_drugs_CID, TTD_drug_target, TTD_targID_dict

def get_geneset_targets(geneset_up, geneset_dn, TTD_targID_dict, TTD_drug_target, TTD_drugs_CID, pert_iname2CID):

    # getting the direct targets of the input geneset from TTD
    geneset_targetIDs = [TTD_targID_dict[i] for i in list(set(TTD_targID_dict.keys()) & (geneset_dn | geneset_up))]
    geneset_drugID_targetID = TTD_drug_target[TTD_drug_target['TargetID'].isin(geneset_targetIDs)]
    geneset_drugCIDs = TTD_drugs_CID[TTD_drugs_CID['D00AAN'].isin(geneset_drugID_targetID['DrugID'])]
    geneset_drugID_targetCID = pd.merge(geneset_drugID_targetID, geneset_drugCIDs, left_on='DrugID', right_on='D00AAN', how='left')


    # Get the subset of drugs targeting the given pathway above that is within the QUIZC drugs by
    # (1) first getting the lookup table between CIDs, InCHI Keys, SMILES and pert_inames of QUIZC compounds;
    # (2) then subsetting the drugs targeting the given pathway by the ones that are in QUIZC, in terms of their pert_inames.
    # The resulting number of compounds is naturally small since the drugs in QUIZC networks are a small subset of all the drugs 
    # in the TTD database.

    geneset_targets_pert_iname = pd.merge(geneset_drugID_targetCID, pert_iname2CID, 
                                                    left_on=1, right_on='pubchem_cid_x', how='left')
    geneset_targets_pert_iname = geneset_targets_pert_iname[~pd.isnull(geneset_targets_pert_iname['Pert_iname'])]

    return set(geneset_targets_pert_iname['Pert_iname'])

def get_MSigDB_geneset_targets(cgp_updn_labels, cgp_updn, TTD_targID_dict, TTD_drug_target, TTD_drugs_CID, pert_iname2CID):
    
    geneset_pert_iname_dict = {}
    for geneset_name in tqdm(cgp_updn_labels, position=0, leave=True):

        geneset_up = set(cgp_updn[geneset_name + '_UP'])
        geneset_dn = set(cgp_updn[geneset_name + '_DN'])        
        geneset_pert_iname_dict[geneset_name] = get_geneset_targets(geneset_up, geneset_dn, TTD_targID_dict, TTD_drug_target, 
                                                                         TTD_drugs_CID, pert_iname2CID)
        
    return geneset_pert_iname_dict

def import_Enrichr_GEO(Enrichr_GEO_up_path, Enrichr_GEO_dn_path):

    Enrichr_GEO_disease_human_up = {}
    f = open(Enrichr_GEO_up_path)
    for line in f:
        line = line.rstrip()
        d = line.split('\t')[0]
        g = line.split('\t')[2:]
        if 'human' in d:
            Enrichr_GEO_disease_human_up[d] = g       
    f.close()

    Enrichr_GEO_disease_human_dn = {}
    f = open(Enrichr_GEO_dn_path)
    for line in f:
        line = line.rstrip()
        d = line.split('\t')[0]
        g = line.split('\t')[2:]
        if 'human' in d:
            Enrichr_GEO_disease_human_dn[d] = g       
    f.close()
    
    return Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn

# Calculate the Signature Congruity Score (SCS) between two genesets, typically an input signature (or perturbation) and a disease signature.
def SCS(geneset_up, geneset_dn, disease_up, disease_dn):
        
        if len((geneset_up | geneset_dn) & (disease_up | disease_dn)) > 0:
            SCS_g_d = 1.0 * ((len(disease_up & geneset_up) + len(disease_dn & geneset_dn) 
                       - len(disease_dn & geneset_up) - len(disease_up & geneset_dn)) / len(disease_up | disease_dn))
            return SCS_g_d
        else:
            return np.nan        

def PCP_geneset(geneset_up, geneset_dn, disease_sig_up_dict, disease_sig_dn_dict, geneset_name='Geneset'):

    PCP_geneset_df = pd.DataFrame(index=[geneset_name], columns=sorted(disease_sig_up_dict.keys()), dtype='float')
    
    for d in disease_sig_up_dict.keys():
        disease_up = set(disease_sig_up_dict[d])
        disease_dn = set(disease_sig_dn_dict[d])    
        
        PCP_geneset_df.at[geneset_name, d] = SCS(geneset_up, geneset_dn, disease_up, disease_dn)
        
    return PCP_geneset_df

def PCP_geneset_allMsigDB(cgp_updn, cgp_updn_labels, disease_sig_up_dict, disease_sig_dn_dict):
    
    PCP_geneset_allMsigDB_df = pd.DataFrame(index=cgp_updn_labels, columns=sorted(disease_sig_up_dict.keys()), dtype='float')
        
    for geneset_name in tqdm(cgp_updn_labels, position=0, leave=True):
        geneset_up = set(cgp_updn[geneset_name + '_UP'])
        geneset_dn = set(cgp_updn[geneset_name + '_DN'])        
        
        for d in disease_sig_up_dict.keys():
            disease_up = set(disease_sig_up_dict[d])
            disease_dn = set(disease_sig_dn_dict[d])    

            PCP_geneset_allMsigDB_df.at[geneset_name, d] = SCS(geneset_up, geneset_dn, disease_up, disease_dn)
            
    return PCP_geneset_allMsigDB_df    

def PCP_perturbation(perturbation, allcells, edgelist_df_dict, disease_sig_up_dict, disease_sig_dn_dict):

    PCP_perturbation_df = pd.DataFrame(index=allcells['Cell_type'], columns=sorted(disease_sig_up_dict.keys()), dtype='float')
    
    for c in allcells['Cell_type']:

        cell_drugs = sorted(list(set(edgelist_df_dict[c]['Drug'])))

        if perturbation in cell_drugs:

            cell_sorted_df = edgelist_df_dict[c].sort_values(['Drug', 'Target']) 
            drug_up = set(cell_sorted_df[(cell_sorted_df['Drug']==perturbation) & (cell_sorted_df['Direction']=='Up')]['Target'].values)
            drug_dn = set(cell_sorted_df[(cell_sorted_df['Drug']==perturbation) & (cell_sorted_df['Direction']=='Down')]['Target'].values)

            for d in disease_sig_up_dict.keys():
                disease_up = set(disease_sig_up_dict[d])
                disease_dn = set(disease_sig_dn_dict[d])    

                PCP_perturbation_df.at[c, d] = SCS(drug_up, drug_dn, disease_up, disease_dn)

        else:
            PCP_perturbation_df = PCP_perturbation_df.drop(c)
            
    return PCP_perturbation_df

def PCP_perturbation_alldrugs(alldrugs, allcells, edgelist_df_dict, disease_sig_up_dict, disease_sig_dn_dict, 
                              proj_path, method_name, return_output=False):
    
    if not os.path.exists(proj_path + method_name + '_pcp_perturbation_df_dict.pickle'):
        
        print('Calculating PCPs for all perturbations...', flush=True)
        PCP_perturbation_df_dict = {}

        for drug in tqdm(alldrugs['Pert_iname'], position=0, leave=True):

            PCP_perturbation_df_dict[drug] = PCP_perturbation(drug, allcells, edgelist_df_dict, disease_sig_up_dict, disease_sig_dn_dict)

        with open(proj_path + method_name + '_pcp_perturbation_df_dict.pickle', 'wb') as fp:
            pickle.dump(PCP_perturbation_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
        if return_output==True:  
            return PCP_perturbation_df_dict
                
    else:
        
        print('Reading PCPs for all perturbations...', flush=True)
        with open(proj_path + method_name + '_pcp_perturbation_df_dict.pickle', 'rb') as fp:
            PCP_perturbation_df_dict = pickle.load(fp)      
    
        if return_output==True:  
            return PCP_perturbation_df_dict

def calculate_auroc(TPR, FPR):
    
    dTPR = np.concatenate((np.ediff1d(TPR), [0]))
    dFPR = np.concatenate((np.ediff1d(FPR), [0]))
    
    return sum(TPR * dFPR) + sum(dTPR * dFPR)/2

def AUROC_AUPRC(ranked_array, positives, auc_threshold=5, binsize=1):

    ranked_array_positives = set(ranked_array) & positives                 
    
    if len(ranked_array_positives) >= auc_threshold:
        
        bins = np.arange(1, len(ranked_array), binsize)      
        TPR = np.zeros(len(bins))
        FPR = np.zeros(len(bins))
        precision = np.zeros(len(bins))
        recall = np.zeros(len(bins))    

        for i, n in enumerate(bins):

            topN = ranked_array[0:n]

            overlap = set(topN) & ranked_array_positives

            TP = 1.0 * len(overlap)
            FP = len(topN) - TP
            FN = len(ranked_array_positives) - TP
            TN = len(ranked_array) - (TP + FP + FN)
            TPR[i] = TP / (TP + FN)
            FPR[i] = FP / (FP + TN)
            precision[i] = TP / (TP + FP)
            recall[i] = TP / (TP + FN)

        auroc = calculate_auroc(TPR, FPR)
        auprc = calculate_auroc(precision, recall)
    
        return auroc, auprc
    
    else:
        return np.nan, np.nan

def PACOS(PCP_geneset_df, PCP_perturbation_df_dict, alldrugs, allcells, tool_scores, 
          proj_path, method_name, geneset_name, r=2.0, threshold=10, tqdm_off=False, messages=False):

    if not os.path.exists(proj_path +  '%s_%s_r%s.csv' % (method_name, geneset_name, int(r))):
        
        if messages: print('Running Pathopticon...', flush=True)
        PACOS_spearman_rho_df = pd.DataFrame(index=alldrugs['Pert_iname'], columns=allcells['Cell_type'], dtype='float')
        PACOS_spearman_pval_df = pd.DataFrame(index=alldrugs['Pert_iname'], columns=allcells['Cell_type'], dtype='float')
        
        for drug in tqdm(alldrugs['Pert_iname'], position=0, leave=True, disable=tqdm_off):

            for c in PCP_perturbation_df_dict[drug].index:
                common_diseases = set(PCP_geneset_df.loc[geneset_name][~pd.isnull(PCP_geneset_df.loc[geneset_name])].index) &\
                                    set(PCP_perturbation_df_dict[drug].loc[c][~pd.isnull(PCP_perturbation_df_dict[drug].loc[c])].index)
                if len(common_diseases) >= threshold:
                    PACOS_spearman_rho_df.at[drug, c], PACOS_spearman_pval_df.at[drug, c] = st.spearmanr(PCP_geneset_df.loc[geneset_name][common_diseases], 
                                                                      PCP_perturbation_df_dict[drug].loc[c][common_diseases])

        PACOS_spearman_rho_df[pd.isnull(PACOS_spearman_rho_df)] = -666
        sorted_ix = np.argsort(PACOS_spearman_rho_df.values.flatten())
        nonan_len = len(PACOS_spearman_rho_df.values.flatten()[PACOS_spearman_rho_df.values.flatten()!=-666])
        nonan_ranked = PACOS_spearman_rho_df.values.flatten()[np.flip(sorted_ix)[0:nonan_len]]
        nonan_ranked_ix = np.flip(sorted_ix)[0:nonan_len]
        nonan_ranked_i, nonan_ranked_j = np.unravel_index(nonan_ranked_ix, PACOS_spearman_rho_df.shape)

        PACOS_spearman_rho_all_ranked = pd.DataFrame()
        PACOS_spearman_rho_all_ranked['Pert_iname'] = PACOS_spearman_rho_df.index[nonan_ranked_i]
        PACOS_spearman_rho_all_ranked['Cell_type'] = PACOS_spearman_rho_df.columns[nonan_ranked_j]
        PACOS_spearman_rho_all_ranked['PACOS_Spearman_rho'] = np.array([PACOS_spearman_rho_df.at[PACOS_spearman_rho_df.index[i], 
                                                                                                 PACOS_spearman_rho_df.columns[j]] 
                                                                        for i, j in zip(nonan_ranked_i, nonan_ranked_j)])
        PACOS_spearman_rho_all_ranked['PACOS_Spearman_pval'] = np.array([PACOS_spearman_pval_df.at[PACOS_spearman_rho_df.index[i], 
                                                                                                   PACOS_spearman_rho_df.columns[j]] 
                                                                        for i, j in zip(nonan_ranked_i, nonan_ranked_j)]) 


        PACOS_tool_merged_df = pd.merge(PACOS_spearman_rho_all_ranked, tool_scores, left_on='Pert_iname', right_on='Pert_iname', how='left')
        PACOS_tool_merged_df['PACOS_Spearman_rho_reverse'] = -1.0*PACOS_tool_merged_df['PACOS_Spearman_rho']
        PACOS_tool_merged_df['tool_score_imputed'] = PACOS_tool_merged_df['tool score scaled'].fillna(tool_scores['tool score scaled'].median())   
        PACOS_tool_merged_df['PACOS_tool_combined'] = (r*PACOS_tool_merged_df['PACOS_Spearman_rho'] + 
                                                       PACOS_tool_merged_df['tool_score_imputed']) / (r + 1.0)
        PACOS_tool_merged_df['PACOS_tool_combined_reverse'] = (-1.0*r*PACOS_tool_merged_df['PACOS_Spearman_rho'] + 
                                                       PACOS_tool_merged_df['tool_score_imputed']) / (r + 1.0)


        PACOS_tool_merged_df.to_csv(proj_path +  '%s_%s_r%s.csv' % (method_name, geneset_name, int(r)), index=False)
        
        return PACOS_tool_merged_df
        
    else:
        
        PACOS_tool_merged_df = pd.read_csv(proj_path +  '%s_%s_r%s.csv' % (method_name, geneset_name, int(r)))
        
        return PACOS_tool_merged_df

def PACOS_cell_AUC(PACOS_tool_merged_df, allcells, positives, 
                   models=['PACOS_Spearman_rho', 'PACOS_Spearman_rho_reverse', 'PACOS_tool_combined', 'PACOS_tool_combined_reverse'],
                   auc_params={'auc_threshold':5, 'binsize':1}):
    
    model_auroc_df = pd.DataFrame(index=allcells['Cell_type'], columns=['%s_AUROC' % m for m in models])
    model_auprc_df = pd.DataFrame(index=allcells['Cell_type'], columns=['%s_AUPRC' % m for m in models])
    
    for m in models:
        
        allcells_sorted = PACOS_tool_merged_df.sort_values(['Cell_type', m], ascending=[True, False])
        
        for c in allcells['Cell_type']:  
            cell_sorted = allcells_sorted[allcells_sorted['Cell_type']==c]['Pert_iname'].values
            model_auroc_df.at[c, '%s_AUROC' % m], model_auprc_df.at[c, '%s_AUPRC' % m] = AUROC_AUPRC(cell_sorted, positives, 
                                                                                                                   auc_threshold=auc_params['auc_threshold'], 
                                                                                                                   binsize=auc_params['binsize'])
        
    return model_auroc_df, model_auprc_df

def PACOS_cell_AUC_randomize(PACOS_tool_merged_df, allcells, positives, 
                             models=['PACOS_Spearman_rho', 'PACOS_Spearman_rho_reverse', 'PACOS_tool_combined', 'PACOS_tool_combined_reverse'],
                             auc_params={'auc_threshold':5, 'binsize':1}, 
                             Nrand=200, tqdm_off=False):

    rand_model_auroc_df_dict = {}
    rand_model_auprc_df_dict = {}
    
    for m in models:
        
        cell_geneset_auroc_df = pd.DataFrame(index=allcells['Cell_type'], columns=np.arange(Nrand))
        cell_geneset_auprc_df = pd.DataFrame(index=allcells['Cell_type'], columns=np.arange(Nrand))        

        allcells_sorted = PACOS_tool_merged_df.sort_values(['Cell_type', m], ascending=[True, False])
        
        for nrand in tqdm(np.arange(Nrand), position=0, leave=True, disable=tqdm_off):
            for c in allcells['Cell_type']:  
                cell_sorted = allcells_sorted[allcells_sorted['Cell_type']==c]['Pert_iname'].values
                np.random.shuffle(cell_sorted)
                cell_geneset_auroc_df.at[c, nrand], cell_geneset_auprc_df.at[c, nrand] = AUROC_AUPRC(cell_sorted, positives, 
                                                                                                     auc_threshold=auc_params['auc_threshold'], 
                                                                                                     binsize=auc_params['binsize'])
            
        rand_model_auroc_df_dict[m] = cell_geneset_auroc_df
        rand_model_auprc_df_dict[m] = cell_geneset_auprc_df 
        
    return rand_model_auroc_df_dict, rand_model_auprc_df_dict

def PACOS_nested_prioritization(PACOS_tool_merged_df, model_auroc_df, rand_model_auroc_df_dict):
    
    emp_pval_df = pd.DataFrame(index=model_auroc_df.index, columns=['%s_AUROC_emp_pval' % m for m in list(rand_model_auroc_df_dict.keys())])
    
    for m in rand_model_auroc_df_dict.keys():
    
        nonan_cells = model_auroc_df['%s_AUROC' % m][~pd.isnull(model_auroc_df['%s_AUROC' % m])].index
        for c in nonan_cells:
            emp_pval_df.at[c, '%s_AUROC_emp_pval' % m] = ((rand_model_auroc_df_dict[m].loc[c].values > 
                                                             model_auroc_df['%s_AUROC' % m].loc[c]).sum() /
                                                            float(np.shape(rand_model_auroc_df_dict[m])[1]))
    
    model_auroc_pval_df = pd.merge(model_auroc_df, emp_pval_df, left_index=True, right_index=True)
    PACOS_nested_df = pd.merge(PACOS_tool_merged_df, model_auroc_pval_df, left_on='Cell_type', 
                               right_index=True).dropna(subset=['%s_AUROC' % m for m in list(rand_model_auroc_df_dict.keys())])        
        
    return emp_pval_df, PACOS_nested_df

def run_benchmark(proj_path, benchmark_path, method_name, 
                  PACOS_Spearman_threshold=10, r=2.0,  Nrand=100, auc_params={'auc_threshold':5, 'binsize':1}, 
                  models=['PACOS_Spearman_rho', 'PACOS_Spearman_rho_reverse', 'PACOS_tool_combined', 'PACOS_tool_combined_reverse']):
    
    print('Loading input data...', flush=True)
    input_paths_dict = input_paths(proj_path)
    
    L1000_gene_info, L1000_cell_info, L1000_inst_info, L1000_pert_info, L1000_sig_info = import_L1000_metadata(input_paths_dict['L1000_gene_info_path'], 
                                                                                                               input_paths_dict['L1000_cell_info_path'],
                                                                                                               input_paths_dict['L1000_inst_info_path'], 
                                                                                                               input_paths_dict['L1000_pert_info_path'], 
                                                                                                               input_paths_dict['L1000_sig_info_path'])    
    
    QUIZC_cids_inchi_smiles = pd.read_csv(input_paths_dict['QUIZC_cids_inchi_smiles_path'])
    QUIZC_activityStats_nooutliers_df_besttool = pd.read_csv(input_paths_dict['tool_path'])
        
    cgp_updn, cgp_updn_labels, cgp_updn_allgenes = process_MSigDB_CGP(input_paths_dict['cgp_dir'])
    Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn = import_Enrichr_GEO(input_paths_dict['Enrichr_GEO_up_path'], 
                                                                                    input_paths_dict['Enrichr_GEO_dn_path'])
    TTD_drugs, TTD_InChI2CID, TTD_drugs_CID, TTD_drug_target, TTD_targID_dict = import_TDD(input_paths_dict['TTD_drugs_path'], 
                                                                                           input_paths_dict['TTD_InChI2CID_path'],
                                                                                           input_paths_dict['TTD_drug_target_path'], 
                                                                                           input_paths_dict['TTD_target_path'])

    geneset_pert_iname_dict = get_MSigDB_geneset_targets(cgp_updn_labels, cgp_updn, TTD_targID_dict, TTD_drug_target, 
                                                         TTD_drugs_CID, QUIZC_cids_inchi_smiles)   
    
    
    if method_name == 'QUIZ-C':
        print('Processing QUIZ-C networks...', flush=True)
        edgelist_df_dict, nodelist_df_dict, allcells, allgenes, alldrugs, allnodes = process_QUIZC_output(input_paths_dict['pos_edges_dict_path'], 
                                                                                                         input_paths_dict['neg_edges_dict_path'],
                                                                                                         input_paths_dict['drugs_dict_path'],
                                                                                                         L1000_gene_info)
    elif (method_name == 'MODZ') | (method_name == 'CD'):
        print('Importing gene-perturbation networks...', flush=True)
        edgelist_df_dict, nodelist_df_dict, allcells, allgenes, alldrugs, allnodes = import_CD_MODZ_networks(input_paths_dict['%s_networks_path' % method_name])  

    
    pcp_perturbation_df_dict = PCP_perturbation_alldrugs(alldrugs, allcells, edgelist_df_dict, 
                                                         Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn,
                                                         benchmark_path, method_name=method_name, return_output=True) 
                                                                                                               
    print('Calculating PCPs for all benchmark gene sets...', flush=True)                                                                                                                
    pcp_geneset_allMsigDB_df = PCP_geneset_allMsigDB(cgp_updn, cgp_updn_labels, Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn)
    
    
    print('Running the benchmark...', flush=True)      
    benchmark_genesets = pd.read_csv(input_paths_dict['benchmark_genesets_path'], header=None).rename(columns={0: 'Geneset_name'})
                                                                                                               
    for geneset_name in tqdm(benchmark_genesets['Geneset_name'].to_numpy(), position=0, leave=True):
        
        if not os.path.exists(benchmark_path +  'Nested_%s_%s_r%s.csv' % (method_name, geneset_name, int(r))):
        
            geneset_up = set(cgp_updn[geneset_name + '_UP'])
            geneset_dn = set(cgp_updn[geneset_name + '_DN'])   
            pcp_geneset_df = PCP_geneset(geneset_up, geneset_dn, Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn, geneset_name=geneset_name)
            
            pacos_tool_merged_df = PACOS(pcp_geneset_df, pcp_perturbation_df_dict, alldrugs, allcells,
                                         QUIZC_activityStats_nooutliers_df_besttool, benchmark_path, method_name=method_name, geneset_name=geneset_name, 
                                         r=r, threshold=PACOS_Spearman_threshold, tqdm_off=True)

            model_auroc_df, model_auprc_df = PACOS_cell_AUC(pacos_tool_merged_df, allcells, geneset_pert_iname_dict[geneset_name], 
                                                            models=models,
                                                            auc_params={'auc_threshold':5, 'binsize':1})   


            rand_model_auroc_df_dict, rand_model_auprc_df_dict = PACOS_cell_AUC_randomize(pacos_tool_merged_df, allcells, 
                                                                                          geneset_pert_iname_dict[geneset_name], 
                                                                                          models=models,
                                                                                          auc_params={'auc_threshold':5, 'binsize':1}, 
                                                                                          Nrand=Nrand, tqdm_off=True) 

            emp_pval_df, pacos_nested_df = PACOS_nested_prioritization(pacos_tool_merged_df, model_auroc_df, rand_model_auroc_df_dict)    

            pacos_nested_df.to_csv(benchmark_path +  'Nested_%s_%s_r%s.csv' % (method_name, geneset_name, int(r)), index=False)
    
            
    print('Benchmark complete.', flush=True)                     

def get_topN_drugs(ranked_df_path, benchmark_genesets, method_name='QUIZ-C', model='PACOS_tool_combined', topN=50):
    
    top50_alldrugs = set()
    if method_name in ['QUIZ-C', 'MODZ', 'CD']:
        for geneset_name in tqdm(benchmark_genesets['Geneset_name'].to_numpy(), position=0, leave=True):      
            ranked_df = pd.read_csv(ranked_df_path + '%s_%s_r2.csv' % (method_name, geneset_name))
            combined_topN = ranked_df.sort_values(model, ascending=False)['Pert_iname'].unique()[0:topN]
            top50_alldrugs.update(set(combined_topN))
                        
    elif method_name == 'L1000CDS2':
        with open(ranked_df_path + 'L1000CDS2_result_df_dict.pickle', 'rb') as handle:
            L1000CDS2_result_df_dict = pickle.load(handle)
        for geneset_name in tqdm(benchmark_genesets['Geneset_name'].to_numpy(), position=0, leave=True):      
            top50_alldrugs.update(set(L1000CDS2_result_df_dict[geneset_name]['pert_iname']))
           
    return top50_alldrugs

def get_topN_tanimoto(ranked_df_path, benchmark_genesets, apfp_Tanimoto_pert_iname, 
                      method_name='QUIZ-C', model='PACOS_tool_combined', topN=50):
    
    top50_tanimoto_dict = {} 
    if method_name in ['QUIZ-C', 'MODZ', 'CD']:
        for geneset_name in tqdm(benchmark_genesets['Geneset_name'].to_numpy(), position=0, leave=True):
            ranked_df = pd.read_csv(ranked_df_path + '%s_%s_r2.csv' % (method_name, geneset_name))
            combined_topN = ranked_df.sort_values(model, ascending=False)['Pert_iname'].unique()[0:topN]
            temp_tanimoto_df = apfp_Tanimoto_pert_iname[(apfp_Tanimoto_pert_iname['Pert_iname_x'].isin(set(combined_topN))) & 
                                                              (apfp_Tanimoto_pert_iname['Pert_iname_y'].isin(set(combined_topN)))].drop_duplicates()
            top50_tanimoto_dict[geneset_name] = temp_tanimoto_df['Tanimoto coefficient'].to_numpy()
            
    elif method_name == 'L1000CDS2':
        with open(ranked_df_path + 'L1000CDS2_result_df_dict.pickle', 'rb') as handle:
            L1000CDS2_result_df_dict = pickle.load(handle)        
        for geneset_name in tqdm(benchmark_genesets['Geneset_name'].to_numpy(), position=0, leave=True): 
            rank_df = L1000CDS2_result_df_dict[geneset_name]
            combined_topN = rank_df['pert_iname'].unique()[0:topN]

            temp_tanimoto_df = apfp_Tanimoto_pert_iname[(apfp_Tanimoto_pert_iname['Pert_iname_x'].isin(set(combined_topN))) & 
                                                          (apfp_Tanimoto_pert_iname['Pert_iname_y'].isin(set(combined_topN)))].drop_duplicates()
            top50_tanimoto_dict[geneset_name] = temp_tanimoto_df['Tanimoto coefficient'].to_numpy()
            
    return top50_tanimoto_dict
    
def get_activity_stats(chembl27_QUIZC_CIDs):
    
    activityStats_df = pd.DataFrame()
    for cid in tqdm(chembl27_QUIZC_CIDs['PubChem CID'].sort_values().drop_duplicates(), position=0, leave=True):

        temp_cid_df = chembl27_QUIZC_CIDs[(chembl27_QUIZC_CIDs['Standard Units'] == 'nM') &
                                         (chembl27_QUIZC_CIDs['Standard Relation'] == "'='") &
                                         (chembl27_QUIZC_CIDs['PubChem CID']==cid)]
        cid_targets = set(temp_cid_df['NCBI Gene ID(supplied by NCBI)'])
        temp_on_nM = {}
        temp_off_nM = {}
        for on_t in cid_targets:
            off_t = cid_targets - set([on_t])
            temp_on_nM[on_t] = temp_cid_df[temp_cid_df['NCBI Gene ID(supplied by NCBI)']==on_t]['Standard Value'].values
            temp_off_nM[on_t] = temp_cid_df[temp_cid_df['NCBI Gene ID(supplied by NCBI)'].isin(off_t)]['Standard Value'].values

            # remove measurements with negative nM values, if any   
            temp_on_nM[on_t] = temp_on_nM[on_t][temp_on_nM[on_t] > 0]
            temp_off_nM[on_t] = temp_off_nM[on_t][temp_off_nM[on_t] > 0]

            if (len(temp_on_nM[on_t]) > 0) & (len(temp_off_nM[on_t]) > 0):

                if (len(temp_on_nM[on_t][temp_on_nM[on_t] <= 100.0]) > 1):
                    strength = 7
                elif (len(temp_on_nM[on_t][temp_on_nM[on_t] <= 1000.0]) > 4):
                    strength = 4
                else:
                    strength = 1

                data_bias = len(temp_on_nM[on_t]) / (len(temp_on_nM[on_t]) + len(temp_off_nM[on_t]))
                Q1_diff = np.quantile(np.log10(temp_off_nM[on_t]), 0.25) - np.quantile(np.log10(temp_on_nM[on_t]), 0.25)
                Q3_diff = np.quantile(np.log10(temp_off_nM[on_t]), 0.75) - np.quantile(np.log10(temp_on_nM[on_t]), 0.75)

                KS_pval = st.ks_2samp(temp_on_nM[on_t], temp_off_nM[on_t])[1]
                if (KS_pval<1e-16): KS_pval = 1e-16

                selectivity = ((Q1_diff/3) + (-np.log10(KS_pval)/16.0) + (1 - data_bias)) / 3
                tool_score = selectivity * strength

                activityStats_df = activityStats_df.append(pd.DataFrame([cid, on_t, KS_pval, Q1_diff, Q3_diff, data_bias, 
                                                                         selectivity, strength, tool_score]).T)

    activityStats_df = activityStats_df.reset_index(drop=True).rename(columns={0: 'cid', 1: 'Target Entrez ID', 2: 'K-S p-value', 
                                                                               3: 'Q1_diff', 4: 'Q3_diff',  5: 'data_bias', 
                                                                               6: 'selectivity score', 7: 'strength', 8: 'tool score'})
    
    return activityStats_df

def calculate_selectivity_tool(activityStats_df, QUIZC_cids_inchi_smiles, iqr_cutoff, tool_outlier_cutoff, proj_intermediary_outputs_path):
     
    # Based on the distributions, remove outliers.
    selectivity_Q1 =  activityStats_df['selectivity score'].quantile(0.25)
    selectivity_Q3 =  activityStats_df['selectivity score'].quantile(0.75)
    selectivity_IQR = selectivity_Q3 - selectivity_Q1
    selectivity_outliers_ix = activityStats_df[((activityStats_df['selectivity score'] < (selectivity_Q1 - (selectivity_IQR  * iqr_cutoff))) | 
                                   (activityStats_df['selectivity score'] > (selectivity_Q3 + (selectivity_IQR  * iqr_cutoff))))].index

    tool_outliers_ix = activityStats_df[activityStats_df['tool score'] >= tool_outlier_cutoff].index

    nooutliers_ix = set(activityStats_df.index) - (set(selectivity_outliers_ix) | set(tool_outliers_ix))

    activityStats_df_nooutliers = activityStats_df.loc[nooutliers_ix]
    
    
    # Scale selectivity and tool scores
    sel_scale_len = (activityStats_df_nooutliers['selectivity score'].max() - activityStats_df_nooutliers['selectivity score'].min()) / 2
    sel_scale_diff = activityStats_df_nooutliers['selectivity score'].max() - sel_scale_len
    activityStats_df_nooutliers['selectivity score scaled'] = (activityStats_df_nooutliers['selectivity score'] - sel_scale_diff) / sel_scale_len

    tool_scale_len = (activityStats_df_nooutliers['tool score'].max() - activityStats_df_nooutliers['tool score'].min()) / 2
    tool_scale_diff = activityStats_df_nooutliers['tool score'].max() - tool_scale_len
    activityStats_df_nooutliers['tool score scaled'] = (activityStats_df_nooutliers['tool score'] - tool_scale_diff) / tool_scale_len
    
    
    # Merge with QUIZC_cids_inchi_smiles to get the pert_iname. Note that when we are writing to file we are collapsing drug-target pairs onto
    # the drug space. When we do so, we take the drug-target pairs with either the highest tool scores (first one below) or selectivity scores
    # (second one below). The resulting file name clarifies this choice.
    QUIZC_activityStats_nooutliers_df_besttool = pd.merge(activityStats_df_nooutliers.sort_values(['cid', 'tool score scaled'], ascending=[True, False])
                                                              .drop_duplicates('cid'), 
                                                              QUIZC_cids_inchi_smiles[['Pert_iname', 'pubchem_cid_x']], 
                                                              left_on='cid', right_on='pubchem_cid_x', how='inner').drop_duplicates()

    QUIZC_activityStats_nooutliers_df_besttool.to_csv(proj_intermediary_outputs_path + 'QUIZC_activityStats_nooutliers_df_besttool.csv', 
                                                     index=False)

    QUIZC_activityStats_nooutliers_df_bestselectivity = pd.merge(activityStats_df_nooutliers.sort_values(['cid', 'selectivity score scaled'], 
                                                                                                        ascending=[True, False])
                                                                .drop_duplicates('cid'), 
                                                                QUIZC_cids_inchi_smiles[['Pert_iname', 'pubchem_cid_x']], 
                                                                left_on='cid', right_on='pubchem_cid_x', how='inner').drop_duplicates()

    QUIZC_activityStats_nooutliers_df_bestselectivity.to_csv(proj_intermediary_outputs_path + 'QUIZC_activityStats_nooutliers_df_bestselectivity.csv', index=False)
    
    return QUIZC_activityStats_nooutliers_df_besttool, QUIZC_activityStats_nooutliers_df_bestselectivity
    
    
# Fisher's Exact test-based overlap
def two_geneset_overlap(genesetA_dict, genesetB_dict, universe_geneset, outfile):

    overlap_df = pd.DataFrame(index=genesetA_dict.keys(), columns=genesetB_dict.keys())
    fisher_pval_df = pd.DataFrame(index=genesetA_dict.keys(), columns=genesetB_dict.keys())
    fisher_OR_df = pd.DataFrame(index=genesetA_dict.keys(), columns=genesetB_dict.keys())

    for dt1, dt2 in product(genesetA_dict.keys(), genesetB_dict.keys()):

        A = len(set(genesetA_dict[dt1]) & set(genesetB_dict[dt2]))
        B = len(set(genesetA_dict[dt1]) - set(genesetB_dict[dt2]))
        C = len(set(genesetB_dict[dt2]) - set(genesetA_dict[dt1]))
        D = len(universe_geneset) - len(set(genesetA_dict[dt1]) | set(genesetB_dict[dt2]))

        overlap_df.at[dt1, dt2] = A
        fisher_OR_df.at[dt1, dt2] = st.fisher_exact([[A, B], [C, D]])[0]
        fisher_pval_df.at[dt1, dt2] = st.fisher_exact([[A, B], [C, D]])[1]

    overlap_df.to_csv(outfile + '_overlap.csv')
    fisher_OR_df.to_csv(outfile + '_fisher_OR.csv')
    fisher_pval_df.to_csv(outfile + '_fisher_pval.csv')
   
    return overlap_df, fisher_OR_df, fisher_pval_df
    
    
# Fisher's Exact test-based overlap
def fisher_pathway_enrichment(pathway_geneset_dict, input_geneset_dict, universe_geneset, method, outfile):

    overlap_df = pd.DataFrame(index=pathway_geneset_dict.keys(), columns=input_geneset_dict.keys())
    fisher_pval_df = pd.DataFrame(index=pathway_geneset_dict.keys(), columns=input_geneset_dict.keys())
    fisher_adjpval_df = pd.DataFrame(index=pathway_geneset_dict.keys(), columns=input_geneset_dict.keys())
    fisher_OR_df = pd.DataFrame(index=pathway_geneset_dict.keys(), columns=input_geneset_dict.keys())

    for dt1, dt2 in product(pathway_geneset_dict.keys(), input_geneset_dict.keys()):

        A = len(set(pathway_geneset_dict[dt1]) & set(input_geneset_dict[dt2]))
        B = len(set(pathway_geneset_dict[dt1]) - set(input_geneset_dict[dt2]))
        C = len(set(input_geneset_dict[dt2]) - set(pathway_geneset_dict[dt1]))
        D = len(universe_geneset) - len(set(pathway_geneset_dict[dt1]) | set(input_geneset_dict[dt2]))

        overlap_df.at[dt1, dt2] = A
        fisher_OR_df.at[dt1, dt2] = st.fisher_exact([[A, B], [C, D]])[0]
        fisher_pval_df.at[dt1, dt2] = st.fisher_exact([[A, B], [C, D]])[1]

    for dt2 in input_geneset_dict.keys():
        fisher_adjpval_df[dt2] = pd.Series(sm.stats.multipletests(fisher_pval_df[dt2], alpha=0.05, method=method)[1], index=fisher_pval_df.index)
        
    overlap_df.to_csv(outfile + '_overlap.csv')
    fisher_OR_df.to_csv(outfile + '_fisher_OR.csv')
    fisher_pval_df.to_csv(outfile + '_fisher_pval.csv')
    fisher_adjpval_df.to_csv(outfile + '_fisher_adjpval.csv')
   
    return overlap_df, fisher_OR_df, fisher_pval_df, fisher_adjpval_df
    
    