from MuXTalk_master import *

if __name__ == "__main__":  

	run_MuXTalk(proj_path=sys.argv[1], 
            input_filenames_dict = {'HUGO_symb_entrez_uniprot': 'HugoGene_20200528.txt', 'PPI_Cheng_2019_data': 'PPI_Cheng_NatComms2019.csv',
                       'KEGG_all_nodes_df': 'KEGG_expanded_all_nodes.csv', 'KEGG_all_edges_df': 'KEGG_expanded_all_edges.csv',
                        'df_motinf' : 'cisbpall_motinf.txt', 'XTalk_DB': 'XTalkDB_crosstalk.csv', 
                        'almen_etal': 'almen_etal_12915_2009_258_MOESM1_ESM.csv', 
                        'lambert_etal': 'lambert_etal_1-s2.0-S0092867418301065-mmc2_TableS1.csv'},
            input_GRN=sys.argv[2], 
            MuXTalk_method=sys.argv[3], 
            sp_threshold=sys.argv[4], 
            parquet=sys.argv[5])