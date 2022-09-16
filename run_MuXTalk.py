import argparse, sys
from MuXTalk_master import *

if __name__ == "__main__":

	parser=argparse.ArgumentParser()
	
	parser.add_argument("--proj_path", help="Path to the project file where all the input/output files will be located.")
	parser.add_argument("--input_GRN", help="Name of the input gene regulatory network (GRN). Please note that this must be the same as the name preceding '_edges.csv' in the GRN edgelist.")
	parser.add_argument("--MuXTalk_method", help="MuXTalk method to be used: The two valid options are 'MuXTalk_between' and 'MuXTalk_shortest'.")
	parser.add_argument("--get_n", help="Number of randomized instances to be retrieved from the ensemble. For memory efficiency, get_n=100 is recommended."
	parser.add_argument("--get_randomly", help="Option to retrieve randomized networks in order or randomly. Default value is True."
	parser.add_argument("--sp_threshold", help="The shortest path threshold to be used in MuXTalk. The valid options are 'None', '1', or '2'.")
	parser.add_argument("--parquet", help="Option to output the MuXTalk results as a .parquet file. False (default) outputs .csv files.")

	args=parser.parse_args()

	run_MuXTalk_npz(proj_path=args.proj_path, 
            input_filenames_dict = {'HUGO_symb_entrez_uniprot': 'HugoGene_20200528.txt', 'PPI_Cheng_2019_data': 'PPI_Cheng_NatComms2019.csv',
                       'KEGG_all_nodes_df': 'KEGG_expanded_all_nodes.csv', 'KEGG_all_edges_df': 'KEGG_expanded_all_edges.csv',
                        'df_motinf' : 'cisbpall_motinf.txt', 'XTalk_DB': 'XTalkDB_crosstalk.csv', 
                        'almen_etal': 'almen_etal_12915_2009_258_MOESM1_ESM.csv', 
                        'lambert_etal': 'lambert_etal_1-s2.0-S0092867418301065-mmc2_TableS1.csv'},
            npz_filenames_dict = {'GRN': 'A_GRN_sparr_rand_npz_files/', 'KEGG_e': 'A_KEGG_e_sparr_rand_npz_files/', 
										  'KEGGPPI': 'A_KEGGPPI_sparr_rand_npz_files/'},
            input_GRN=args.input_GRN, 
            MuXTalk_method=args.MuXTalk_method, 
            get_n=args.get_n, 
			get_randomly=args.get_randomly, 
            sp_threshold=args.sp_threshold, 
            parquet=args.parquet)