import streamlit as st
import streamlit.components.v1 as components
from MuXTalk_functions_Streamlit_Docker import *
import networkx as nx
from pyvis.network import Network

st.set_page_config(layout='wide')

st.sidebar.title(
	'''
	MuXTalk
	***	
	'''
)

input_GRN = st.sidebar.selectbox('Select GRN', ['HumanGRN10e6', 'HumanGRN10e5', 'HumanGRN10e4'])	
MuXTalk_method = st.sidebar.selectbox('Select MuXTalk type', ['MuXTalk_shortest', 'MuXTalk_between'])
sp_threshold = st.sidebar.selectbox('Select shortest path threshold (MuXTalk_shortest only)', [2, 1, 'None'])
st.sidebar.write(''' *** ''')
custom_GRN_edges = st.sidebar.file_uploader('Custom GRN edges')
custom_GRN_parquet = st.sidebar.file_uploader('Custom GRN parquet file')

with st.container():
	'''
	# **MuXTalk:** Detecting and dissecting signaling crosstalk via multilayer networks	
	***
	'''

proj_path = '/MuXTalk_Streamlit_app/'
input_filenames_dict = {'HUGO_symb_entrez_uniprot': 'HugoGene_20200528.txt', 'PPI_Cheng_2019_data': 'PPI_Cheng_NatComms2019.csv',
					   'KEGG_all_nodes_df': 'KEGG_expanded_all_nodes.csv', 'KEGG_all_edges_df': 'KEGG_expanded_all_edges.csv',
						'df_motinf' : 'cisbpall_motinf.txt', 'XTalk_DB': 'XTalkDB_crosstalk.csv', 
						'almen_etal': 'almen_etal_12915_2009_258_MOESM1_ESM.csv', 
						'lambert_etal': 'lambert_etal_1-s2.0-S0092867418301065-mmc2_TableS1.csv'
					   }
					   
multilink_annot = {}
multilink_annot['signaling'] = {0: 'No edge', 1: 'Activation', 2: 'Binding/association', 3: 'Compound',
                                4: 'Dephosphorylation', 5: 'Dissociation', 6: 'Indirect effect', 7: 'Inhibition', 
                                8: 'Phosphorylation', 9: 'PPI', 10: 'State change', 11: 'Ubiquitination'}
multilink_annot['regulatory'] = {-1: '(B->A)', 0: 'No edge', 1: '(A->B)'}
multilink_annot_dict = {'%s%s' % (i, j): 'Signaling edge type: %s; Regulatory edge direction: %s' % (multilink_annot['signaling'][i], multilink_annot['regulatory'][j]) for i, j in 
                        product(multilink_annot['signaling'].keys(), multilink_annot['regulatory'].keys())}				
					   
expander_readme = st.expander(label='Read me first!', expanded=False)
with expander_readme:
	readme_msg = st.container()	
	readme_msg.write('''
	
	**MuXTalk: Detecting and dissecting signaling crosstalk via multilayer networks.** A network-based statistical 
	framework to explore crosstalk between signaling pathways.\n
	**Reference:** \n	
	This is the visualization component of the [MuXTalk framework](https://github.com/r-duh/MuXTalk), 
	where users can explore the crosstalk between 60 pairs of KEGG pathways. There are two options to use 
	this MuXTalk visualization app: \n
	1) Using the **default** human gene regulatory networks (GRNs) presented in the manuscript as input. 
	For this option, you can simply select the desired GRN (from the least dense (HumanGRN10e-6) to the densest 
	(HumanGRN10e-6) network) and choose which MuXTalk method to use. \n
	2) Using **custom GRNs** as input. For this  option, you first need to run MuXTalk with the custom GRN following the 
	instructions on the [GitHub page](https://github.com/r-duh/MuXTalk). After running MuXTalk, you can use its output 
	for visualization. You will need to provide two files: \n
	- The custom GRN edgelist file. This file should consist of two columns, without headers, the first one for 
	the source gene (or transcription factor) and the second one for the target gene, 
	as described on the GitHub page. This is also the file you used with MuXTalk as the input_GRN 
	edgelist file ("customGRN_edges.csv"). \n
	- The output of MuXTalk in .parquet format, which includes the list of prioritized pathway pairs predicted by MuXTalk
	("..._detected_discovery.parquet").\n	
	Once the MuXTalk options are selected or custom GRN files are uploaded, this Streamlit app will **initialize** the 
	networks for visualization and display the ranked list of pathway pairs. You can inspect this list for the pathway
	pairs that are the most likely to crosstalk as per the MuXTalk algorithm. \n
	After Pathway A and Pathway B are chosen, the corresponding **significant multilink types** that MuXTalk predicted to 
	mediate the crosstalk will auto-populate the "Select multilink type" dropdown menu. \n
	Choosing the multilink type will create a table that shows the edges of that multilink type and the Gene Symbols and
	Entrez IDs of the nodes that are connected by those edges.\n
	Once the **network visualization** is generated, you can toggle edge descriptions (this feature is off by default for
	visual clarity), zoom in and out and adjust the physics parameters and layout options for the visualization for the 
	desired output.\n
	For the multilink definitions and color codings, refer to the **Legend** box.\n
		
	''')

legend_readme = st.expander(label='Legend', expanded=False)
with legend_readme:
	legend_readme_msg = st.container()	
	legend_readme_msg.write('''
	Multilinks are coded in the form (S, R) where S and R are integers encoding the type of signaling and regulatory interaction.
	R can be 1, -1 or 0, corresponding to a transcriptional regulatory edge (i.e., from a transcription factor to its target) 
	from Pathway A to Pathway B, Pathway B to Pathway A, and no edge, respectively. Signaling edge types S are encoded 
	as below: \n 
	
	| Signaling edge index S | Signaling edge type|
	| --- | --- |
	| 0 | No edge |
	| 1 | Activation |
	| 2 | Binding/association |
	| 3 | Compound |
	| 4 | Dephosphorylation |
	| 5 | Dissociation |
	| 6 | Indirect effect |
	| 7 | Inhibition |
	| 8 | Phosphorylation |
	| 9 | PPI |
	| 10 | State change |
	| 11 | Ubiquitination | \n
	In the visualization, the node and edge color coding is as below:\n
	
	| Color | Edge |
	| --- | --- |
	| Dark grey | Signaling edges from KEGG |
	| Grey |PPI edges |
	| Magenta | Edges of the selected significant multilink type |
	| Pale pink | Edges connecting the significant multilinks to either pathway (for MuXTalk_shortest only) |\n
	
	| Color | Node |
	| --- | --- |
	| Red | Genes in Pathway A |
	| Blue | Genes in Pathway B |
	| Orange |  Genes in both Pathway A and Pathway B |\n
	
	''')

if (custom_GRN_edges is None) & (custom_GRN_parquet is None):
				   
	expander_status = st.expander(label='MuXTalk progress status', expanded=True)
	with expander_status:
		status_msg = st.container()	

	muxtalk_status = status_msg.text('Initializing MuXTalk for the default GRNs...')	

	(KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, PPI_all_edges_entrez, KEGG_PPI_all_edges_entrez, KEGG_PPI_allnodes_entrez_df, 
	 KEGG_interaction_types_dict, KEGG_all_edges_entrez, KEGG_all_paths, KEGG_path_nodes_entrez_dict, 
	 all_motif_types, all_motif_types_list) = process_data(proj_path, input_GRN, input_filenames_dict)

	KEGG_all_paths_sansIL17 = set(KEGG_all_paths) - set(['IL-17 signaling pathway'])
	st.session_state['KEGG_all_paths_sansIL17'] = KEGG_all_paths_sansIL17

	(A_GRN_sparr, A_PPI_sparr, A_KEGGPPI_sparr, A_KEGG_e_sparr_dict) = sparse_layers(KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, 
																					 PPI_all_edges_entrez, 
																					 KEGG_all_edges_entrez, KEGG_PPI_all_edges_entrez, 
																					 KEGG_interaction_types_dict)

	G_PPI =  nx.from_scipy_sparse_array(A_PPI_sparr.astype(bool).astype(int))
	G_PPI_entrez = nx.relabel_nodes(G_PPI, KEGG_PPI_allnodes_entrez_df[[0, 'ix']].to_dict()[0])
	G_KEGGPPI = nx.from_scipy_sparse_array(A_KEGGPPI_sparr.astype(bool).astype(int))
	G_KEGGPPI_entrez = nx.relabel_nodes(G_KEGGPPI, KEGG_PPI_allnodes_entrez_df[[0, 'ix']].to_dict()[0])

	st.session_state['G_PPI_entrez'] = G_PPI_entrez
	st.session_state['G_KEGGPPI_entrez'] = G_KEGGPPI_entrez
	st.session_state['KEGG_path_nodes_entrez_dict'] = KEGG_path_nodes_entrez_dict
	st.session_state['KEGG_PPI_allnodes_entrez_df'] = KEGG_PPI_allnodes_entrez_df
	st.session_state['A_GRN_sparr'] = A_GRN_sparr
	st.session_state['A_KEGG_e_sparr_dict'] = A_KEGG_e_sparr_dict
	st.session_state['A_KEGGPPI_sparr'] = A_KEGGPPI_sparr
	st.session_state['KEGG_interaction_types_dict'] = KEGG_interaction_types_dict

	if MuXTalk_method == 'MuXTalk_between':

		st.session_state['detected_ranked_pathway_pairs_df'] = pd.read_parquet(proj_path + input_GRN + '_between_detected_discovery.parquet')
				
	elif MuXTalk_method == 'MuXTalk_shortest':

		st.session_state['detected_ranked_pathway_pairs_df'] = pd.read_parquet(proj_path + input_GRN + '_shortest_sp%s_detected_discovery.parquet' % sp_threshold)

elif (custom_GRN_edges is not None) & (custom_GRN_parquet is not None):

	expander_status = st.expander(label='MuXTalk progress status', expanded=True)
	with expander_status:
		status_msg = st.container()	

	muxtalk_status = status_msg.text('Initializing MuXTalk for the user-provided GRN...')	

	(KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, PPI_all_edges_entrez, KEGG_PPI_all_edges_entrez, KEGG_PPI_allnodes_entrez_df, 
	 KEGG_interaction_types_dict, KEGG_all_edges_entrez, KEGG_all_paths, KEGG_path_nodes_entrez_dict, 
	 all_motif_types, all_motif_types_list) = process_data_forStreamlit_customGRN(proj_path, custom_GRN_edges, input_filenames_dict)

	KEGG_all_paths_sansIL17 = set(KEGG_all_paths) - set(['IL-17 signaling pathway'])
	st.session_state['KEGG_all_paths_sansIL17'] = KEGG_all_paths_sansIL17

	(A_GRN_sparr, A_PPI_sparr, A_KEGGPPI_sparr, A_KEGG_e_sparr_dict) = sparse_layers(KEGG_PPI_allnodes_entrez, GRN_KEGG_PPI_edges, 
																					 PPI_all_edges_entrez, 
																					 KEGG_all_edges_entrez, KEGG_PPI_all_edges_entrez, 
																					 KEGG_interaction_types_dict)

	G_PPI =  nx.from_scipy_sparse_array(A_PPI_sparr.astype(bool).astype(int))
	G_PPI_entrez = nx.relabel_nodes(G_PPI, KEGG_PPI_allnodes_entrez_df[[0, 'ix']].to_dict()[0])
	G_KEGGPPI = nx.from_scipy_sparse_array(A_KEGGPPI_sparr.astype(bool).astype(int))
	G_KEGGPPI_entrez = nx.relabel_nodes(G_KEGGPPI, KEGG_PPI_allnodes_entrez_df[[0, 'ix']].to_dict()[0])

	st.session_state['G_PPI_entrez'] = G_PPI_entrez
	st.session_state['G_KEGGPPI_entrez'] = G_KEGGPPI_entrez
	st.session_state['KEGG_path_nodes_entrez_dict'] = KEGG_path_nodes_entrez_dict
	st.session_state['KEGG_PPI_allnodes_entrez_df'] = KEGG_PPI_allnodes_entrez_df
	st.session_state['A_GRN_sparr'] = A_GRN_sparr
	st.session_state['A_KEGG_e_sparr_dict'] = A_KEGG_e_sparr_dict
	st.session_state['A_KEGGPPI_sparr'] = A_KEGGPPI_sparr
	st.session_state['KEGG_interaction_types_dict'] = KEGG_interaction_types_dict

	st.session_state['detected_ranked_pathway_pairs_df'] = pd.read_parquet(custom_GRN_parquet)
	
else:

	st.error('Please ensure that both files are provided for the custom GRN.')
	st.stop()

expander_results = st.expander(label='MuXTalk prioritization results', expanded=True)
with expander_results:
	results_cont = st.container()

results_cont.write(st.session_state['detected_ranked_pathway_pairs_df'].reset_index()[['Pathway A', 'Pathway B', 'sig multilink types', 'MuXTalk_score']])  		
																			
form_cols_selectpath = st.columns(2)
p1 = form_cols_selectpath[0].selectbox('Select Pathway A (will be shown as red nodes in the visualization)', sorted(list(st.session_state['KEGG_all_paths_sansIL17'])))
p2 = form_cols_selectpath[1].selectbox('Select Pathway B (will be shown as blue nodes in the visualization)', sorted(list(st.session_state['KEGG_all_paths_sansIL17'])))

if p1 == p2:
	st.error('Pathways A and B must be different.')
	st.stop()

elif '%s-%s' % (p1, p2) not in st.session_state['detected_ranked_pathway_pairs_df'].index:

	st.error('MuXTalk did not find any statistically significant multilink types between these two pathways. Please choose another pair of pathways to visualize.')
	st.stop()	
		
form_cols_selectmultilink = st.columns(1)
multilink_annot_list = ['%s (%s)' %(i, multilink_annot_dict[i]) for i in  st.session_state['detected_ranked_pathway_pairs_df'].loc['%s-%s' % (p1, p2)]['sig multilink types']]
multilink_type = form_cols_selectmultilink[0].selectbox('Select multilink type', multilink_annot_list).split(' (')[0]
	
if MuXTalk_method == 'MuXTalk_shortest':	

	multilink_edges_entrez_df = return_shortest_multilink_edges_forStreamlit(proj_path, sp_threshold, p1, p2, multilink_type, 
																		st.session_state['KEGG_PPI_allnodes_entrez_df'], 
																		st.session_state['KEGG_path_nodes_entrez_dict'], 
																		st.session_state['A_GRN_sparr'], 
																		st.session_state['A_KEGGPPI_sparr'], 
																		st.session_state['KEGG_interaction_types_dict'], 
																		st.session_state['A_KEGG_e_sparr_dict'])		
elif MuXTalk_method == 'MuXTalk_between':	

	multilink_edges_entrez_df = return_between_multilink_edges(p1, p2, multilink_type, 
																		st.session_state['KEGG_PPI_allnodes_entrez_df'], 
																		st.session_state['KEGG_path_nodes_entrez_dict'], 
																		st.session_state['A_GRN_sparr'], 
																		st.session_state['A_KEGGPPI_sparr'], 
																		st.session_state['KEGG_interaction_types_dict'], 
																		st.session_state['A_KEGG_e_sparr_dict'])		
	
st.write(multilink_edges_entrez_df[['Node 1 Gene Symbol', 'Node 2 Gene Symbol', 'Node 1 Entrez ID', 'Node 2 Entrez ID']])

show_edge_labels = st.checkbox('Show edge labels')

subG_PPI_p1 = nx.Graph(nx.subgraph(st.session_state['G_PPI_entrez'], st.session_state['KEGG_path_nodes_entrez_dict'][p1]))
subG_PPI_p2 = nx.Graph(nx.subgraph(st.session_state['G_PPI_entrez'], st.session_state['KEGG_path_nodes_entrez_dict'][p2]))

KEGG_all_edges_entrez = KEGG_all_edges_entrez[~pd.isnull(KEGG_all_edges_entrez['Edge_subtype'])]
subG_KEGG_p1 = nx.DiGraph()
for e1, e2, e in KEGG_all_edges_entrez[KEGG_all_edges_entrez['Path_label']==p1][['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y', 'Edge_subtype']].values:
	if show_edge_labels:
		subG_KEGG_p1.add_edge(e1, e2, label=e)
	else:
		subG_KEGG_p1.add_edge(e1, e2, label='')

subG_KEGG_p2 = nx.DiGraph()
for e1, e2, e in KEGG_all_edges_entrez[KEGG_all_edges_entrez['Path_label']==p2][['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y', 'Edge_subtype']].values:
	if show_edge_labels:
		subG_KEGG_p2.add_edge(e1, e2, label=e) 
	else:
		subG_KEGG_p2.add_edge(e1, e2, label='')
		
nt = Network('500px', '800px', directed=True, bgcolor='k', font_color='m')
master_G = nx.DiGraph()

if MuXTalk_method == 'MuXTalk_between':	

	for node in st.session_state['KEGG_path_nodes_entrez_dict'][p1]:
		master_G.add_node(node, color='tomato', size=10)
	for e1, e2, in list(subG_PPI_p1.edges()):
		if (e1, e2) not in subG_KEGG_p1.edges():
			if show_edge_labels:
				master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
				master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
			else:
				master_G.add_edge(e1, e2, color='lightgrey', label='')
				master_G.add_edge(e2, e1, color='lightgrey', label='')			
	for e1, e2, in list(subG_KEGG_p1.edges()):
		if show_edge_labels:
			master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p1[e1][e2]['label'])
		else:
			master_G.add_edge(e1, e2, color='grey', label='')
			
	for node in st.session_state['KEGG_path_nodes_entrez_dict'][p2]:
		master_G.add_node(node, color='dodgerblue', size=10)
	for e1, e2, in list(subG_PPI_p2.edges()):
		if (e1, e2) not in subG_KEGG_p2.edges():
			if show_edge_labels:
				master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
				master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
			else:
				master_G.add_edge(e1, e2, color='lightgrey', label='')
				master_G.add_edge(e2, e1, color='lightgrey', label='')			
			
	for e1, e2, in list(subG_KEGG_p2.edges()):
		if show_edge_labels:
			master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p2[e1][e2]['label'])
		else:
			master_G.add_edge(e1, e2, color='grey', label='')
			
	for node in set(st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & set(st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
		master_G.add_node(node, color='gold', size=10)	

	for i, j, ie, je in multilink_edges_entrez_df[['Node 1 Gene Symbol', 'Node 2 Gene Symbol', 'Node 1 Entrez ID', 'Node 2 Entrez ID']].values:
		if '-' not in multilink_type:
			if show_edge_labels:
				master_G.add_edge(ie, je, color='magenta', value=40, label=multilink_type)
			else:
				master_G.add_edge(ie, je, color='magenta', value=40, label='')
		else:
			if show_edge_labels:
				master_G.add_edge(je, ie, color='magenta', value=40, label=multilink_type)
			else:
				master_G.add_edge(je, ie, color='magenta', value=40, label='')
			
elif MuXTalk_method == 'MuXTalk_shortest':

	if sp_threshold == 1:

		for node in st.session_state['KEGG_path_nodes_entrez_dict'][p1]:
			master_G.add_node(node, color='tomato', size=10)
		for e1, e2, in list(subG_PPI_p1.edges()):
			if (e1, e2) not in subG_KEGG_p1.edges():
				if show_edge_labels:
					master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
					master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
				else:
					master_G.add_edge(e1, e2, color='lightgrey', label='')
					master_G.add_edge(e2, e1, color='lightgrey', label='')
									
		for e1, e2, in list(subG_KEGG_p1.edges()):
			if show_edge_labels:
				master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p1[e1][e2]['label'])
			else:
				master_G.add_edge(e1, e2, color='grey', label='')
		
		for node in st.session_state['KEGG_path_nodes_entrez_dict'][p2]:
			master_G.add_node(node, color='dodgerblue', size=10)
		for e1, e2, in list(subG_PPI_p2.edges()):
			if (e1, e2) not in subG_KEGG_p2.edges():
				if show_edge_labels:
					master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
					master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
				else:
					master_G.add_edge(e1, e2, color='lightgrey', label='')
					master_G.add_edge(e2, e1, color='lightgrey', label='')				
				
		for e1, e2, in list(subG_KEGG_p2.edges()):
			if show_edge_labels:
				master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p2[e1][e2]['label'])
			else:
				master_G.add_edge(e1, e2, color='grey', label='')
		
		for node in set(st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & set(st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
			master_G.add_node(node, color='gold', size=10)	
			
		for i, j, ie, je in multilink_edges_entrez_df[['Node 1 Gene Symbol', 'Node 2 Gene Symbol', 'Node 1 Entrez ID', 'Node 2 Entrez ID']].values:
			if (ie not in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & (ie not in st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
				master_G.add_node(ie, size=10, color='darkgrey')
			if (je not in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & (je not in st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
				master_G.add_node(je, size=10, color='darkgrey')

		
			if ((ie in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) == False) & ((je in st.session_state['KEGG_path_nodes_entrez_dict'][p2]) == True):
				extension_edges_df = st.session_state['KEGG_PPI_allnodes_entrez_df'][st.session_state['KEGG_PPI_allnodes_entrez_df'][0].isin(set(nx.neighbors(st.session_state['G_KEGGPPI_entrez'], ie)) & 
																			(st.session_state['KEGG_path_nodes_entrez_dict'][p1] - st.session_state['KEGG_path_nodes_entrez_dict'][p2]))]
				for ii in extension_edges_df[0]:
					if (ii not in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & (i not in st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
						master_G.add_node(ii, color='plum')
					master_G.add_edge(ii, ie, color='plum', value=40)

			elif ((je in st.session_state['KEGG_path_nodes_entrez_dict'][p2]) == False) & ((ie in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) == True):
				extension_edges_df = st.session_state['KEGG_PPI_allnodes_entrez_df'][st.session_state['KEGG_PPI_allnodes_entrez_df'][0].isin(set(nx.neighbors(st.session_state['G_KEGGPPI_entrez'], je)) & 
																			(st.session_state['KEGG_path_nodes_entrez_dict'][p2] - st.session_state['KEGG_path_nodes_entrez_dict'][p1]))]
				for jj in extension_edges_df[0]:
					if (jj not in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & (i not in st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
						master_G.add_node(jj, color='plum')
					master_G.add_edge(je, jj, color='plum', value=40)
								
			if '-' not in multilink_type:
				if show_edge_labels:
					master_G.add_edge(ie, je, color='magenta', value=40, label=multilink_type)
				else:
					master_G.add_edge(ie, je, color='magenta', value=40, label='')
		
			else:
				if show_edge_labels:
					master_G.add_edge(je, ie, color='magenta', value=40, label=multilink_type)	
				else:
					master_G.add_edge(je, ie, color='magenta', value=40, label='')
						
	else:

		for node in st.session_state['KEGG_path_nodes_entrez_dict'][p1]:
			master_G.add_node(node, color='tomato', size=10)
		for e1, e2, in list(subG_PPI_p1.edges()):
			if (e1, e2) not in subG_KEGG_p1.edges():
				if show_edge_labels:
					master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
					master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
				else:
					master_G.add_edge(e1, e2, color='lightgrey', label='')
					master_G.add_edge(e2, e1, color='lightgrey', label='')	
								
		for e1, e2, in list(subG_KEGG_p1.edges()):
			if show_edge_labels:
				master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p1[e1][e2]['label'])
			else:
				master_G.add_edge(e1, e2, color='grey', label='')
				
		for node in st.session_state['KEGG_path_nodes_entrez_dict'][p2]:
			master_G.add_node(node, color='dodgerblue', size=10)
		for e1, e2, in list(subG_PPI_p2.edges()):
			if (e1, e2) not in subG_KEGG_p2.edges():
				if show_edge_labels:
					master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
					master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
				else:
					master_G.add_edge(e1, e2, color='lightgrey', label='')
					master_G.add_edge(e2, e1, color='lightgrey', label='')	
								
		for e1, e2, in list(subG_KEGG_p2.edges()):
			if show_edge_labels:
				master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p2[e1][e2]['label'])
			else:
				master_G.add_edge(e1, e2, color='grey', label='')
		
		for node in set(st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & set(st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
			master_G.add_node(node, color='gold', size=10)	
			
		for i, j, ie, je in multilink_edges_entrez_df[['Node 1 Gene Symbol', 'Node 2 Gene Symbol', 'Node 1 Entrez ID', 'Node 2 Entrez ID']].values:
			if (ie not in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & (ie not in st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
				master_G.add_node(ie, size=10, color='darkgrey')
			if (je not in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & (je not in st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
				master_G.add_node(je, size=10, color='darkgrey')
					
			if ((ie in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) == False) & ((je in st.session_state['KEGG_path_nodes_entrez_dict'][p2]) == True):
			
				pp1_ie_sp = []
				for pp1 in (st.session_state['KEGG_path_nodes_entrez_dict'][p1] - st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
					if nx.has_path(G_KEGGPPI_entrez, pp1, ie):
						pp1_ie_sp.append(nx.shortest_path(G_KEGGPPI_entrez, pp1, ie))
				extension_edges = set([(e1, e2)  for temp_pp1_ie_sp in pp1_ie_sp for e1, e2 in zip(temp_pp1_ie_sp, temp_pp1_ie_sp[1:])])

			elif ((je in st.session_state['KEGG_path_nodes_entrez_dict'][p2]) == False) & ((ie in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) == True):
			
				je_pp2_sp = []
				for pp2 in (st.session_state['KEGG_path_nodes_entrez_dict'][p2] - st.session_state['KEGG_path_nodes_entrez_dict'][p1]):
					if nx.has_path(G_KEGGPPI_entrez, je, pp2):
						je_pp2_sp.append(nx.shortest_path(G_KEGGPPI_entrez, je, pp2))
				extension_edges = set([(e1, e2)  for temp_je_pp2_sp in je_pp2_sp for e1, e2 in zip(temp_je_pp2_sp, temp_je_pp2_sp[1:])])

			elif ((ie in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) == False) & ((je in st.session_state['KEGG_path_nodes_entrez_dict'][p2]) == False):
		
				pp1_ie_sp = []
				for pp1 in (st.session_state['KEGG_path_nodes_entrez_dict'][p1] - st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
					if nx.has_path(G_KEGGPPI_entrez, pp1, ie):
						pp1_ie_sp.append(nx.shortest_path(G_KEGGPPI_entrez, pp1, ie))
				extension_edges_pp1_ie = set([(e1, e2)  for temp_pp1_ie_sp in pp1_ie_sp for e1, e2 in zip(temp_pp1_ie_sp, temp_pp1_ie_sp[1:])])	

				je_pp2_sp = []
				for pp2 in (st.session_state['KEGG_path_nodes_entrez_dict'][p2] - st.session_state['KEGG_path_nodes_entrez_dict'][p1]):
					if nx.has_path(G_KEGGPPI_entrez, je, pp2):
						je_pp2_sp.append(nx.shortest_path(G_KEGGPPI_entrez, je, pp2))
				extension_edges_je_pp2 = set([(e1, e2)  for temp_je_pp2_sp in je_pp2_sp for e1, e2 in zip(temp_je_pp2_sp, temp_je_pp2_sp[1:])])			
			
				extension_edges = extension_edges_pp1_ie | extension_edges_je_pp2
				
			for e1, e2 in extension_edges:
				if (e1 not in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & (e1 not in st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
					master_G.add_node(e1, size=5, color='lightgrey')
				if (e2 not in st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & (e2 not in st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
					master_G.add_node(e2, size=5, color='lightgrey')
				master_G.add_edge(e1, e2, color='plum', value=20)
				
			if '-' not in multilink_type:
				if show_edge_labels:
					master_G.add_edge(ie, je, color='magenta', value=40, label=multilink_type)
				else:
					master_G.add_edge(ie, je, color='magenta', value=40, label='')				
			else:
				if show_edge_labels:
					master_G.add_edge(je, ie, color='magenta', value=40, label=multilink_type)
				else:
					master_G.add_edge(je, ie, color='magenta', value=40, label='')

master_G = nx.relabel_nodes(master_G, st.session_state['KEGG_PPI_allnodes_entrez_df'].set_index(0).to_dict()['Approved symbol'])
nt.from_nx(master_G)					
nt.show_buttons(filter_=["physics"])

nt.show('muxtalk_nx.html')
f_html = open('muxtalk_nx.html', 'r', encoding='utf-8')
source_html =  f_html.read()
components.html(source_html, height=900, width=900)
	
	
tips_expander = st.expander(label='Tips for adjusting network visulaziation parameters')
with tips_expander:
	tips_msg = st.container()	
	tips_msg.write('''
	- For fast convergence and visual clarity, we recommend to use the solver "forceAtlas2Based."
	- To be able to move nodes freely, check off the "enable" box.
	- To "cool down" the networks, i.e. keep the nodes from moving, you can decrease gravitational constant (e.g. to ~ -20,000)\
	and increase damping (e.g. to ~ 0.3). 
	''')			
																					
						