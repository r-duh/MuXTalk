import streamlit as st
import streamlit.components.v1 as components
from MuXTalk_functions_Streamlit_Docker import *
import networkx as nx
from pyvis.network import Network

st.set_page_config(layout='wide')

st.sidebar.title(
	'''
	MuXTalk: Detecting and dissecting signaling crosstalk via multilayer networks
	***
	
	'''
)

input_GRN = st.sidebar.selectbox('Select GRN', ['HumanGRN10e6', 'HumanGRN10e5', 'HumanGRN10e4'])	
MuXTalk_method = st.sidebar.selectbox('Select MuXTalk type', ['MuXTalk_shortest', 'MuXTalk_between'])
sp_threshold = st.sidebar.selectbox('Select shortest path threshold (MuXTalk_shortest only)', [2, 1, 'None'])
custom_GRN_edges = st.sidebar.file_uploader('Custom GRN edges')
custom_GRN_parquet = st.sidebar.file_uploader('Custom GRN parquet file')

with st.container():
	'''
	# **MuXTalk:** Detecting and dissecting signaling crosstalk via multilayer networks
	
	***
	
	'''

proj_path = '/Volumes/Partition1/Pathway_crosstalk_files/MuXTalk_for_Docker_final/'
input_filenames_dict = {'HUGO_symb_entrez_uniprot': 'HugoGene_20200528.txt', 'PPI_Cheng_2019_data': 'PPI_Cheng_NatComms2019.csv',
					   'KEGG_all_nodes_df': 'KEGG_expanded_all_nodes.csv', 'KEGG_all_edges_df': 'KEGG_expanded_all_edges.csv',
						'df_motinf' : 'cisbpall_motinf.txt', 'XTalk_DB': 'XTalkDB_crosstalk.csv', 
						'almen_etal': 'almen_etal_12915_2009_258_MOESM1_ESM.csv', 
						'lambert_etal': 'lambert_etal_1-s2.0-S0092867418301065-mmc2_TableS1.csv'
					   }
					   
expander_readme = st.expander(label='Read me first!', expanded=False)
with expander_readme:
	readme_msg = st.container()	

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

	G_PPI =  nx.from_scipy_sparse_matrix(A_PPI_sparr.astype(bool).astype(int))
	G_PPI_entrez = nx.relabel_nodes(G_PPI, KEGG_PPI_allnodes_entrez_df[[0, 'ix']].to_dict()[0])
	G_KEGGPPI = nx.from_scipy_sparse_matrix(A_KEGGPPI_sparr.astype(bool).astype(int))
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

	G_PPI =  nx.from_scipy_sparse_matrix(A_PPI_sparr.astype(bool).astype(int))
	G_PPI_entrez = nx.relabel_nodes(G_PPI, KEGG_PPI_allnodes_entrez_df[[0, 'ix']].to_dict()[0])
	G_KEGGPPI = nx.from_scipy_sparse_matrix(A_KEGGPPI_sparr.astype(bool).astype(int))
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
		
form_cols_selectmultilink = st.columns(2)
multilink_type = form_cols_selectmultilink[0].selectbox('Select multilink type', st.session_state['detected_ranked_pathway_pairs_df'].loc['%s-%s' % (p1, p2)]['sig multilink types'])
	
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

subG_PPI_p1 = nx.Graph(nx.subgraph(st.session_state['G_PPI_entrez'], st.session_state['KEGG_path_nodes_entrez_dict'][p1]))
subG_PPI_p2 = nx.Graph(nx.subgraph(st.session_state['G_PPI_entrez'], st.session_state['KEGG_path_nodes_entrez_dict'][p2]))
subG_KEGG_p1 = nx.DiGraph()
for e1, e2, e in KEGG_all_edges_entrez[KEGG_all_edges_entrez['Path_label']==p1][['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y', 'Edge_subtype']].values:
	subG_KEGG_p1.add_edge(e1, e2, label=e) 
subG_KEGG_p2 = nx.DiGraph()
for e1, e2, e in KEGG_all_edges_entrez[KEGG_all_edges_entrez['Path_label']==p2][['NCBI Gene ID(supplied by NCBI)_x', 'NCBI Gene ID(supplied by NCBI)_y', 'Edge_subtype']].values:
	subG_KEGG_p2.add_edge(e1, e2, label=e) 


nt = Network('500px', '800px', directed=True, bgcolor='k', font_color='m')
master_G = nx.DiGraph()

if MuXTalk_method == 'MuXTalk_between':	

	for node in st.session_state['KEGG_path_nodes_entrez_dict'][p1]:
		master_G.add_node(node, color='tomato', size=10)
	for e1, e2, in list(subG_PPI_p1.edges()):
		if (e1, e2) not in subG_KEGG_p1.edges():
			master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
			master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
	for e1, e2, in list(subG_KEGG_p1.edges()):
		master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p1[e1][e2]['label'])
		
	for node in st.session_state['KEGG_path_nodes_entrez_dict'][p2]:
		master_G.add_node(node, color='dodgerblue', size=10)
	for e1, e2, in list(subG_PPI_p2.edges()):
		if (e1, e2) not in subG_KEGG_p1.edges():
			master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
			master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
	for e1, e2, in list(subG_KEGG_p2.edges()):
		master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p2[e1][e2]['label'])
		
	for node in set(st.session_state['KEGG_path_nodes_entrez_dict'][p1]) & set(st.session_state['KEGG_path_nodes_entrez_dict'][p2]):
		master_G.add_node(node, color='gold', size=10)	

	for i, j, ie, je in multilink_edges_entrez_df[['Node 1 Gene Symbol', 'Node 2 Gene Symbol', 'Node 1 Entrez ID', 'Node 2 Entrez ID']].values:
		if '-' not in multilink_type:
			master_G.add_edge(ie, je, color='magenta', value=40, label=multilink_type)
		else:
			master_G.add_edge(je, ie, color='magenta', value=40, label=multilink_type)
			
elif MuXTalk_method == 'MuXTalk_shortest':

	if sp_threshold == 1:

		for node in st.session_state['KEGG_path_nodes_entrez_dict'][p1]:
			master_G.add_node(node, color='tomato', size=10)
		for e1, e2, in list(subG_PPI_p1.edges()):
			if (e1, e2) not in subG_KEGG_p1.edges():
				master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
				master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
		for e1, e2, in list(subG_KEGG_p1.edges()):
			master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p1[e1][e2]['label'])
		
		for node in st.session_state['KEGG_path_nodes_entrez_dict'][p2]:
			master_G.add_node(node, color='dodgerblue', size=10)
		for e1, e2, in list(subG_PPI_p2.edges()):
			if (e1, e2) not in subG_KEGG_p1.edges():
				master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
				master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
		for e1, e2, in list(subG_KEGG_p2.edges()):
			master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p2[e1][e2]['label'])
		
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
				master_G.add_edge(ie, je, color='magenta', value=40, label=multilink_type)
			else:
				master_G.add_edge(je, ie, color='magenta', value=40, label=multilink_type)					
						
	else:

		for node in st.session_state['KEGG_path_nodes_entrez_dict'][p1]:
			master_G.add_node(node, color='tomato', size=10)
		for e1, e2, in list(subG_PPI_p1.edges()):
			if (e1, e2) not in subG_KEGG_p1.edges():
				master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
				master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
		for e1, e2, in list(subG_KEGG_p1.edges()):
			master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p1[e1][e2]['label'])
		
		for node in st.session_state['KEGG_path_nodes_entrez_dict'][p2]:
			master_G.add_node(node, color='dodgerblue', size=10)
		for e1, e2, in list(subG_PPI_p2.edges()):
			if (e1, e2) not in subG_KEGG_p1.edges():
				master_G.add_edge(e1, e2, color='lightgrey', label='PPI')
				master_G.add_edge(e2, e1, color='lightgrey', label='PPI')
		for e1, e2, in list(subG_KEGG_p2.edges()):
			master_G.add_edge(e1, e2, color='grey', label=subG_KEGG_p2[e1][e2]['label'])
		
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
				master_G.add_edge(ie, je, color='magenta', value=40, label=multilink_type)
			else:
				master_G.add_edge(je, ie, color='magenta', value=40, label=multilink_type)

master_G = nx.relabel_nodes(master_G, st.session_state['KEGG_PPI_allnodes_entrez_df'].set_index(0).to_dict()['Approved symbol'])
nt.from_nx(master_G)					
nt.show_buttons(filter_=["physics"])

nt.show('nx.html')
f_html = open('nx.html', 'r', encoding='utf-8')
source_html =  f_html.read()
components.html(source_html, height=900, width=900)
	
	
expander1 = st.expander(label='Tips for adjusting network physics parameters')
with expander1:
	':bulb: To "cool down" the networks, i.e. keep the nodes from moving, you can decrease gravitational constant (e.g. to ~ -20,000)\
	and increase damping (e.g. to ~ 0.3). To stop physics from being implemented altogether and to be able to move nodes freely\
	simply check off the "enable" box.'			
																					
			
		
			