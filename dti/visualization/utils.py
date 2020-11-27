import pandas as pd
import pickle
import matplotlib.pyplot as plt
import networkx as nx

from collections import defaultdict

def get_dtis(dti_path):
    dti_df = pd.read_csv(dti_path, usecols=('UniProtID', 'DrugBankID', 'pIC50'))
    targets = dti_df['UniProtID'].unique().tolist()
    drugs = dti_df['DrugBankID'].unique().tolist()
    
    dtis = defaultdict(list)
    for i in dti_df.index:
        tid = dti_df.loc[i, 'UniProtID']
        did = dti_df.loc[i, 'DrugBankID']
        weight = dti_df.loc[i, 'pIC50']
        dtis[did].append((tid, weight))
        
    return dtis, targets, drugs

def get_ddis(ddi_path, drug):
    with open(ddi_path, 'rb') as f:
        ddis_dict = pickle.load(f)
    
    sim_drugs_with_sim = ddis_dict[drug]
    sim_drugs = [sim_drug for sim_drug, _ in sim_drugs_with_sim]
    ddis = [(drug, sim_drug, round(sim, 4)) for sim_drug, sim in sim_drugs_with_sim]
    return ddis, sim_drugs

def get_prev_dtis(ddis, prev_dtis_path, prev_cov_path):
    with open(prev_dtis_path, 'rb') as f:
        dtis_dict = pickle.load(f)

    prev_cov_df = pd.read_csv(prev_cov_path)
    prev_cov = prev_cov_df['prev_cov'].tolist()
    
    prev_dtis = []
    for _, sim_drug, _ in ddis:
        prev_targets = dtis_dict[sim_drug]
        if prev_targets:
            for prev_target in prev_targets:
                if prev_target in prev_cov:
                    prev_dtis.append((prev_target, sim_drug, 1))
    
    prev_sars_targets = [t for t, _, _ in prev_dtis]
    return prev_dtis, prev_sars_targets

def multilayered_graph(targets, drug, sim_drugs, prev_sars_targets, dtis, ddis, prev_dtis):
    G = nx.Graph()
    G.add_nodes_from(targets, layer=0)
    G.add_nodes_from([drug], layer=1)
    G.add_nodes_from(sim_drugs, layer=2)
    G.add_nodes_from(prev_sars_targets, layer=3)
    
    for drug_in_dti, target_list in dtis.items():
        for target, weight in target_list:
            if drug_in_dti == drug:
                G.add_edge(target, drug_in_dti, weight=weight)

    for d1, d2, weight in ddis:
        G.add_edge(d1, d2, weight=weight)

    for t, d, weight in prev_dtis:
        G.add_edge(t, d)
        
    return G

def create_and_save_dti_graph(targets, drug, sim_drugs, prev_sars_targets, dtis, ddis, prev_dtis):
    G = multilayered_graph(targets, drug, sim_drugs, prev_sars_targets, dtis, ddis, prev_dtis)
    # color = [subset_color[data["layer"]] for v, data in G.nodes(data=True)]
    pos = nx.multipartite_layout(G, subset_key="layer")
    plt.figure(figsize=(8, 8))
    # nx.draw(G, pos, node_color=color, with_labels=False)
    nx.draw(G, pos, with_labels=False)

    t_options = {'node_size': 300, 'node_color': 'g', 'node_shape': 's'}
    nx.draw_networkx_nodes(G, pos, nodelist=targets, **t_options)
    d_options = {'node_size': 300, 'node_color': 'r', 'node_shape': 'H'}
    nx.draw_networkx_nodes(G, pos, nodelist=[drug], **d_options)
    pt_options = {'node_size': 300, 'node_color': 'gold', 'node_shape': 's'}
    nx.draw_networkx_nodes(G, pos, nodelist=prev_sars_targets, **pt_options)


    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)

    for p in pos:  # raise text positions
        pos[p][1] += 0.065
    nx.draw_networkx_labels(G, pos)

    # plt.axis("equal")
    plt.savefig('dti_graph.png')
    plt.show()