import pandas as pd
import pickle
import os
import matplotlib.pyplot as plt
import networkx as nx

from collections import defaultdict

def get_dtis(dti_path, target_col='UniProtID', drug_col='DrugBankID',
             aff_score_col='pIC50', smiles_dict=None):
    dti_df = pd.read_csv(dti_path, usecols=(target_col, drug_col, aff_score_col))
    targets = dti_df[target_col].unique().tolist()
    drugs = dti_df[drug_col].unique().tolist()
    if target_col == 'protein name':
        targets = ['\n'.join(target.split()) for target in targets]
    if drug_col == 'SMILES':
        drugs = [smiles_dict[smiles] for smiles in drugs]
    
    dtis = defaultdict(list)
    for i in dti_df.index:
        tid = dti_df.loc[i, target_col]
        did = dti_df.loc[i, drug_col]
        weight = round(dti_df.loc[i, aff_score_col], 4)
        if target_col == 'protein name':
            tid = '\n'.join(tid.split())
        if drug_col == 'SMILES':
            did = smiles_dict[did]
            
        dtis[did].append((tid, weight))
        
    return dtis, targets, drugs

def get_ddis(ddi_path, drug, smiles2dname=False, smiles_dict=None):
    with open(ddi_path, 'rb') as f:
        ddis_dict = pickle.load(f)
    
    if ddis_dict.get(drug):
        sim_drugs_with_sim = ddis_dict[drug]
        sim_drugs = [sim_drug for sim_drug, _ in sim_drugs_with_sim]
        if smiles2dname:
            with open('data/dbid2dname_dict.pkl', 'rb') as f:
                dbid2dname_dict = pickle.load(f)
                
            drug = dbid2dname_dict[drug]
            sim_drugs_with_sim = [(dbid2dname_dict[sim_drug_with_sim], sim)
                                   for sim_drug_with_sim, sim in sim_drugs_with_sim]
            sim_drugs = [dbid2dname_dict[sim_drug]
                         for sim_drug in sim_drugs]
            
        ddis = [(drug, sim_drug, round(sim, 4)) for sim_drug, sim in sim_drugs_with_sim]
        return ddis, sim_drugs
    return [], []

def get_prev_dtis(ddis, prev_dtis_path, prev_cov_path, smiles2dname=False, prev_cov=True):
    with open(prev_dtis_path, 'rb') as f:
        dtis_dict = pickle.load(f)
        
    if smiles2dname:
        with open('data/dname2dbid_dict.pkl', 'rb') as f:
            dname2dbid_dict = pickle.load(f)
            
        with open('data/dbid2dname_dict.pkl', 'rb') as f:
            dbid2dname_dict = pickle.load(f)

    prev_cov_df = pd.read_csv(prev_cov_path)
    prev_covs = prev_cov_df['prev_cov'].tolist()
    
    prev_dtis = []
    for _, sim_drug, _ in ddis:
        if smiles2dname:
            sim_drug = dname2dbid_dict[sim_drug]
            
        if dtis_dict.get(sim_drug):
            prev_targets = dtis_dict[sim_drug]
            if smiles2dname:
                sim_drug = dbid2dname_dict[sim_drug]
                
            if prev_targets:
                for prev_target in prev_targets:
                    if prev_target in prev_covs or not prev_cov:
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

def create_and_save_dti_graph(targets, drug, sim_drugs, prev_sars_targets,
                              dtis, ddis, prev_dtis, output_dir):
    G = multilayered_graph(targets, drug, sim_drugs, prev_sars_targets, dtis, ddis, prev_dtis)
    # color = [subset_color[data["layer"]] for v, data in G.nodes(data=True)]
    pos = nx.multipartite_layout(G, subset_key="layer")
    plt.figure(figsize=(10, 8))
    # nx.draw(G, pos, node_color=color, with_labels=False)
    nx.draw(G, pos, with_labels=False)

    t_options = {'node_size': 300, 'node_color': 'aquamarine', 'node_shape': 's'}
    nx.draw_networkx_nodes(G, pos, nodelist=targets, **t_options)
    d_options = {'node_size': 300, 'node_color': 'orangered', 'node_shape': 'H'}
    nx.draw_networkx_nodes(G, pos, nodelist=[drug], **d_options)
    sd_options = {'node_size': 300, 'node_color': 'lightsalmon', 'node_shape': 'H'}
    nx.draw_networkx_nodes(G, pos, nodelist=sim_drugs, **sd_options)
    pt_options = {'node_size': 300, 'node_color': 'gold', 'node_shape': 's'}
    nx.draw_networkx_nodes(G, pos, nodelist=prev_sars_targets, **pt_options)

    labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)

    for p in pos:  # raise text positions
        pos[p][1] += 0.065
    nx.draw_networkx_labels(G, pos)

    axis = plt.gca()
    axis.set_xlim([1.5*x for x in axis.get_xlim()])
    axis.set_ylim([1.5*y for y in axis.get_ylim()])
    
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    
    graph_fname = f'dti_graph_for_{drug}.png'
    plt.savefig(f'{output_dir}/{graph_fname}')
    print(f'{graph_fname} is saved.')
    plt.show()
    
    return graph_fname
    
def draw_dti_graph(dti_path, drug,
                   target_col, drug_col, aff_score_col,
                   ddi_path='data/ddis.pkl',
                   prev_dtis_path='data/dtis.pkl',
                   prev_cov_path='data/prev_cov.csv',
                   smiles2dname=False, prev_cov=True,
                   output_dir='dti_graphs'):
    '''Draw DTI graph'''
    
    db_id, smiles_dict = drug, None
    if drug_col == 'SMILES':
        if smiles2dname:
            with open('data/smiles2dbid_dict.pkl', 'rb') as f:
                smiles2dbid_dict = pickle.load(f)
            db_id = smiles2dbid_dict[drug]
            
            with open('data/smiles2dname_dict.pkl', 'rb') as f:
                smiles_dict = pickle.load(f)
            drug = smiles_dict[drug]
        else:
            with open('data/smiles2dbid_dict.pkl', 'rb') as f:
                smiles_dict = pickle.load(f)
            db_id = smiles_dict[drug]
            drug = db_id
    
    dtis, targets, drugs = get_dtis(dti_path, target_col,
                                    drug_col, aff_score_col,
                                    smiles_dict=smiles_dict)
    
    if drug not in drugs:
        print(f'{drug} is not exist in predicted DTIs.')
        return

        
    ddis, sim_drugs = get_ddis(ddi_path, db_id, smiles2dname=smiles2dname, smiles_dict=smiles_dict)
    prev_dtis, prev_sars_targets = get_prev_dtis(ddis, prev_dtis_path, prev_cov_path,
                                                 smiles2dname=smiles2dname, prev_cov=prev_cov)

    print(f'# of targets: {len(targets)}, # of drugs: {len(drugs)}')
    print(f'DTI sample: {dtis[list(dtis.keys())[0]]}')

    graph_fname = create_and_save_dti_graph(targets, drug, sim_drugs, prev_sars_targets,
                                            dtis, ddis, prev_dtis, output_dir)
    return graph_fname