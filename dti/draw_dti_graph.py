import argparse
from visualization.utils import get_dtis, get_ddis, get_prev_dtis, \
                                create_and_save_dti_graph

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a DTI graph.')
    parser.add_argument('--dti', type=str, default='data/sample_dtis.csv',
                        help='the path which saved predicted DTIs')
    parser.add_argument('--ddi', type=str, default='data/ddis.pkl',
                        help='the path which saved DDIs from DrugBank')
    parser.add_argument('--pdti', type=str, default='data/dtis.pkl',
                        help='the path which saved DTIs from DrugBank')
    parser.add_argument('--pcov', type=str, default='data/prev_cov.csv',
                        help='the path which saved previous corona viruses')

    args = parser.parse_args()
    print(args)
    
    dti_path, ddi_path = args.dti, args.ddi
    prev_dtis_path, prev_cov_path = args.pdti, args.pcov
    
    dtis, targets, drugs = get_dtis(dti_path)
    ddis, sim_drugs = get_ddis(ddi_path, drugs[0])
    prev_dtis, prev_sars_targets = get_prev_dtis(ddis, prev_dtis_path, prev_cov_path)
    
    print(f'# of targets: {len(targets)}, # of drugs: {len(drugs)}')
    print(f'dtis sample: {dtis[list(dtis.keys())[0]]}')
    print(f'ddis sample: {ddis[0]}')
    print(f'prev_dtis sample: {prev_dtis[0]}')
    
    create_and_save_dti_graph(targets, drugs[0], sim_drugs, prev_sars_targets, dtis, ddis, prev_dtis)
    