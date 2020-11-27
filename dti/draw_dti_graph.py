import argparse
import sys
from visualization.utils import get_dtis, get_ddis, get_prev_dtis, \
                                create_and_save_dti_graph

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a DTI graph.')
    parser.add_argument('--dti', type=str, default='data/sample_dtis.csv',
                        help='the path which saved predicted DTIs')
    parser.add_argument('--drug', type=str, default='DB00503',
                        help='the DrugBank ID which you want to see a DTI graph')
    parser.add_argument('--ddi', type=str, default='data/ddis.pkl',
                        help='the path which saved DDIs from DrugBank')
    parser.add_argument('--pdti', type=str, default='data/dtis.pkl',
                        help='the path which saved DTIs from DrugBank')
    parser.add_argument('--pcov', type=str, default='data/prev_cov.csv',
                        help='the path which saved previous corona viruses')

    args = parser.parse_args()
    
    dti_path, ddi_path = args.dti, args.ddi
    prev_dtis_path, prev_cov_path = args.pdti, args.pcov
    drug = args.drug
    
    dtis, targets, drugs = get_dtis(dti_path)
    if drug not in drugs:
        print(f'{drug} is not exist in predicted DTIs.')
        sys.exit()
        
    ddis, sim_drugs = get_ddis(ddi_path, drug)
    prev_dtis, prev_sars_targets = get_prev_dtis(ddis, prev_dtis_path, prev_cov_path)
    
    print(f'# of targets: {len(targets)}, # of drugs: {len(drugs)}')
    print(f'DTI sample: {dtis[list(dtis.keys())[0]]}')
    
    create_and_save_dti_graph(targets, drug, sim_drugs, prev_sars_targets, dtis, ddis, prev_dtis)
    