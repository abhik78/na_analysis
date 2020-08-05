#!/usr/bin/env python
from ccdc import io
from ccdc import conformer
import pandas as pd
import csv
import multiprocessing
import os
import argparse

os.environ['CSDHOME'] = '/home/amukhopadhyay/not-backed-up/CCDC/CSD_2020'

def open_text_file(filename):
    '''
    read csv file and yiled rows
    :param decoy_filename: read .picked files
    :return: return rows as generator
    '''
    with open(filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        #next(csv_reader) #ignore the first line
        for row in csv_reader:
            yield row

def parse_command_line_args():
    """Return the command line arguments.

    Parse the command line arguments and return them."""

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-d', '--pdb_dir', type=str, default='',
                        help='pdb files dir')

    parser.add_argument('--pdb_ids', '-p',
                        help='single column csv file with target uniprot ids, header uniprot_id')
    parser.add_argument('-o', '--output_directory', type=str, default='screen_data',
                        help='Output directory')

    args = parser.parse_args()

    return args


class MogulProcessor:

    def __init__(self, args):
        self.args = args

    def mogul_engine(self, pdb_id):
        engine = conformer.GeometryAnalyser()
        pdb_dir = self.args.pdb_dir
        pdb_file = os.path.join(pdb_dir, '{}.pdb'.format(pdb_id))
        output_dir = self.args.output_directory

        print(pdb_file)
        #p_6wk7 = 'pdb6wk7.ent'
        mol_reader = io.MoleculeReader(pdb_file)

        mol = mol_reader[0]

        mol_reader.close()

        mol.remove_hydrogens()
        mol.assign_bond_types(which='unknown')
        mol.standardise_aromatic_bonds()
        mol.standardise_delocalised_bonds()
        mol.add_hydrogens()


        geometry_analysed_mol = engine.analyse_molecule(mol)

        torsions_df = pd.DataFrame(
                [(x.fragment_label, engine.fragment_identifier(x), x.value, x.unusual, x.enough_hits, x.d_min, x.local_density, x.z_score, x.standard_deviation) for x in geometry_analysed_mol.analysed_torsions],
                columns=['atom_labels', 'fragment_id', 'value', 'unusual', 'enough_hits', 'd_min', 'local_density', 'z_score', 'std_dev']
            ).sort_values('d_min', ascending=False).reset_index(drop=True)
        torsions_df.to_csv(os.path.join(output_dir, 'torsions_{}.csv'.format(pdb_id)), index=False)

        angles_df = pd.DataFrame(
                [(x.fragment_label, engine.fragment_identifier(x), x.value, x.unusual, x.enough_hits, x.d_min, x.local_density, x.z_score, x.standard_deviation) for x in geometry_analysed_mol.analysed_angles],
                columns=['atom_labels', 'fragment_id', 'value', 'unusual', 'enough_hits', 'd_min', 'local_density', 'z_score', 'std_dev']
            ).sort_values('d_min', ascending=False).reset_index(drop=True)
        angles_df.to_csv(os.path.join(output_dir, 'angles_{}.csv'.format(pdb_id)), index=False)

        bonds_df = pd.DataFrame(
                [(x.fragment_label, engine.fragment_identifier(x), x.value, x.unusual, x.enough_hits, x.d_min, x.local_density, x.z_score, x.standard_deviation) for x in geometry_analysed_mol.analysed_bonds],
                columns=['atom_labels', 'fragment_id', 'value', 'unusual', 'enough_hits', 'd_min', 'local_density', 'z_score', 'std_dev']
            ).sort_values('d_min', ascending=False).reset_index(drop=True)
        bonds_df.to_csv(os.path.join(output_dir, 'bonds_{}.csv'.format(pdb_id)), index=False)
        print("done")


def main():
    p = multiprocessing.Pool(8)
    args = parse_command_line_args()

    pdb_id_list = [row[0] for row in open_text_file(args.pdb_ids)]
    print(pdb_id_list)
    mogul_proc = MogulProcessor(args)
    p.map(mogul_proc.mogul_engine, pdb_id_list)

if __name__=='__main__':
     main()



