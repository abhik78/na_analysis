import csv
import requests
import argparse
import os

def open_text_file(filename):
    '''
    read csv file and yiled rows
    :param decoy_filename: read .picked files
    :return: return rows as generator
    '''
    with open(filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            yield row


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv_file', '-f', help='csv file of pdb entries')
    parser.add_argument('-o', '--output_directory', type=str, default='screen_data',
                        help='Output directory')
    args = parser.parse_args()
    return args


def fetch_pdb(entry_id):
    '''
    fetch pdb file from rcsb
    :param entry_id: entry id
    :return: pdb file
    '''

    args = parse_arguments()
    output_dir = args.output_directory
    req = requests.get('http://files.rcsb.org/download/{}.pdb'.format(entry_id), '{}.pdb'.format(entry_id))

    filename = '{}.pdb'.format(entry_id)
    with open(os.path.join(output_dir, filename), 'wb') as f:
        f.write(req.content)




def main():
    args = parse_arguments()

    pdb_id_list = [row[0] for row in open_text_file(args.csv_file)]

    for entry in pdb_id_list:
        fetch_pdb(entry)

if __name__=='__main__':
     main()