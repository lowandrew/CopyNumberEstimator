#!/usr/bin/env python

import argparse
import requests
import logging
import json
import os


def get_list_of_orthodb_ids():
    # This query gets a list of all genes in bacteria that are present in at least 90 percent of species and single copy
    # in at least 90 percent of species.
    query = 'https://www.orthodb.org//search?query=&level=2&species=2&universal=0.9&singlecopy=0.9'
    response = requests.get(query)
    if response.status_code == '404':
        raise Exception('OrthoDB must have changed something! Query no longer works.')
    response_dict = json.loads(response.content)
    return response_dict['data']


def download_ortho_db_fasta(ortho_db_id):
    query = 'https://www.orthodb.org/fasta?id={}&species=2'.format(ortho_db_id)
    response = requests.get(query)
    return response.content


if __name__ == '__main__':
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                    level=logging.INFO,
                    datefmt='%Y-%m-%d %H:%M:%S')
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_directory',
                        type=str,
                        required=True,  # TODO: Make this default to an env var or something.
                        help='Folder to store universal(ish) single copy(ish) genes. Must not already exist.')
    args = parser.parse_args()
    if os.path.isdir(args.output_directory):
        raise FileExistsError('Output folder specified ({}) already exists. Please choose a folder that does not '
                              'already exist.'.format(args.output_directory))
    os.makedirs(args.output_directory)
    ortho_db_ids = get_list_of_orthodb_ids()
    logging.info('Found {} genes to download.'.format(len(ortho_db_ids)))
    downloaded_count = 1
    for odbid in ortho_db_ids:
        logging.info('Downloading gene {} of {}'.format(downloaded_count, len(ortho_db_ids)))
        fasta_content = download_ortho_db_fasta(odbid)
        with open(os.path.join(args.output_directory, odbid + '.fasta'), 'w') as f:
            f.write(fasta_content)
        downloaded_count += 1
