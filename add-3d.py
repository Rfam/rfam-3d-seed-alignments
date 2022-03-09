#!/usr/bin/env python3
"""
Add secondary structure based on 3D annotations to Rfam seed alignments.

Infernal binaries should be available in PATH.

Usage:
python add-3d.py
python add-3d.py RF00008
"""


import argparse
import os
import re
import subprocess

import collections

from fr3d_2d import fr3d_2d


SKIP_LARGE_ALIGNMENT = 30

PDB_BLACKLIST = [
    '7AS5',  # DNA
    '6LAB',  # DNA
]

FAMILY_BLACKLIST = [
    'RF00029', # group II intron, too large
    'RF00106', # RNAI matches a DNA molecule 7NPN
    'RF00843', # microRNA matching DNA in complex with histones
    'RF02545', # Trypanosomatid SSU
    'RF02546', # Trypanosomatid LSU
]


def download_rfam_files(rfam_acc, nocache):
    """
    Download and uncompress Rfam CM and SEED files for a given family.
    """
    if not os.path.exists('temp'):
        os.system('mkdir temp')
    if not os.path.exists('data/cm/{}.cm'.format(rfam_acc)) or nocache:
        print('Downloading the latest {} CM from SVN'.format(rfam_acc))
        cmd = 'wget -q -O data/cm/{0}.cm  https://xfamsvn.ebi.ac.uk/svn/data_repos/trunk/Families/{0}/CM'.format(rfam_acc)
        subprocess.check_output(cmd, shell=True)
    if not os.path.exists('data/seed/{}.seed'.format(rfam_acc)) or nocache:
        print('Downloading the latest {} SEED from SVN'.format(rfam_acc))
        cmd = 'wget -q -a temp/wget.log -O data/seed/{0}.seed https://xfamsvn.ebi.ac.uk/svn/data_repos/trunk/Families/{0}/SEED'.format(rfam_acc)
        subprocess.check_output(cmd, shell=True)


def get_rfam_3d_mapping():
    """
    Parse Rfam PDB mapping file:

    RF00012	6zqb	D4	2	175	54.50	7.2e-12	1	218	c00f0f	0
    RF00015	2n7m	X	1	74	61.90	7.1e-14	1	139	c00f0f	1
    """
    if not os.path.exists('pdb_full_region.txt'):
        # cmd = 'wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/.preview/pdb_full_region.txt.gz && gunzip pdb_full_region.txt.gz'
        cmd = 'wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/pdb_full_region.txt.gz && gunzip pdb_full_region.txt.gz'
        subprocess.check_output(cmd, shell=True)
    data = collections.defaultdict(list)
    with open('pdb_full_region.txt', 'r') as f_pdb:
        for line in f_pdb:
            if line.startswith('rfam_acc'):
                continue  # skip header
            parts = re.split(r'\s+', line)
            if parts[10] == '0':  # skip redundant clan entries
                continue
            pdb_id = '{}_{}'.format(parts[1].upper(), parts[2])
            data[parts[0]].append(pdb_id)
    # remove duplicates (e.g. palindromic miRNAs matching in 2 directions)
    for rfam_acc, pdbs in data.items():
        data[rfam_acc] = list(dict.fromkeys(pdbs))
    return data


def get_curated_3d_mapping():
    """
    Some PDBs are not found automatically so they can be manually recorded
    in a file to be included in processing.
    """
    data = collections.defaultdict(list)
    with open('pdb_full_region_curated.txt', 'r') as f_pdb:
        for line in f_pdb:
            parts = re.split(r'\s+', line)
            rfam_acc = parts[0]
            pdb_id = '{}_{}'.format(parts[1].upper(), parts[2])
            data[rfam_acc].append(pdb_id)
    return data


def merge_3d_mappings(pdb_data, pdb_curated_data):
    """
    Merge PDB data from automatic Rfam mapping and a curated list.
    If a pdb_id is found in the curated list, but not found in the automatic
    mapping, it is added to a merged list.
    """
    for rfam_acc in pdb_curated_data.keys():
        curated_data = pdb_curated_data[rfam_acc]
        if rfam_acc in pdb_data:
            for pdb_id in curated_data:
                if pdb_id not in pdb_data[rfam_acc]:
                    pdb_data[rfam_acc].append(pdb_id)
        else:
            pdb_data[rfam_acc] = curated_data
    return pdb_data


def get_secondary_structure(pdb_id):
    with open('data/fasta/{}.fasta'.format(pdb_id), 'r') as fasta:
        data = fasta.readlines()
    return list(data[2].strip())


def generate_ss_line(structure, line):
    """
    Generate SS Stockholm line.
    """
    new_line = []
    match = re.match(r'^(\w+)(\s+)(\S+)', line)
    for character in match.group(3):
        if character in ['-', '.']:
            new_line.append('.')
        else:
            if structure:
                new_character = structure.pop(0)
                new_line.append(new_character)

    left_column_width = len(match.group(1)) + len(match.group(2))
    left_column = '#=GR {0} {0}_SS'.format(match.group(1)).ljust(left_column_width)
    return left_column + ''.join(new_line), structure


def get_fasta_file(pdb_id):
    """
    Convert JSON basepairs to a
    Some JSON files do not exist because the NMR structures are not annotated
    (e.g. 1ju7).
    """
    pdb_fasta = 'data/fasta/{}.fasta'.format(pdb_id)
    data = fr3d_2d(pdb_id)
    if not data:
        return None
    with open(pdb_fasta, 'w') as f_out:
        f_out.write(data)
    return pdb_fasta


def align_to_seed(rfam_acc, pdb_fasta):
    pdb_sto = 'temp/{}-with-3d.sto'.format(rfam_acc)
    pdb_sto_new = 'temp/{}-with-3d-new.sto'.format(rfam_acc)
    temp_fasta = pdb_fasta.replace('.fasta', '_no_ss.fasta')
    cmd = 'head -2 {} > {}'.format(pdb_fasta, temp_fasta)
    os.system(cmd)
    cmd = 'cmalign --mapali {0} temp/{1}.cm {2} > {3} && cp {3} {0}'.format(pdb_sto, rfam_acc, temp_fasta, pdb_sto_new)
    print(cmd)
    subprocess.check_output(cmd, shell=True)
    return pdb_sto


def add_structure_to_alignment(pdb_id, pdb_sto, structure):
    structure_lines = []
    sequence_lines = []
    with open(pdb_sto, 'r') as f:
        for line in f.readlines():
            if line.startswith(pdb_id):
                new_line, structure = generate_ss_line(structure, line)
                structure_lines.append(new_line.rstrip())
                sequence_lines.append(line.rstrip())
    return structure_lines, sequence_lines


def generate_new_seed(rfam_acc, new_lines, pdb_id):
    data = []
    block_id = 0
    with open('temp/{}-with-3d.sto'.format(rfam_acc), 'r') as f:
        for line in f.readlines():
            if line.startswith('#=GC SS_cons'):
                for lines in new_lines:
                    try:
                        data.append(lines['2d'][block_id])
                    except:
                        raise Exception('temp/{}-with-3d.sto'.format(rfam_acc))
                block_id += 1
            elif line.startswith('#=GR ' + pdb_id):
                continue
            data.append(line.rstrip())
    return data


def generate_new_cm(rfam_acc, pdb_sto):
    new_cm = 'temp/{}.cm'.format(rfam_acc)
    cmd = 'cmbuild -o temp/cmbuild -F {0} {1}'.format(new_cm, pdb_sto)
    print(cmd)
    os.system(cmd)


def validate_pdb_ids(pdb_ids):
    valid = []
    for pdb_id in pdb_ids:
        pdb_code, _ = pdb_id.split('_')
        if pdb_code not in PDB_BLACKLIST:
            valid.append(pdb_id)
    return valid


def process_family(rfam_acc, pdb_ids):
    pdb_ids.sort()
    print('{} PDB structures: '.format(len(pdb_ids)) + ', '.join(pdb_ids))
    for i, pdb_id in enumerate(pdb_ids):
        if i == 0:
            cmd = 'cp data/cm/{0}.cm temp/{0}.cm'.format(rfam_acc)
            os.system(cmd)
            cmd = 'cp data/seed/{0}.seed temp/{0}-with-3d.sto'.format(rfam_acc)
            os.system(cmd)
        pdb_fasta = get_fasta_file(pdb_id)
        if not pdb_fasta:
            print('No pdb fasta file found %s' % pdb_id)
            continue
        pdb_sto = align_to_seed(rfam_acc, pdb_fasta)
        generate_new_cm(rfam_acc, pdb_sto)

    for pdb_id in pdb_ids:
        pdb_fasta = get_fasta_file(pdb_id)
        if not pdb_fasta:
            continue
        new_lines = []
        structure = get_secondary_structure(pdb_id)
        structure_lines, sequence_lines = add_structure_to_alignment(pdb_id, pdb_sto, structure)
        new_lines.append({
            '2d': structure_lines,
            '1d': sequence_lines
        })
        print(pdb_id)
        lines = generate_new_seed(rfam_acc, new_lines, pdb_id)

        with open('temp/{}-with-3d.sto'.format(rfam_acc), 'w') as f:
            for line in lines:
                f.write(line + '\n')

    cmd = 'esl-reformat pfam temp/{0}-with-3d.sto > temp/{0}-final.sto'.format(rfam_acc)
    os.system(cmd)
    rename_accessions(rfam_acc, pdb_ids)
    print('Done: output/{}.sto'.format(rfam_acc))


def rename_accessions(rfam_acc, pdb_ids):
    with open('temp/{}-final.sto'.format(rfam_acc), 'r') as f_in:
        with open('data/output/{}.sto'.format(rfam_acc), 'w') as f_out:
            for line in f_in:
                for pdb_id in pdb_ids:
                    if line.startswith(pdb_id):
                        rnacentral_id = map_pdb_id_to_rnacentral(pdb_id)
                        if not rnacentral_id:
                            continue
                        seq_length = len(re.split(r'\s+', line.strip())[-1].replace('-', '').replace('.', ''))
                        accession = '{rnacentral_id}/1-{length}'.format(rnacentral_id=rnacentral_id, length=seq_length)
                        line = line.replace(pdb_id.ljust(len(accession)), accession)
                        break
                    if line.startswith('#=GR ' + pdb_id):
                        rnacentral_id = map_pdb_id_to_rnacentral(pdb_id)
                        if not rnacentral_id:
                            continue
                        line = line.replace('#=GR ' + pdb_id.ljust(len(rnacentral_id)), '#=GR ' + rnacentral_id)
                        break

                f_out.write(line)


def map_pdb_id_to_rnacentral(pdb_id):
    if not os.path.exists('pdb.tsv'):
        cmd = 'wget http://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/pdb.tsv'
        os.system(cmd)
    rnacentral_id = ''
    with open('pdb.tsv', 'r') as f_pdb:
        for line in f_pdb.readlines():
            if pdb_id in line:
                #URS000080E05C	PDB	3CW1_w	9606	snRNA
                (urs, _, _, taxid, _) = line.strip().split('\t')
                rnacentral_id = urs + '_' + taxid
                break
    if not rnacentral_id:
        # RNAcentral IDs don't exist for mRNA fragments, e.g. 6v4x chain Y
        print('RNAcentral ID not found for {}'.format(pdb_id))
    return rnacentral_id


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('rfam_acc', nargs="+", help='Rfam accession', action='store')
    parser.add_argument("--nocache", help="Recompute output and redownload CM/SEED", action="store_true", default=False)
    args = parser.parse_args()
    rfam_accs = args.rfam_acc
    nocache = args.nocache

    pdb_data = merge_3d_mappings(get_rfam_3d_mapping(), get_curated_3d_mapping())

    if not rfam_accs:
        rfam_accs = pdb_data.keys()

    for rfam_acc in sorted(rfam_accs):
        if rfam_acc in FAMILY_BLACKLIST:
            print('Skipping blacklisted ID {}'.format(rfam_acc))
            continue
        if not re.match(r'RF\d{5}', rfam_acc):
            print('Invalid Rfam accession: {}'.format(rfam_acc))
            continue
        print(rfam_acc)
        if os.path.exists('data/output/{}.sto'.format(rfam_acc)) and not nocache:
            print('Output already exists')
            continue
        if len(pdb_data[rfam_acc]) > SKIP_LARGE_ALIGNMENT:
            print('Skipping alignment with >{} structures'.format(SKIP_LARGE_ALIGNMENT))
            continue
        valid_pdb_ids = validate_pdb_ids(pdb_data[rfam_acc])
        if not valid_pdb_ids:
            print('No valid PDB ids found for {}'.format(rfam_acc))
            continue
        download_rfam_files(rfam_acc, nocache)
        process_family(rfam_acc, valid_pdb_ids)


if __name__ == '__main__':
    main()
