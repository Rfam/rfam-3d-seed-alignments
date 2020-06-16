"""

Add secondary structure based on 3D annotations to Rfam seed alignments.

Infernal binaries should be available in PATH.

Usage:
python add-3d.py

"""


import os
import re

import collections


def download_rfam_files(rfam_acc):
    if not os.path.exists('cm/{}.cm'.format(rfam_acc)):
        cmd = 'wget -O cm/{0}.cm  http://rfam.org/family/{0}/cm'.format(rfam_acc)
        os.system(cmd)
    if not os.path.exists('seed/{}.seed'.format(rfam_acc)):
        cmd = 'wget -O seed/{0}.seed.gz http://rfam.org/family/{0}/alignment/stockholm?gzip=1&download=1 && gunzip seed/{0}.seed.gz'.format(rfam_acc)
        os.system(cmd)


def get_rfam_3d_mapping():
    if not os.path.exists('pdb_full_region.txt'):
        cmd = 'wget ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/pdb_full_region.txt.gz && gunzip pdb_full_region.txt.gz'
        os.system(cmd)
    data = collections.defaultdict(list)
    with open('pdb_full_region.txt', 'r') as f:
        for line in f.readlines():
            if line.startswith('rfam_acc'):
                continue  # skip header
            parts = re.split('\s+', line)
            pdb_id = '{}_{}'.format(parts[1].upper(), parts[2])
            data[parts[0]].append(pdb_id)
    # remove duplicates (e.g. palindromic miRNAs matching in 2 directions)
    for rfam_acc, pdbs in data.iteritems():
        data[rfam_acc] = list(dict.fromkeys(pdbs))
    return data


def get_secondary_structure(pdb_id):
    with open('fasta/{}.fasta'.format(pdb_id), 'r') as f:
        data = f.readlines()
    return list(data[2].strip())


def get_ss_line(structure, line):
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
    left_column = '#=GC {}_SS'.format(match.group(1)).ljust(left_column_width)
    return left_column + ''.join(new_line), structure


def get_fasta_file(pdb_id):
    if not os.path.exists('/Users/apetrov/Desktop/basepairs/{}.json'.format(pdb_id)):
        return None
    pdb_fasta = 'fasta/{}.fasta'.format(pdb_id)
    cmd = 'python json2dotbracket.py {} > {}'.format(pdb_id, pdb_fasta)
    os.system(cmd)
    return pdb_fasta


def align_to_seed(rfam_acc, pdb_fasta):
    pdb_sto = 'temp/{}-with-3d.sto'.format(rfam_acc)
    pdb_sto_new = 'temp/{}-with-3d-new.sto'.format(rfam_acc)
    temp_fasta = pdb_fasta.replace('.fasta', '_no_ss.fasta')
    cmd = 'head -2 {} > {}'.format(pdb_fasta, temp_fasta)
    os.system(cmd)
    cmd = 'cmalign --mapali {0} temp/{1}.cm {2} > {3} && cp {3} {0}'.format(pdb_sto, rfam_acc, temp_fasta, pdb_sto_new)
    print(cmd)
    os.system(cmd)
    return pdb_sto


def add_structure_to_alignment(pdb_id, pdb_sto, structure):
    structure_lines = []
    sequence_lines = []
    with open(pdb_sto, 'r') as f:
        for line in f.readlines():
            if line.startswith(pdb_id):
                new_line, structure = get_ss_line(structure, line)
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
                    data.append(lines['2d'][block_id])
                block_id += 1
            elif line.startswith('#=GR ' + pdb_id):
                continue
            data.append(line.rstrip())
    return data


def generate_new_cm(rfam_acc, pdb_sto):
    new_cm = 'temp/{}.cm'.format(rfam_acc)
    cmd = 'cmbuild -o temp/cmbuild -F {0} {1}'.format(new_cm, pdb_sto)
    os.system(cmd)


def process_family(rfam_acc, pdb_ids):
    for i, pdb_id in enumerate(pdb_ids):
        if i == 0:
            cmd = 'cp cm/{0}.cm temp/{0}.cm'.format(rfam_acc)
            os.system(cmd)
            cmd = 'cp seed/{0}.seed temp/{0}-with-3d.sto'.format(rfam_acc)
            os.system(cmd)
        pdb_fasta = get_fasta_file(pdb_id)
        if not pdb_fasta:
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
        lines = generate_new_seed(rfam_acc, new_lines, pdb_ids[-1])

        with open('temp/{}-with-3d.sto'.format(rfam_acc), 'w') as f:
            for line in lines:
                f.write(line + '\n')

    cmd = 'esl-reformat pfam temp/{0}-with-3d.sto > temp/{0}-final.sto'.format(rfam_acc)
    os.system(cmd)


def main():
    pdb_data = get_rfam_3d_mapping()
    for rfam_acc in pdb_data.keys():
        if rfam_acc not in ['RF00050', 'RF00162']:
            continue
        # if len(pdb_data[rfam_acc]) < 100:
        #     continue
        print(rfam_acc)
        download_rfam_files(rfam_acc)
        process_family(rfam_acc, pdb_data[rfam_acc])


main()
