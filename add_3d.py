#!/usr/bin/env python3
"""
Add secondary structure based on 3D annotations to Rfam seed alignments.

Infernal binaries should be available in PATH.

Usage:
python add_3d.py all --nocache
python add_3d.py RF00008
"""


import argparse
import os
import re
import shutil
import subprocess

import collections

from colorama import Fore, Style, init

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
    'RF00957', # microRNA matching an rRNA
    'RF02545', # Trypanosomatid SSU
    'RF02546', # Trypanosomatid LSU
]

TEMPDIR = 'temp'


def get_temp_cm_filename(rfam_acc):
    """
    Get the location of a temporary covariance model file.
    """
    return os.path.join(TEMPDIR, f'{rfam_acc}.cm')


def get_rfam_cm_filename(rfam_acc):
    """
    Get the location of an officially released covariance model file.
    """
    return os.path.join('data', 'cm', f'{rfam_acc}.cm')


def get_rfam_seed_filename(rfam_acc):
    """
    Get the location of an official seed file.
    """
    return os.path.join('data', 'seed', f'{rfam_acc}.seed')


def get_temp_3d_seed_filename(rfam_acc):
    """
    Get the location of a temporary seed alignment with 3D structures.
    """
    return os.path.join(TEMPDIR, f'{rfam_acc}-with-3d.sto')


def get_rfam_cm(rfam_acc, nocache):
    """
    Build a covariance model for a given family using the official seed file
    and the --hand option. This is required to be able to transfer the GR and GC
    lines from the official seed file, as the official covariance model cannot
    be used here.
    """
    if not os.path.exists(TEMPDIR):
        os.mkdir(TEMPDIR)
    seed_file = get_rfam_seed_filename(rfam_acc)
    cm_file = get_rfam_cm_filename(rfam_acc)
    if not os.path.exists(cm_file) or nocache:
        print('Building a new CM with the --hand option')
        cmd = f'rm -f {cm_file}.i1*' # remove old cmpress files if found
        subprocess.check_output(cmd, shell=True)
        cmd = f'cmbuild --hand -F {cm_file} {seed_file}'
        subprocess.check_output(cmd, shell=True)


def cmalign_keep_annotations(cm_file, seed_file, fasta, output_file):
    """
    # $1: CM file
    # $2: seed file (used to build CM with --hand)
    # $3: fasta file of seqs to align
    # $4: name of alignment file to create
    """
    # cmd = f'esl-seqstat -a {fasta} | grep ^\\= | awk "{{ print $2 }}" > {output_file}.list'
    cmd = f"grep '>' {fasta} | sed 's/>//' > {output_file}.list"
    subprocess.check_output(cmd, shell=True)

    # align sequences using --mapstr --mapali to keep SS_cons from SEED
    cmd = f'cmalign --sub --notrunc -g --mapstr --mapali {seed_file} -o {output_file}.tmp ' \
          f'{cm_file} {fasta} > {output_file}.cmalign'
    subprocess.check_output(cmd, shell=True)

    # remove orignal seed seqs (to keep SS_cons we need --mapstr --mapali above)
    cmd = f'esl-alimanip --seq-k {output_file}.list {output_file}.tmp > {output_file}.tmp2'
    subprocess.check_output(cmd, shell=True)

    # merge original seed with new seqs to get original annotation from the seed
    cmd = f'esl-alimerge {seed_file} {output_file}.tmp2 > {output_file}'
    subprocess.check_output(cmd, shell=True)


def download_rfam_seed(rfam_acc, nocache):
    """
    Download and uncompress Rfam SEED for a given family.
    """
    if not os.path.exists(TEMPDIR):
        os.mkdir(TEMPDIR)
    filename = get_rfam_seed_filename(rfam_acc)
    if not os.path.exists(filename) or nocache:
        print(f'Downloading the latest {rfam_acc} SEED from SVN')
        url = 'https://xfamsvn.ebi.ac.uk/svn/data_repos/trunk/Families'
        cmd = f'wget -q -a {TEMPDIR}/wget.log -O {filename} ' \
              f'{url}/{rfam_acc}/SEED'
        subprocess.check_output(cmd, shell=True)
        shutil.copyfile(filename, f'{TEMPDIR}/{rfam_acc}-original.seed')


def get_rfam_3d_mapping():
    """
    Parse Rfam PDB mapping file:

    RF00012	6zqb	D4	2	175	54.50	7.2e-12	1	218	c00f0f	0
    RF00015	2n7m	X	1	74	61.90	7.1e-14	1	139	c00f0f	1
    """
    filename = 'pdb_full_region.txt'
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files'
    # url = 'ftp://ftp.ebi.ac.uk/pub/databases/Rfam/.preview'
    if not os.path.exists(filename):
        cmd = f'wget {url}/{filename}.gz && gunzip {filename}.gz'
        subprocess.check_output(cmd, shell=True)
    data = collections.defaultdict(list)
    with open(filename, 'r', encoding='UTF-8') as f_pdb:
        for line in f_pdb:
            if line.startswith('rfam_acc'):
                continue  # skip header
            parts = re.split(r'\s+', line)
            if parts[10] == '0':  # skip redundant clan entries
                continue
            pdb_id = f'{parts[1].upper()}_{parts[2]}'
            data[parts[0]].append(pdb_id)
    # remove duplicates (e.g. palindromic miRNAs matching in 2 directions)
    for rfam_acc, pdbs in data.items():
        data[rfam_acc] = list(dict.fromkeys(pdbs))
    return data


def get_curated_3d_mapping():
    """
    Parse a mamually curated 3D mapping file which supplements the main Rfam
    PDB mapping. This step is needed because some PDBs are not found
    automatically and they need to be recorded manually.
    """
    data = collections.defaultdict(list)
    with open('pdb_full_region_curated.txt', 'r', encoding='UTF-8') as f_pdb:
        for line in f_pdb:
            parts = re.split(r'\s+', line)
            rfam_acc = parts[0]
            pdb_id = f'{parts[1].upper()}_{parts[2]}'
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
    """
    Get secondary structure line from a FASTA file.
    """
    with open(f'data/fasta/{pdb_id}.fasta', 'r', encoding='UTF-8') as fasta:
        data = fasta.readlines()
    return list(data[2].strip())


def generate_ss_line(structure, line):
    """
    Generate SS Stockholm line.
    """
    new_line = []
    match = re.match(r'^(\w+)(\s+)(\S+)', line)
    if not match:
        # URS000080E2F0_93929/1-52
        match = re.match(r'(URS\w{10}_\d+\/\d+-\d+)(\s+)(\S+)', line)

    accession = match.group(1)
    spacer = match.group(2)
    data = match.group(3)

    for character in data:
        if character in ['-', '.']:
            new_line.append('.')
        elif structure:
            new_character = structure.pop(0)
            new_line.append(new_character)

    left_column_width = len(accession) + len(spacer)
    left_column = f'#=GR {accession} {accession}_SS'.ljust(left_column_width)
    return left_column + ''.join(new_line), structure


def get_pdb_fasta_file(pdb_id):
    """
    Convert JSON basepairs to a
    Some JSON files do not exist because the NMR structures are not annotated
    (e.g. 1ju7).
    """
    pdb_fasta = f'data/fasta/{pdb_id}.fasta'
    data = fr3d_2d(pdb_id)
    if not data:
        return None
    with open(pdb_fasta, 'w', encoding='UTF-8') as f_out:
        f_out.write(data)
    return pdb_fasta


def align_to_seed(rfam_acc, pdb_fasta):
    """
    Add PDB sequence to the SEED alignment.
    """
    pdb_sto = get_temp_3d_seed_filename(rfam_acc)
    pdb_sto_new = os.path.join(TEMPDIR, f'{rfam_acc}-with-3d-new.sto')
    temp_fasta = pdb_fasta.replace('.fasta', '_no_ss.fasta')
    cmd = f'head -2 {pdb_fasta} > {temp_fasta}'
    subprocess.check_output(cmd, shell=True)
    print('\tRunning cmalign')
    cmalign_keep_annotations(get_temp_cm_filename(rfam_acc),
                             get_temp_3d_seed_filename(rfam_acc), temp_fasta,
                             pdb_sto_new)
    shutil.copyfile(pdb_sto_new, pdb_sto)


def add_structure_to_alignment(pdb_id, pdb_sto, structure):
    """
    """
    structure_lines = []
    with open(pdb_sto, 'r', encoding='UTF-8') as f_sto:
        for line in f_sto.readlines():
            if line.startswith(pdb_id):
                new_line, structure = generate_ss_line(structure, line)
                structure_lines.append(new_line.rstrip())
    return structure_lines


def add_structure_to_alignment_and_rename(aligned_pdb_id, pdb_id, pdb_sto, structure):
    """
    """
    structure_lines = []
    with open(pdb_sto, 'r', encoding='UTF-8') as f_sto:
        for line in f_sto.readlines():
            if line.startswith(aligned_pdb_id):
                new_line, _ = generate_ss_line(structure, line)
                new_line = new_line.rstrip().replace(f'{aligned_pdb_id}_SS',
                                                     f'{pdb_id}_SS')
                structure_lines.append(new_line)
    return structure_lines


def generate_new_seed(rfam_acc, new_lines, pdb_id):
    """
    Add 2D structure lines to a seed alignment.
    """
    data = []
    block_id = 0
    filename = get_temp_3d_seed_filename(rfam_acc)
    with open(filename, 'r', encoding='UTF-8') as f_seed:
        for line in f_seed.readlines():
            if line.startswith('#=GC SS_cons'):
                data.append(new_lines[block_id])
                block_id += 1
            elif line.startswith('#=GR ' + pdb_id):
                continue
            data.append(line.rstrip())
    with open(filename, 'w', encoding='UTF-8') as f_sto:
        for line in data:
            f_sto.write(line + '\n')


def generate_new_cm(rfam_acc, pdb_sto):
    """
    Build a new covariance model using cmbuild.
    """
    print('\tRunning cmbuild')
    cmd = f'cmbuild --hand -o {TEMPDIR}/cmbuild.log -F {get_temp_cm_filename(rfam_acc)} {pdb_sto}'
    subprocess.check_output(cmd, shell=True)


def validate_pdb_ids(pdb_ids):
    """
    Return a list of PDB ids that excludes blacklisted PDBs.
    """
    valid = set()
    for pdb_id in pdb_ids:
        pdb_code, _ = pdb_id.split('_')
        if pdb_code not in PDB_BLACKLIST:
            valid.add(pdb_id)
        else:
            print(f'Skipping {pdb_id} because it is found in the blacklist')
    return valid


def get_structured_pdb_ids(pdb_ids):
    """
    Return a list of PDB ids with basepairs.
    """
    structured = set()
    for pdb_id in pdb_ids:
        filename = get_pdb_fasta_file(pdb_id)
        with open(filename, 'r', encoding='UTF-8') as f_fasta:
            lines = f_fasta.readlines()
            structure_line = lines[2]
            if '(' in structure_line and ')' in structure_line:
                structured.add(pdb_id)
            else:
                print(f'No basepairs found for {pdb_id} in {filename}')
    return structured


def parse_desc(rfam_acc):
    """
    Fetch and parse DESC file to get basic family metadata using the latest
    committed version of the family.

    AC   RF00008
    ID   Hammerhead_3
    PI   Hammerhead
    DE   Hammerhead ribozyme (type III)
    """
    data = {}
    svn_url = 'https://xfamsvn.ebi.ac.uk/svn/data_repos/trunk/Families'
    desc_file = os.path.join(TEMPDIR, f'{rfam_acc}.desc')
    cmd = f'wget -q -O {desc_file} {svn_url}/{rfam_acc}/DESC'
    subprocess.check_output(cmd, shell=True)
    with open(desc_file, 'r', encoding='latin-1') as f_desc:
        for line in f_desc:
            parts = re.split(r'\s+', line)
            if parts[0] in ['ID', 'DE']:
                data[parts[0]] = ' '.join(parts[1:])
    return data


def add_metadata_gf_lines(rfam_acc, pdb_ids, output_file):
    """
    Include ID and DE lines in the output file to make it easier to explore
    output alignments.

    # STOCKHOLM 1.0
    #=GF ID Hammerhead_3
    #=GF DE Hammerhead ribozyme (type III)
    """
    desc_data = parse_desc(rfam_acc)
    with open(output_file, 'r', encoding='UTF-8') as f_out:
        lines = f_out.readlines()
    unique = []
    with open(output_file, 'w', encoding='UTF-8') as f_out:
        for line_num, line in enumerate(lines):
            if line not in unique:
                f_out.write(line)
                unique.append(line)
            else:
                continue
            if line_num == 0:
                gf_line = f'#=GF ID {desc_data["ID"]}\n'
                f_out.write(gf_line)
                gf_line = f'#=GF DE {desc_data["DE"]}\n'
                f_out.write(gf_line)
                for pdb_id in pdb_ids:
                    gf_line = f'#=GF CC New structure {pdb_id}\n'
                    f_out.write(gf_line)


def get_sequence_length(line):
    """
    Calculate the length of the sequence line excluding gap characters.
    """
    parts = re.split(r'\s+', line.strip())
    sequence = parts[-1]
    sequence = sequence.replace('-', '').replace('.', '')
    return len(sequence)


def rename_accessions(rfam_acc, pdb_ids, rnacentral_ids):
    """
    - Change PDB IDs like 2QUS_A to URS IDs like URS000080DFDA_31504/1-57.
    - Change GR lines with secondary structure from PDB ID to:
        #=GR URS_TAXID/start-stop PDB_ID
    """
    accessions = {}
    left_column_width = 0
    input_file = os.path.join(TEMPDIR, f'{rfam_acc}-final.sto')
    output_file = os.path.join('data', 'output', f'{rfam_acc}.sto')
    with open(input_file, 'r', encoding='UTF-8') as f_in:
        with open(output_file, 'w', encoding='UTF-8') as f_out:
            for line in f_in:
                if line.startswith('#=GS'):
                    continue
                for pdb_id in pdb_ids:
                    if line.startswith(pdb_id):
                        if pdb_id in rnacentral_ids:
                            rnacentral_id = rnacentral_ids[pdb_id]
                        else:
                            continue
                        seq_length = get_sequence_length(line)
                        accession = f'{rnacentral_id}/1-{seq_length}'
                        line = line.replace(pdb_id.ljust(len(accession)), accession)
                        if rnacentral_id not in accessions:
                            accessions[rnacentral_id] = accession
                        elif accessions[rnacentral_id] != accession:
                            print('WARNING: Different start-stop for the same URS')
                        if not left_column_width:
                            left_column_width = len(line) - seq_length
                        break
                    if line.startswith('#=GR ' + pdb_id):
                        if pdb_id in rnacentral_ids:
                            rnacentral_id = rnacentral_ids[pdb_id]
                        else:
                            continue
                        accession = accessions.get(rnacentral_id, f'WARNING {rnacentral_id}')
                        match = re.search(r'^#=GR\s+\S+\s+(\S+_SS)', line)
                        left_column = match.group(0)
                        pdb_id_ss = match.group(1)
                        new_left_column = f'#=GR {accession} {pdb_id_ss}'
                        new_left_column = new_left_column.ljust(len(left_column))
                        line = line.replace(left_column, new_left_column)
                        break
                f_out.write(line)
    return output_file


def map_pdb_id_to_rnacentral(pdb_id):
    """
    Find a URS_taxid corresponding to a specific PDB and chain ID.
    """
    if not os.path.exists('pdb.tsv'):
        cmd = 'wget http://ftp.ebi.ac.uk/pub/databases/RNAcentral/' \
              'current_release/id_mapping/database_mappings/pdb.tsv'
        subprocess.check_output(cmd, shell=True)
    rnacentral_id = ''
    with open('pdb.tsv', 'r', encoding='UTF-8') as f_pdb:
        for line in f_pdb.readlines():
            if pdb_id in line:
                #URS000080E05C	PDB	3CW1_w	9606	snRNA
                (urs, _, _, taxid, _) = line.strip().split('\t')
                rnacentral_id = urs + '_' + taxid
                break
    if not rnacentral_id:
        # RNAcentral IDs don't exist for mRNA fragments, e.g. 6v4x chain Y
        print(f'RNAcentral ID not found for {pdb_id}')
    return rnacentral_id


def get_rfam_family_rnacentral_ids(rfam_acc):
    """
    Get a list of RNAcentral IDs that are already included in the seed alignment
    in order to avoid re-including them.
    """
    seed_file = get_rfam_seed_filename(rfam_acc)
    rnacentral_ids = set()
    with open(seed_file, 'r', encoding='UTF-8') as f_seed:
        for line in f_seed:
            if line.startswith('#'):
                continue
            match = re.search(r'^(URS\w{10}_\d+)\/\d+-\d+', line)
            if match:
                rnacentral_ids.add(match.group(1))
    return rnacentral_ids


def get_rfam_family_pdb_ids(rfam_acc):
    """
    Find PDB ids like 5VTO_R that have already been added to seed alignments
    to prevent re-adding them.
    """
    seed_file = get_rfam_seed_filename(rfam_acc)
    pdb_ids = set()
    with open(seed_file, 'r', encoding='UTF-8') as f_seed:
        for line in f_seed:
            if not line.startswith('#=GR'):
                continue
            match = re.search(r'^#=GR\s+URS\w{10}_\d+\/\d+\-\d+\s+(\w+_\w+)', line)
            if match:
                pdb_ids.add(match.group(1).replace('_SS', ''))
    return pdb_ids


def skip_family(rfam_acc, nocache):
    """
    Return True if a family should be skipped.
    """
    skip = False
    if rfam_acc in FAMILY_BLACKLIST:
        print(f'Skipping blacklisted ID {rfam_acc}')
        return True
    if not re.match(r'RF\d{5}', rfam_acc):
        print(f'Invalid Rfam accession {rfam_acc}')
        return True
    if os.path.exists(f'data/output/{rfam_acc}.sto') and not nocache:
        print('Output already exists')
        return True
    return skip


def get_rnacentral_ids(pdb_ids):
    """
    Get a mapping between PDB IDs and RNAcentral IDs (if available).
    """
    rnacentral_ids = {}
    for pdb_id in pdb_ids:
        rnacentral_id = map_pdb_id_to_rnacentral(pdb_id)
        if rnacentral_id:
            rnacentral_ids[pdb_id] = rnacentral_id
    return rnacentral_ids


def align_pdbs_to_seed(rfam_acc, pdb_ids, rnacentral_ids):
    """
    Align PDB sequences to seed if it has not been done already.
    """
    print(f'{len(pdb_ids)} new PDB structure(s): ' + ', '.join(pdb_ids))

    aligned_rnacentral_ids = []
    aligned_pdbs_ids = {}

    seed_rnacentral_ids = get_rfam_family_rnacentral_ids(rfam_acc)

    aligned_counter = 0
    files_copied = False
    for pdb_id in pdb_ids:
        if pdb_id in rnacentral_ids:
            print(f'{pdb_id} {rnacentral_ids[pdb_id]}')
        else:
            print(pdb_id)
        if pdb_id in rnacentral_ids and rnacentral_ids[pdb_id] in seed_rnacentral_ids:
            msg = f'\tSkipping {pdb_id} because {rnacentral_ids[pdb_id]} ' \
                  f'is already in the original seed alignment'
            print(msg)
            aligned_pdbs_ids[rnacentral_ids[pdb_id]] = pdb_id
            continue
        if pdb_id in rnacentral_ids and rnacentral_ids[pdb_id] not in aligned_rnacentral_ids:
            aligned_rnacentral_ids.append(rnacentral_ids[pdb_id])
            aligned_pdbs_ids[rnacentral_ids[pdb_id]] = pdb_id
        elif pdb_id in rnacentral_ids and rnacentral_ids[pdb_id] in aligned_rnacentral_ids:
            msg = f'\tSkipping {pdb_id} because {rnacentral_ids[pdb_id]} ' \
                  f'is already in the alignment'
            print(msg)
            continue
        else:
            print('\tNo RNAcentral ID but aligning anyway')
        if aligned_counter == 0:
            print('\tStarting with an official seed and cm')
            shutil.copyfile(get_rfam_cm_filename(rfam_acc), get_temp_cm_filename(rfam_acc))
            shutil.copyfile(get_rfam_seed_filename(rfam_acc), get_temp_3d_seed_filename(rfam_acc))
            files_copied = True
        pdb_fasta = get_pdb_fasta_file(pdb_id)
        if not pdb_fasta:
            print(f'No pdb fasta file found {pdb_id}')
            continue
        align_to_seed(rfam_acc, pdb_fasta)
        pdb_sto = get_temp_3d_seed_filename(rfam_acc)
        generate_new_cm(rfam_acc, pdb_sto)
        aligned_counter += 1
        shutil.copyfile(pdb_sto, f'temp/{rfam_acc}-with-3d-{aligned_counter}.sto')

    if not files_copied:
        print('\tStarting with an official seed and cm')
        shutil.copyfile(get_rfam_cm_filename(rfam_acc), get_temp_cm_filename(rfam_acc))
        shutil.copyfile(get_rfam_seed_filename(rfam_acc), get_temp_3d_seed_filename(rfam_acc))
    return aligned_pdbs_ids


def add_secondary_structure(rfam_acc, pdb_ids, rnacentral_ids, aligned_pdbs_ids):
    """
    Add secondary structure GR lines into the alignment.
    """
    pdb_sto = get_temp_3d_seed_filename(rfam_acc)
    for pdb_id in pdb_ids:
        print(f'Adding secondary structure GR lines for {pdb_id}')
        pdb_fasta = get_pdb_fasta_file(pdb_id)
        if not pdb_fasta:
            print('PDB fasta file not found')
            continue
        structure = get_secondary_structure(pdb_id)

        if pdb_id in rnacentral_ids:
            rnacentral_id = rnacentral_ids[pdb_id]
            aligned_pdb_id = aligned_pdbs_ids[rnacentral_id]
        else:
            aligned_pdb_id = pdb_id
        if pdb_id == aligned_pdb_id:
            structure_lines = add_structure_to_alignment(pdb_id, pdb_sto, structure)
        else:
            structure_lines = add_structure_to_alignment_and_rename(aligned_pdb_id, pdb_id, pdb_sto, structure)
        generate_new_seed(rfam_acc, structure_lines, pdb_id)


def finalise_alignment(rfam_acc, pdb_ids, rnacentral_ids):
    """
    Rename accessions, add metadata, and reformat the final alignment.
    """
    pdb_sto = get_temp_3d_seed_filename(rfam_acc)
    cmd = f'esl-reformat pfam {pdb_sto} > temp/{rfam_acc}-final.sto'
    subprocess.check_output(cmd, shell=True)
    output_file = rename_accessions(rfam_acc, pdb_ids, rnacentral_ids)
    add_metadata_gf_lines(rfam_acc, pdb_ids, output_file)
    transfer_gc_annotations(rfam_acc)
    fix_stockholm_whitespace(rfam_acc)
    cmd = f'esl-alistat {output_file}'
    try:
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT,
                                         shell=True, universal_newlines=True)
    except subprocess.CalledProcessError as exc:
        print(f'{Fore.RED}Status : FAIL', exc.returncode, exc.output)
    else:
        print(f'esl-alistat:\n{output}')
    print(f'{Style.BRIGHT}Created {output_file}')


def fix_stockholm_whitespace(rfam_acc):
    """
    Make sure that the Stockholm file is properly formatted and there is no
    misalignment.
    """
    max_width = 0
    new_lines = {}
    filename = os.path.join('data', 'output', f'{rfam_acc}.sto')
    with open(filename, 'r', encoding='UTF-8') as f_sto:
        lines = f_sto.readlines()
        for line in lines:
            # #=GR URS000012D749_4932/1-568 6N7R_R_SS   AUCGCGCG
            match = re.match(r'^(#=GR\s+\S+\s+\S+)\s+(\S+)$', line)
            if match:
                left_column = match.group(1)
                left_column = re.sub(r'\s+', ' ', left_column)
                if len(left_column) > max_width:
                    max_width = len(left_column)
                new_lines[line] = (left_column, match.group(2))
                continue
            # U03476.1/1-572    AUCGCGCG
            match = re.match(r'^(\S+)\s+(\S+)$', line)
            if match:
                left_column = match.group(1)
                if len(left_column) > max_width:
                    max_width = len(left_column)
                new_lines[line] = (left_column, match.group(2))
                continue
            # #=GC RF   AUCGCGCG
            match = re.match(r'^(#=GC\s+\S+)\s+(\S+)$', line)
            if match:
                left_column = match.group(1)
                if len(left_column) > max_width:
                    max_width = len(left_column)
                new_lines[line] = (left_column, match.group(2))
                continue
            new_lines[line] = line
    with open(filename, 'w', encoding='UTF-8') as f_out:
        for line in lines:
            if isinstance(new_lines[line], str):
                f_out.write(new_lines[line])
            else:
                left_column, right_column = new_lines[line]
                new_line = left_column.ljust(max_width + 3) + right_column + '\n'
                f_out.write(new_line)


def transfer_gc_annotations(rfam_acc):
    """
    Transfer all manually curated comment lines from the official Rfam seeds
    into the newly created alignments.

    To find all manually added GC lines:
    egrep -h -o '#=GC\s+\w+' data/seed/*.seed | sort | uniq

    As of 14.7 the following GC lines were present:
    #=GC RNA_ligand_AdoCbl
    #=GC RNA_ligand_AqCbl
    #=GC RNA_ligand_FMN
    #=GC RNA_ligand_Guanidinium
    #=GC RNA_ligand_SAM
    #=GC RNA_ligand_THF_1
    #=GC RNA_ligand_THF_2
    #=GC RNA_ligand_TPP
    #=GC RNA_ligand_ZMP
    #=GC RNA_ligand_fluoride
    #=GC RNA_ligand_guanidine
    #=GC RNA_ligand_preQ1
    #=GC RNA_motif_k_turn
    """
    # convert to pfam format
    pfam_format_seed = f'{get_rfam_seed_filename(rfam_acc)}.pfam'
    cmd = f'esl-reformat pfam {get_rfam_seed_filename(rfam_acc)} > ' \
          f'{pfam_format_seed}'
    subprocess.check_output(cmd, shell=True)

    # find a sequence with the least number of gaps
    data = {}
    reference_sequence = {}
    gc_lines = []
    with open(pfam_format_seed, 'r', encoding='UTF-8') as f_seed:
        for line in f_seed:
            if line.startswith('#=GC') and not line.startswith('#=GC SS_cons') \
               and not line.startswith('#=GC RF'):
                match = re.match(r'^(#=GC\s+\S+)\s+(\S+)$', line)
                label = match.group(1)
                annotation = match.group(2)
                gc_lines.append((label, annotation))
            if line.startswith('#') or len(line) < 10:
                continue
            if not reference_sequence:
                match = re.match(r'^(\S+)\s+(\S+)$', line)
                accession = match.group(1)
                sequence = match.group(2)
                gap_count = sequence.count('.') + sequence.count('-')
                data[accession] = {
                    'sequence': sequence,
                    'gap_count': gap_count,
                }
                # if there are no gaps, use as a reference
                if gap_count == 0:
                    reference_sequence[accession] = data[accession]
    if not gc_lines:
        # nothing to transfer
        return
    if not reference_sequence:
        # loop and select an entry with least gaps
        min_gap = 100000000
        min_accession = ''
        for accession, metadata in data.items():
            if metadata['gap_count'] < min_gap:
                min_accession = accession
        reference_sequence[min_accession] = data[min_accession]

    # map GC symbols to the selected sequence
    new_gc_lines = []
    new_seed = os.path.join('data', 'output', f'{rfam_acc}.sto')
    lines = []
    with open(new_seed, 'r', encoding='UTF-8') as f_new_seed:
        reference_accession = list(reference_sequence.keys())[0]
        ref_seq = reference_sequence[reference_accession]['sequence']
        for line in f_new_seed:
            lines.append(line)
            if not line.startswith(reference_accession):
                continue
            match = re.match(r'^(\S+)\s+(\S+)$', line)
            sequence = match.group(2)
            # generate a new GC line for the same sequence in the new alignment
            for gc_line in gc_lines:
                label, annotation = gc_line
                if len(sequence) == len(annotation):
                    new_gc_lines.append((label, annotation))
                    continue
                new_gc_annotation = []
                # compare with ref_seq
                index = 0
                for symbol in sequence:
                    if symbol == ref_seq[index]:
                        new_gc_annotation.append(annotation[index])
                        index += 1
                    else:
                        new_gc_annotation.append('*')
                new_gc_lines.append((label, ''.join(new_gc_annotation)))

    # rewrite the seed with the new GC line
    with open(new_seed, 'w', encoding='UTF-8') as f_new_seed:
        for line in lines:
            f_new_seed.write(line)
            if line.startswith('#=GC RF'):
                for new_line in new_gc_lines:
                    label, annotation = new_line
                    f_new_seed.write(f'{label}  {annotation}\n')
                    print(f'Transferred a {label} line')


def delete_cached_files():
    """
    Trigger download of fresh Rfam-PDB and RNAcentral-PDB mappings files.
    """
    os.remove('pdb_full_region.txt')
    os.remove('pdb.tsv')


def main():
    """
    Main entrypoint.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('rfam_acc',
                        nargs='+',
                        help='Rfam accession',
                        action='store')
    parser.add_argument('--nocache',
                        help='Recompute output and redownload CM/SEED',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    rfam_accs = args.rfam_acc
    nocache = args.nocache

    if nocache:
        delete_cached_files()
    pdb_data = merge_3d_mappings(get_rfam_3d_mapping(), get_curated_3d_mapping())

    if rfam_accs[0] == 'all':
        rfam_accs = pdb_data.keys()

    init(autoreset=True)

    for rfam_acc in sorted(rfam_accs):
        # if rfam_acc <= 'RF01734':
        print(f'{Fore.MAGENTA}{rfam_acc}')
        #     continue
        if skip_family(rfam_acc, nocache):
            continue
        if len(pdb_data[rfam_acc]) > SKIP_LARGE_ALIGNMENT:
            print(f'Skipping alignment with >{SKIP_LARGE_ALIGNMENT} structures')
            continue
        valid_pdb_ids = validate_pdb_ids(pdb_data[rfam_acc])
        if not valid_pdb_ids:
            print('No valid PDB ids found')
            continue
        structured_valid_pdb_ids = get_structured_pdb_ids(valid_pdb_ids)
        if not structured_valid_pdb_ids:
            print('No PDB structures with basepairs found')
            continue
        download_rfam_seed(rfam_acc, nocache)
        get_rfam_cm(rfam_acc, nocache)
        rfam_pdb_ids = get_rfam_family_pdb_ids(rfam_acc)
        new_pdb_ids = structured_valid_pdb_ids - rfam_pdb_ids
        if new_pdb_ids:
            print(f'{len(rfam_pdb_ids)} PDB already in SEED')
        else:
            print('No new PDB ids found')
            continue
        pdb_ids = list(new_pdb_ids)
        pdb_ids.sort()
        rnacentral_ids = get_rnacentral_ids(pdb_ids)
        aligned_pdbs_ids = align_pdbs_to_seed(rfam_acc, pdb_ids, rnacentral_ids)
        add_secondary_structure(rfam_acc, pdb_ids, rnacentral_ids, aligned_pdbs_ids)
        finalise_alignment(rfam_acc, pdb_ids, rnacentral_ids)
        print(f'Added {len(pdb_ids)} new PDB ID(s) from {len(rnacentral_ids)} RNAcentral ID(s)')


if __name__ == '__main__':
    main()
