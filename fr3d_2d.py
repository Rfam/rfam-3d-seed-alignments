#!/usr/bin/env python3
"""
Generate a FASTA file with secondary structure in dot bracket notation
based on FR3D annotations of an RNA 3D structure.

Usage:
python fr3d_2d.py 2QUS_B

>2QUS_B
GGGAGCCCUGUCACCGGAUGUGCUUUCCGGUCUGAUGAGUCCGUGAGGACAAAACAGGGCUCCCGAAUU
.((((((((((.((((((.....{.))))))(....).((((...}))))...))))))))))......
"""

import json
import os
import sys

BASEPAIRS_ALL = 'data/basepairs/all'
BASEPAIRS_NESTED = 'data/basepairs/nested'


def get_rna3dhub_annotations(pdb_id):
    """
    Get FR3D annotations in JSON format.
    """
    bps_all = os.path.join(BASEPAIRS_ALL, '{}.json'.format(pdb_id))
    bps_nested = os.path.join(BASEPAIRS_NESTED, '{}.json'.format(pdb_id))
    pdb_code, chain = pdb_id.split('_')
    if not os.path.exists(bps_all):
        cmd = "curl -s -o {} 'http://rna.bgsu.edu/rna3dhub/rest/getSequenceBasePairs?pdb_id={}&chain={}&only_nested=False'".format(bps_all, pdb_code, chain)
        os.system(cmd)
    if not os.path.exists(bps_nested):
        cmd = "curl -s -o {} 'http://rna.bgsu.edu/rna3dhub/rest/getSequenceBasePairs?pdb_id={}&chain={}&only_nested=True'".format(bps_nested, pdb_code, chain)
        os.system(cmd)


def parse_json(filename):
    """
    Some JSON files do not exist because the NMR structures are not annotated
    (e.g. 1ju7).
    """
    if not os.path.exists(filename):
        print('File {} does not exist'.format(filename))
        return None

    with open(filename, 'r') as f:
        data = json.load(f)

    sequence = data['sequence']
    structure = list('.' * len(sequence))

    for annotation in data['annotations']:
        if annotation['bp'] == 'cWW':
            structure[int(annotation['seq_id1'])-1] = '('
            structure[int(annotation['seq_id2'])-1] = ')'

    return {
        '1d': sequence,
        '2d': structure,
    }


def generate_output(pdb_id, bp_all, bp_nested):
    bp_all_list = list(bp_all['2d'])
    bp_nested_list = list(bp_nested['2d'])
    final = []
    for i, character in enumerate(bp_all_list):
        if character in ['(', ')'] and bp_nested_list[i] == '.':
            if character == '(':
                final.append('{')
            elif character == ')':
                final.append('}')
        else:
            final.append(character)
    output = """>{}
{}
{}
""".format(pdb_id, bp_all['1d'], ''.join(final))
    return output


def fr3d_2d(pdb_id):
    get_rna3dhub_annotations(pdb_id)
    bp_all = parse_json(os.path.join(BASEPAIRS_ALL, pdb_id + '.json'))
    bp_nested = parse_json(os.path.join(BASEPAIRS_NESTED, pdb_id + '.json'))
    if not bp_all or not bp_nested:
        return None
    return generate_output(pdb_id, bp_all, bp_nested)


if __name__ == '__main__':
    print(fr3d_2d(sys.argv[1]))
