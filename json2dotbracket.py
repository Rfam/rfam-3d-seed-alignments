
import sys
import json
import os


json_database_all = '/Users/apetrov/Desktop/basepairs'
json_database_nested = '/Users/apetrov/Desktop/basepairs-no-crossing'

def parse_json(filename):
    if not os.path.exists(filename):
        print('Does not exist')
        return ('', '')

    with open(filename, 'r') as f:
        data = json.load(f)

    sequence = data['sequence']
    structure = list('.' * len(sequence))

    for annotation in data['annotations']:
        if annotation['bp'] == 'cWW':
            structure[annotation['seq_id1']-1] = '('
            structure[annotation['seq_id2']-1] = ')'

    return {
        '1d': sequence,
        '2d': structure,
    }


def generate_output(pdb_id, bp_all, bp_nested):
    all = list(bp_all['2d'])
    nested = list(bp_nested['2d'])
    final = []
    for i, character in enumerate(all):
        if character in ['(', ')'] and nested[i] == '.':
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
    print(output)


pdb_id = sys.argv[1]
filename = os.path.join(json_database_all, pdb_id + '.json')
bp_all = parse_json(filename)
filename = os.path.join(json_database_nested, pdb_id + '.json')
bp_nested = parse_json(filename)

generate_output(pdb_id, bp_all, bp_nested)
