
import sys
import json
import os


json_database = '/Users/apetrov/Desktop/basepairs'

def parse_json(filename):
    location = os.path.join(json_database, filename + '.json')
    if not os.path.exists(location):
        print('Does not exist')
        return ('', '')

    with open(location, 'r') as f:
        data = json.load(f)

    sequence = data['sequence']
    structure = list('.' * len(sequence))

    for annotation in data['annotations']:
        if annotation['bp'] == 'cWW':
            structure[annotation['seq_id1']-1] = '('
            structure[annotation['seq_id2']-1] = ')'

    output = """>{}
{}
{}
""".format(filename, sequence, ''.join(structure))
    print(output)
    return (sequence, structure)


parse_json(sys.argv[1])
