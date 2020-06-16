"""
Copyright [2009-present] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


Export basepair annotations from the RNA 3D Hub database with nucleotide numbering
corresponding to the PDB sequences.

Usage:

1) Create an SSH tunnel to the RNA 3D Hub server:
    ssh -L 3306:127.0.0.1:3306 <username>@rna.bgsu.edu
2) Create and activate a virtual environment:
    virtualenv env
    source env/bin/activate
3) Install all requirements:
    pip install -r requirements.txt
4) Run export:
    python export-basepairs.py

Good test cases:
    pdb_id = '3HL2'
    chain_id = 'E'
"""

import json
import os

from rna3dhub import RNA3DHUB_PASSWORD, RNA3DHUB_DATABASE, RNA3DHUB_USER, RNA3DHUB_LOCAL_HOST


import MySQLdb


def get_basepair_annotations(cursor, pdb_id, chain_id, only_nested):

    query = """
    select A.index as seq_id1, concat(A.number, coalesce(A.ins_code, '')) as 3d_id1, A.nucleotide as nt1, f_lwbp as bp, C.index as seq_id2, C.nucleotide as nt2, concat(C.number, coalesce(C.ins_code, '')) as 3d_id2
    from
    (
        select t3.index + 1 as `index`, t3.`normalized_unit` as `nucleotide`, t2.unit_id, t1.number, t1.ins_code
        from unit_info t1, exp_seq_unit_mapping t2, exp_seq_position t3
        where t1.pdb_id = '{pdb_id}'
        and t1.chain = '{chain_id}'
        and t1.model = 1
        and t1.unit_id = t2.unit_id
        and t2.exp_seq_position_id = t3.exp_seq_position_id
        and (t1.alt_id = 'A' OR t1.alt_id is null)
    ) as A
    JOIN
    (
        select unit_id_1, unit_id_2, f_lwbp, t10.pdb_id
        from unit_pairs_interactions t10, unit_info t11, unit_info t12
        where
        t10.pdb_id = '{pdb_id}'
        and f_lwbp is not null
        and t10.unit_id_1 = t11.unit_id
        and t10.unit_id_2 = t12.unit_id
        and t11.number < t12.number
        and t11.chain = '{chain_id}'
        and t12.chain = '{chain_id}'
        {only_nested}
    ) as B
	JOIN
    (
        select t3.index + 1 as `index`, t3.`normalized_unit` as `nucleotide`, t2.unit_id, t1.number, t1.ins_code
        from unit_info t1, exp_seq_unit_mapping t2, exp_seq_position t3
        where
        t1.pdb_id = '{pdb_id}'
        and t1.chain = '{chain_id}'
        and t1.model = 1
        and t1.unit_id = t2.unit_id
        and t2.exp_seq_position_id = t3.exp_seq_position_id
        and (t1.alt_id = 'A' OR t1.alt_id is null)
    ) as C
    on A.unit_id = B.unit_id_1
    and B.unit_id_2 = C.unit_id
    order by B.pdb_id, A.index
    """

    annotations = []
    cursor.execute(query.format(pdb_id=pdb_id, chain_id=chain_id, only_nested='and f_crossing = 0' if only_nested else '' ))

    rows = cursor.fetchall()

    for row in rows:
        annotations.append(row)

    return annotations


def get_sequence(cursor, pdb_id, chain_id):
    query = """
    select sequence from chain_info where pdb_id = '{pdb_id}' and chain_name = '{chain_id}';
    """
    cursor.execute(query.format(pdb_id=pdb_id, chain_id=chain_id))
    data = cursor.fetchone()
    return data['sequence']


def test_indexing(sequence, annotations):
    for annotation in annotations:
        assert annotation['nt1'] == sequence[annotation['seq_id1']-1], 'Problem with id1'
        assert annotation['nt2'] == sequence[annotation['seq_id2']-1], 'Problem with id2'


def get_structures(cursor):
    query = """
    select pdb_id, chain_name as chain_id
    from chain_info
    where `entity_macromolecule_type` = 'Polyribonucleotide (RNA)'
    and chain_length > 40
    order by pdb_id
    """
    data = []
    cursor.execute(query)
    rows = cursor.fetchall()

    for row in rows:
        data.append(row)

    return data


def run():
    db = MySQLdb.connect(passwd=RNA3DHUB_PASSWORD, db=RNA3DHUB_DATABASE, user=RNA3DHUB_USER, host=RNA3DHUB_LOCAL_HOST)
    cursor = db.cursor(MySQLdb.cursors.DictCursor)

    for structure in get_structures(cursor):
        pdb_id = structure['pdb_id']
        chain_id = structure['chain_id']
        filename_all = os.path.join('basepairs/all', '%s_%s.json' % (pdb_id, chain_id))
        filename_nested = os.path.join('basepairs/nested', '%s_%s.json' % (pdb_id, chain_id))
        if os.path.exists(filename_all) and os.path.exists(filename_nested):
            continue
        print 'Exporting %s %s' % (pdb_id, chain_id)
        nested_bps = get_basepair_annotations(cursor, pdb_id, chain_id, True)
        all_bps = get_basepair_annotations(cursor, pdb_id, chain_id, False)
        data = {
            'pdb_id': pdb_id,
            'chain_id': chain_id,
            'sequence': get_sequence(cursor, pdb_id, chain_id),
        }
        try:
            test_indexing(data['sequence'], nested_bps)
            test_indexing(data['sequence'], all_bps)
        except:
            print('Error: tests failed, skipping')
            continue
        with open(filename_all, 'w') as outfile:
            data['annotations'] = all_bps
            json.dump(data, outfile)
        print 'Created file %s' % filename_all
        with open(filename_nested, 'w') as outfile:
            data['annotations'] = nested_bps
            json.dump(data, outfile)
        print 'Created file %s' % filename_nested


if __name__ == "__main__":
    run()
