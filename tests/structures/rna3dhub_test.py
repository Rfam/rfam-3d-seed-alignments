# -*- coding: utf-8 -*-

# Copyright [2009-2025] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest
from diskcache import Cache

from rfam_3d.structures import rna3dhub as r3d
from rfam_3d.structures.structure_id import PdbChainId, PdbId


@pytest.fixture(scope="module")
def api() -> r3d.Rna3dHubApi:
    return r3d.Rna3dHubApi.build(Cache())


@pytest.mark.parametrize("nested,complete,expected", [])
def test_basepairing_can_create_expected_pk_strings(nested, complete, expected):
    bp = r3d.Basepairing(
        chain_id=PdbChainId(pdb_id=PdbId(pdb_id="1S72"), chain_id="9"),
        sequence="A" * len(nested),
        nested=nested,
        complete=complete,
    )
    assert bp.psuedoknotted() == expected


@pytest.mark.parametrize(
    "chain_id,expected",
    [
        (
            PdbChainId(pdb_id=PdbId(pdb_id="1S72"), chain_id="9"),
            r3d.Basepairing(
                chain_id=PdbChainId(pdb_id=PdbId(pdb_id="1S72"), chain_id="9"),
                sequence="UUAGGCGGCCACAGCGGUGGGGUUGCCUCCCGUACCCAUCCCGAACACGGAAGAUAAGCCCACCAGCGUUCCGGGGAGUACUGGAGUGCGCGAGCCUCUGGGAAACCCGGUUCGCCGCCACC",
                nested="...((((((....((((((((.....(((((((...(.....)...))))..)))...)))))).)).(((((((.....((((((.((....))))))))....))))))).))))))...",
                complete="...((((((....((((((((.....(((((((...(.....)...))))..)))...)))))).)).(((((((.....((((((.((....))))))))....))))))).))))))...",
            ),
        )
    ],
)
def test_fetches_all_basepairs(api, chain_id, expected):
    assert api.basepairing(chain_id) == expected
