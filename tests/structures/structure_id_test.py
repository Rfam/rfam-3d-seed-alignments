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

from rfam_3d.structures.structure_id import PdbChainId, PdbId


def test_pdb_id_correctly_cases_pdb_ids():
    assert PdbId(pdb_id="1J5E").pdb_id == "1j5e"
    assert PdbId(pdb_id="1j5e").pdb_id == "1j5e"


def test_pdb_chain_id_builds_correct_ss_id():
    chain_id = PdbChainId(pdb_id=PdbId(pdb_id="1J5e"), chain_id="A")
    assert chain_id.secondary_structure_id() == "1J5E_A_SS"
