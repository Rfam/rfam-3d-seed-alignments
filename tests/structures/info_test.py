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

from rfam_3d.structures import info
from rfam_3d.structures.structure_id import PdbChainId


@pytest.fixture(scope="module")
def api() -> info.StructureLookup:
    return info.StructureLookup.build(Cache())


@pytest.mark.parametrize("chain_id,expected", [])
def test_can_fetch_structure_info(
    api: info.StructureLookup, chain_id: PdbChainId, expected: info.ChainInfo
):
    assert api.fetch_structure_info(chain_id) == expected
