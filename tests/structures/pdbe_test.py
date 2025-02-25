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

from rfam_3d.structures import pdbe
from rfam_3d.structures.structure_id import PdbChainId, PdbId


@pytest.fixture(scope="module")
def api() -> pdbe.PdbeApi:
    return pdbe.PdbeApi.build(Cache())


@pytest.mark.parametrize(
    "pdb_id,expected",
    [
        (
            PdbId(pdb_id="1S72"),
            pdbe.ExperimentalInfo(
                pdb_id=PdbId(pdb_id="1S72"),
                resolution=2.4,
                method="X-ray diffraction",
            ),
        ),
        (
            PdbId(pdb_id="1cbs"),
            pdbe.ExperimentalInfo(
                pdb_id=PdbId(pdb_id="1cbs"),
                resolution=1.8,
                method="X-ray diffraction",
            ),
        ),
    ],
)
def test_fetches_expected_experimental_info(
    api: pdbe.PdbeApi, pdb_id: PdbId, expected: pdbe.ExperimentalInfo
):
    key = f"pdbe/experiment-{pdb_id.pdb_id}"
    assert api.cache.get(key) is None
    assert api.experiment(pdb_id) == expected
    assert api.cache.get(key) == expected


@pytest.mark.parametrize(
    "chain_id,expected",
    [
        ("1S72_0", 2238),
        ("3cw1_V", 9606),
    ],
)
def test_can_find_expected_taxid(api, chain_id, expected):
    chain_id = PdbChainId.build(chain_id)
    key = f"pdbe/experiment-{chain_id.pdb_id}_{chain_id.chain_id}"
    assert api.cache.get(key) is None
    assert api.molecules(chain_id) == expected
    assert api.cache.get(key) == expected
