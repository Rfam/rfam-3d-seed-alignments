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

from rfam_3d.rfam import api as rfam


@pytest.fixture(scope="module")
def api() -> rfam.RfamApi:
    return rfam.RfamApi.build()


@pytest.mark.parametrize(
    "accession,expected",
    [
        (
            "RF00177",
            rfam.FamilyInfo(
                id="SSU_rRNA_bacteria",
                accession="RF00177",
                rna_type="Gene; rRNA;",
                num_seed=99,
                num_full=42628,
                description="Bacterial small subunit ribosomal RNA",
            ),
        )
    ],
)
def test_can_fetch_expected_info(
    api: rfam.RfamApi, accession: str, expected: rfam.FamilyInfo
):
    assert api.info(accession) == expected


@pytest.mark.parametrize(
    "accession",
    [
        "Missing",
        "RF00000",
    ],
)
def test_fails_given_bad_accession(api, accession):
    with pytest.raises(Exception):
        api.info(accession)
