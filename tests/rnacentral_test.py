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
from Bio.Seq import Seq
from diskcache import Cache

from rfam_3d import rnacentral as rnac


@pytest.fixture(scope="module")
def api():
    return rnac.RnacentralApi.build(Cache())


@pytest.mark.parametrize(
    "sequence,taxid,expected",
    [
        (
            Seq(
                "AUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGUAGAGAGAAGCUUGCUUCUCUUGAGAGCGGCGGACGGGUGAGUAAUGCCUAGGAAUCUGCCUGGUAGUGGGGGAUAACGCUCGGAAACGGACGCUAAUACCGCAUACGUCCUACGGGAGAAAGCAGGGGACCUUCGGGCCUUGCGCUAUCAGAUGAGC"
            ),
            77133,
            "URS0000000001_77133",
        ),
        (
            Seq(
                "AUUGAACGCUGGCGGCAGGCCUAACACAUGCAAGUCGAGCGGUAGAGAGAAGCUUGCUUCUCUUGAGAGCGGCGGACGGGUGAGUAAUGCCUAGGAAUCUGCCUGGUAGUGGGGGAUAACGCUCGGAAACGGACGCUAAUACCGCAUACGUCCUACGGGAGAAAGCAGGGGACCUUCGGGCCUUGCGCUAUCAGAUGAGC"
            ),
            2,
            "URS0000000001_2",
        ),
        (
            Seq(
                "ATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAGCGGTAGAGAGAAGCTTGCTTCTCTTGAGAGCGGCGGACGGGTGAGTAATGCCTAGGAATCTGCCTGGTAGTGGGGGATAACGCTCGGAAACGGACGCTAATACCGCATACGTCCTACGGGAGAAAGCAGGGGACCTTCGGGCCTTGCGCTATCAGATGAGC"
            ),
            77133,
            "URS0000000001_77133",
        ),
        (
            Seq("N"),
            77133,
            None,
        ),
    ],
)
def test_can_find_rnacentral_id(
    api: rnac.RnacentralApi, sequence: Seq, taxid: int, expected: str | None
):
    assert api.rnacentral_id(sequence, taxid) == expected
