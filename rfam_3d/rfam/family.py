# -*- coding: utf-8 -*-

# Copyright [2009-2024] EMBL-European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""This module contains a representation of an Rfam family. This is meant to
"""


from __future__ import annotations

from attr import frozen
from loguru import logger

from rfam_3d.rfam.alignments import Alignment
from rfam_3d.rfam.api import FamilyInfo, RfamApi
from rfam_3d.rfam.matches import FamilyMatch


@frozen
class Family:
    rfam_accession: str
    info: FamilyInfo
    matches: FamilyMatch
    alignment: Alignment


def build_families(
    alignments: list[Alignment],
    matches: list[FamilyMatch],
) -> list[Family]:
    assert alignments, "Must give alignments"
    assert matches, "Must give matches"
    align_map = {}
    for alignment in alignments:
        align_map[alignment.rfam_accession] = alignment
    assert set(align_map.keys()) == set(m.rfam_accession for m in matches)

    api = RfamApi.build()
    families = []
    for family_match in matches:
        accession = family_match.rfam_accession
        families.append(
            Family(
                rfam_accession=accession,
                info=api.info(accession),
                matches=family_match,
                alignment=align_map[accession],
            )
        )
    return families
