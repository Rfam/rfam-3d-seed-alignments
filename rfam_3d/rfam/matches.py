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

from __future__ import annotations

"""This modules deals with loading information from the mapping file.
"""

import csv
import typing as ty

import cattrs
from attr import frozen
from loguru import logger

from rfam_3d.structures.structure_id import PdbChainId, PdbId


@frozen
class Match:
    """This represents a single match of a PDB sequence to a Rfam family."""

    rfam_accession: str
    chain_id: PdbChainId
    sequence_start_position: int
    sequence_stop_position: int
    bit_score: float
    e_value: float
    cm_start_position: int
    cm_end_position: int
    hex_color: str
    is_significant: bool


@frozen
class FamilyMatch:
    rfam_accession: str
    matches: list[Match]

    @classmethod
    def empty(cls, accession: str) -> FamilyMatch:
        return cls(rfam_accession=accession, matches=[])

    def add_match(self, instance: Match):
        """Add a match to this family. It must be for the same Rfam
        accession."""
        assert instance.rfam_accession == self.rfam_accession
        self.matches.append(instance)

    def __iter__(self) -> ty.Iterator[Match]:
        return iter(self.matches)

    def __bool__(self) -> bool:
        """Determine if this FamilyMatch is empty."""
        return len(self.matches) == 0


def load_all(handle: ty.TextIO) -> ty.List[FamilyMatch]:
    """Loads all matches in the given file handle.
    The file must be a TSV with no header with the following fields:

    - rfam_accession
    - pdb_id
    - chain_id
    - sequence_start_position
    - sequence_stop_position
    - bit_score
    - e_value
    - cm_start_position
    - cm_end_position
    - hex_color
    - is_significant
    """

    logger.trace("Loading matches from {}", handle.name)
    reader = csv.DictReader(
        handle,
        delimiter="\t",
        fieldnames=[
            "rfam_accession",
            "pdb_id",
            "chain_id",
            "sequence_start_position",
            "sequence_stop_position",
            "bit_score",
            "e_value",
            "cm_start_position",
            "cm_end_position",
            "hex_color",
            "is_significant",
        ],
    )
    matches = {}
    for row in reader:
        pdb_id = row.pop("pdb_id")
        assert isinstance(pdb_id, str)
        row["chain_id"] = {"pdb_id": {"pdb_id": pdb_id}, "chain_id": row["chain_id"]}
        instance = cattrs.structure(row, Match)
        family = matches.get(
            instance.rfam_accession, FamilyMatch.empty(instance.rfam_accession)
        )
        family.add_match(instance)
        matches[family.rfam_accession] = family
    return list(matches.values())
