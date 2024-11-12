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

"""This module contains the Disallowed class which is what stores and computes
which structures/families should not be updated.
"""

from __future__ import annotations

import enum

from attr import frozen
from loguru import logger

from rfam_3d.rfam.family import Family
from rfam_3d.update.candidate import Candidate


@enum.unique
class CandidateSkipReason(enum.Enum):
    SKIPPED_STRUCTURE = "Skipped Structure"


@enum.unique
class FamilySkipReason(enum.Enum):
    SKIPPED_FAMILY = "Skipped Family"
    TOO_MANY_MATCHES = "Too many matches"


@frozen
class Disallowed:
    min_resolution: float
    max_structures: int
    skip_pdbs: list[str]
    skip_families: list[str]

    @classmethod
    def empty(cls) -> Disallowed:
        return cls(min_resolution=-1, max_structures=-1, skip_pdbs=[], skip_families=[])

    def skipped_family(self, family: Family) -> None | FamilySkipReason:
        if family.rfam_accession in self.skip_families:
            logger.info(
                "No actions for {}, because the family is blacklisted",
                family.rfam_accession,
            )
            return FamilySkipReason.SKIPPED_FAMILY

        if len(family.matches.matches) > self.max_structures:
            return FamilySkipReason.TOO_MANY_MATCHES

        return None

    def skipped_candidate(self, candidate: Candidate) -> None | CandidateSkipReason:
        if str(candidate.pdb_id) in self.skip_pdbs:
            logger.debug("Skipping {} because it is excluded", candidate.pdb_id)
            return CandidateSkipReason.SKIPPED_STRUCTURE
        return None
