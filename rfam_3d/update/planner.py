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

"""This module is what determines what steps need to be taken to update an
alignment.
"""

from __future__ import annotations

from attr import frozen
from loguru import logger

from rfam_3d.disallow import Disallowed
from rfam_3d.rfam.family import Family
from rfam_3d.structures.info import StructureLookup
from rfam_3d.update.actions import (
    Action,
    CompleteCandidateAction,
    CompleteFamily,
    FamilyUpdate,
    SkippedFamily,
)
from rfam_3d.update.candidate import Candidate


class DuplicateSequence(Exception):
    """This is raised if the FamilyActions is asked to align the same sequence twice."""


@frozen
class Planner:
    disallowed: Disallowed
    structure_fetcher: StructureLookup

    @classmethod
    def build(
        cls, structure_fetcher: StructureLookup, disallowed: Disallowed
    ) -> Planner:
        return cls(structure_fetcher=structure_fetcher, disallowed=disallowed)

    def actions(self, family: Family) -> Action:
        """This determines which actions to take given what is already mapped to a
        family and what the family matches. If None is returned then there are no
        actions to take.

        This first filters out all things disallowed, this may lead to nothing
        being done.

        Next this considers all families that match the alignment and for each
        chain checks:

        1) If the structure is already added. This checks if there is any sequence
           with GR annotations for the given structure. This does not require that
           is the PDB sequence because that may have been edited since matches can
           be incomplete. If this is the case then there is nothing to do for that
           structure. If it is not then the structure is add structures_to_add.

        2) The sequence is already aligned. If the PDB is not present then it may
           be part of a sequence which is already aligned. If that is the case then
           there is no need to add the sequence again. This determines if the
           sequence is already aligned by checking if there is a sequence with the
           expected URS if known, or hash otherwise.

        3) If the sequence has already been added to sequences_to_add. If so then
           this need only add the structure

        If none are true then the chain is added to structures_to_add and the
        sequence is added to sequences_to_add.
        """

        update = FamilyUpdate.empty(family)
        for match in family.matches:
            pdb_info = self.structure_fetcher.fetch_structure_info(match.chain_id)
            if not pdb_info:
                logger.warning(
                    "Failed to get information about {}, skipping", match.chain_id
                )
                update.no_structure_info(match)
                continue

            candidate = Candidate(match=match, chain_info=pdb_info)
            if accession := family.alignment.sequence_annotated_with(
                candidate.chain_id
            ):
                logger.info(
                    "Structure {} is already annotated on {}",
                    str(candidate.chain_id),
                    str(accession),
                )
                update.existing_match(accession, candidate)
                continue

            if reason := self.disallowed.skipped_candidate(candidate):
                logger.info(
                    "Skipping candidate {}, reason {}", candidate.chain_id, reason
                )
                update.skip_candidate(candidate, reason)
                continue

            if accession := family.alignment.known_sequence(candidate.chain_info.seq()):
                logger.info(
                    "PDB {} has sequence {}, which is already aligned",
                    str(candidate.chain_id),
                    str(accession),
                )
                update.annotate_structure(accession, candidate)
                continue

            if accession := update.planned_accession(candidate):
                logger.info(
                    "PDB {} has sequence {}, which will be aligned",
                    str(pdb_info.chain_id),
                    str(accession),
                )
                update.annotate_structure(accession, candidate)
                continue

            accession = pdb_info.unique_accession()
            logger.info(
                "PDB {} has sequence {}, which is not aligned",
                str(pdb_info.chain_id),
                str(accession),
            )
            update.align_and_annotate(candidate)

        if reason := self.disallowed.skipped_family(family):
            logger.info(
                "Family {} is skipped, cause: {}", family.rfam_accession, reason
            )
            return SkippedFamily.from_update(update, reason)

        if update.is_complete():
            logger.info("No actions to take for {}", family.rfam_accession)
            actions = []
            for action in update.actions:
                assert isinstance(action, CompleteCandidateAction)
                actions.append(action)
            return CompleteFamily(family=family, actions=actions)
        assert update.actions
        return update
