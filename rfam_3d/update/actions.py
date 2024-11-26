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

import json
import typing as ty

import cattrs
from attr import frozen
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger

from rfam_3d.disallow import CandidateSkipReason, FamilySkipReason
from rfam_3d.rfam.accessions import SequenceAccession
from rfam_3d.rfam.family import Family
from rfam_3d.rfam.matches import Match
from rfam_3d.structures.structure_id import PdbChainId, PdbId
from rfam_3d.update.candidate import Candidate
from rfam_3d.update.structure_actions import StructureActions
from rfam_3d.utils import assert_never


@frozen
class StructurelessCandidate:
    """This indicates a candidate where no structure could be found for the match.
    This is an error case."""

    match: Match

    @property
    def pdb_id(self) -> PdbId:
        return self.match.chain_id.pdb_id

    @property
    def chain_id(self) -> PdbChainId:
        return self.match.chain_id


@frozen
class SkippedCandidate:
    """This indicates that a candidate is being skipped and why."""

    candidate: Candidate
    reason: CandidateSkipReason

    @property
    def pdb_id(self) -> PdbId:
        return self.candidate.pdb_id

    @property
    def chain_id(self) -> PdbChainId:
        return self.candidate.chain_id

    @property
    def match(self) -> Match:
        return self.candidate.match


@frozen
class AlreadyPresentCandidate:
    """This indicates that a candidate has already been added."""

    candidate: Candidate
    accession: SequenceAccession

    @property
    def pdb_id(self) -> PdbId:
        return self.candidate.pdb_id

    @property
    def chain_id(self) -> PdbChainId:
        return self.candidate.chain_id

    @property
    def match(self) -> Match:
        return self.candidate.match


@frozen
class AlignCandidateSequence:
    """This represents the need to align and annotate a candidate in an
    alignment."""

    candidate: Candidate
    accession: SequenceAccession

    @property
    def pdb_id(self) -> PdbId:
        return self.candidate.pdb_id

    @property
    def chain_id(self) -> PdbChainId:
        return self.candidate.chain_id

    @property
    def match(self) -> Match:
        return self.candidate.match


@frozen
class AddCandidateStructure:
    """This represents the need to annotate a candidate on a specific
    sequence in the alignment."""

    candidate: Candidate
    accession: SequenceAccession

    @property
    def pdb_id(self) -> PdbId:
        return self.candidate.pdb_id

    @property
    def chain_id(self) -> PdbChainId:
        return self.candidate.chain_id

    @property
    def match(self) -> Match:
        return self.candidate.match


CandidateAction = (
    StructurelessCandidate
    | SkippedCandidate
    | AlreadyPresentCandidate
    | AlignCandidateSequence
    | AddCandidateStructure
)

CompleteCandidateAction = (
    StructurelessCandidate | SkippedCandidate | AlreadyPresentCandidate
)


@frozen
class NoAcceptedMatches:
    """This is the case where all matches within a family are skipped"""

    family: Family
    actions: list[StructurelessCandidate | SkippedCandidate]


@frozen
class CompleteFamily:
    """This represents a family where all candidates have already been added."""

    family: Family
    actions: list[CompleteCandidateAction]


@frozen
class SkippedFamily:
    """This is a family which has been skipped for some reason."""

    family: Family
    reason: FamilySkipReason
    actions: list[CandidateAction]

    @classmethod
    def from_update(
        cls, update: FamilyUpdate, reason: FamilySkipReason
    ) -> SkippedFamily:
        return cls(
            family=update.family,
            reason=reason,
            actions=update.actions,
        )


@frozen
class FamilyUpdate:
    """This represents a family which needs some sort of updates."""

    family: Family
    actions: list[CandidateAction]

    @classmethod
    def empty(cls, family: Family) -> FamilyUpdate:
        """Create an empty FamilyUpdate, which has no actions.

        :family Family: The Family to create updates for.
        """
        return cls(family=family, actions=[])

    def skip_candidate(self, candidate: Candidate, reason: CandidateSkipReason):
        """Mark a match as skipped for the given reason."""
        self.actions.append(SkippedCandidate(candidate=candidate, reason=reason))

    def existing_match(self, accession: SequenceAccession, candidate: Candidate):
        """Mark a match as already annotated in the alignment."""
        # assert self.family.has_accession(accession)
        self.actions.append(
            AlreadyPresentCandidate(candidate=candidate, accession=accession)
        )

    def no_structure_info(self, match: Match):
        """Mark a match as having failed to get structure information."""
        self.actions.append(StructurelessCandidate(match=match))

    def planned_accession(self, candidate: Candidate) -> None | SequenceAccession:
        """Check if the given structure is already planned to be aligned and if
        it is return the accession that will be used.

        :pdb_info: The chain to align.
        """
        for action in self.actions:
            if isinstance(
                action, AlignCandidateSequence
            ) and action.candidate.same_sequence(candidate):
                return action.accession
        return None

    def annotate_structure(self, accession: SequenceAccession, candidate: Candidate):
        """Add a structure which will be annotated as an annotation on the given
        accession.
        """
        assert self.family.alignment.known_accession(accession) or any(
            c.accession == accession
            for c in self.actions
            if isinstance(c, AlignCandidateSequence)
        ), "Cannot annotate onto unknown structure"
        self.actions.append(
            AddCandidateStructure(candidate=candidate, accession=accession)
        )

    def align_candidate_sequence(
        self, accession: SequenceAccession, candidate: Candidate
    ):
        self.actions.append(
            AlignCandidateSequence(
                accession=accession,
                candidate=candidate,
            )
        )

    def align_and_annotate(self, candidate: Candidate):
        accession = candidate.chain_info.unique_accession()
        self.align_candidate_sequence(accession, candidate)
        self.annotate_structure(accession, candidate)

    def is_complete(self) -> bool:
        for action in self.actions:
            if not isinstance(action, CompleteCandidateAction):
                return False
        return True

    def write_data(self, seq_out: ty.TextIO, struct_out: ty.TextIO):
        seen_ids = set()
        sequences: list[SeqRecord] = []
        structures = StructureActions.empty()
        for action in self.actions:
            match action:
                case AlreadyPresentCandidate():
                    logger.debug("Match {} already exists", action.match.chain_id)
                case StructurelessCandidate():
                    logger.debug("Skipping match without structure {}", action.match)
                case SkippedCandidate():
                    logger.debug(
                        "Not writing sequence for skipped candidate {}", action
                    )
                case AddCandidateStructure():
                    logger.debug(
                        "Not writing sequence for candidate {} with known sequence",
                        action,
                    )
                    structures.annotate_structure(
                        action.accession.raw, action.candidate.chain_info
                    )
                case AlignCandidateSequence():
                    logger.debug("Will add sequence and structure from: {}", action)
                    accession = action.accession.raw
                    if accession in seen_ids:
                        logger.trace("Skipping already written id {}", accession)
                        continue
                    candidate = action.candidate
                    sequences.append(SeqRecord(candidate.full_sequence(), id=accession))
                    structures.annotate_structure(accession, candidate.chain_info)
                case _:
                    assert_never(action)

        logger.info(
            "{} sequences to align for {}", len(sequences), self.family.rfam_accession
        )
        SeqIO.write(sequences, seq_out, "fasta")

        logger.info(
            "{} structures to annotate for {}",
            len(structures),
            self.family.rfam_accession,
        )
        raw = cattrs.unstructure(structures)
        json.dump(raw, struct_out)

    def has_aligned(self) -> bool:
        """Check if there is already a structure aligned to this family."""
        return self.family.alignment.has_structure()

    def has_sequences(self) -> bool:
        """Check if this update requires adding any sequences."""
        return any(isinstance(a, AlignCandidateSequence) for a in self.actions)


Action = SkippedFamily | CompleteFamily | NoAcceptedMatches | FamilyUpdate
