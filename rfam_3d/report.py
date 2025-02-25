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

import enum
import json
import typing as ty

from attr import frozen
from attrs import Attribute, fields, has
from cattrs import Converter
from cattrs.gen import make_dict_structure_fn, make_dict_unstructure_fn, override
from loguru import logger

from rfam_3d.rfam.accessions import SequenceAccession
from rfam_3d.rfam.matches import Match
from rfam_3d.structures.info import ChainInfo
from rfam_3d.structures.structure_id import PdbChainId
from rfam_3d.update.actions import (
    Action,
    AddCandidateStructure,
    AlignCandidateSequence,
    AlreadyPresentCandidate,
    CandidateAction,
    CompleteFamily,
    FamilyUpdate,
    NoAcceptedMatches,
    SkippedCandidate,
    SkippedFamily,
    StructurelessCandidate,
)
from rfam_3d.utils import assert_never, to_camel_case


@enum.unique
class SequenceStatus(enum.Enum):
    NEW_SEQUENCE = "new_sequence"
    KNOWN_SEQUENCE = "known_sequence"


@enum.unique
class BasepairingStatus(enum.Enum):
    NEW_STRUCTURE = "new_pairing"
    KNOWN_STRUCTURE = "known_pairing"


@enum.unique
class MatchStatus(enum.Enum):
    """This indicates why a chain matches a family. Most families are
    automatically assigned, but it is possible for them to be manually
    assigned.

    AUTOMATED_MATCH - The chain matches the family above the threshold.
    MANUALLY_FORCED - A curator assigned this match manually.
    """

    AUTOMATED_MATCH = "automated_match"
    MANUALLY_FORCED = "manually_forced"


@enum.unique
class PkStatus(enum.Enum):
    """The PkStatus indiciates if the chain has at least one pseudoknot which is new.

    EXISTING_PK - All pseudoknots already exist in the family
    NEW_PK - If there is at least one pk which is new.
    NO_PK - The chain has now pseudoknotss.
    """

    EXISTING_PK = "existing_pk"
    NEW_PK = "new_pk"
    NO_PK = "no_pk"


@enum.unique
class ChainStatus(enum.Enum):
    """A ChainStatus represents if the chain is new for the family.

    ALREADY_PRESENT - The chain is already a member of the family.
    ERROR - There was an error fetching information about the chain.
    NEW_CHAIN - The chain is not already present in the alignment.
    SKIPPED - The chain was skipped by the disallow list.
    """

    ALREADY_PRESENT = "already_present"
    ERROR = "error"
    NEW_CHAIN = "new_chain"
    SKIPPED = "skipped"

    @classmethod
    def status_for(cls, action: CandidateAction) -> None | ChainStatus:
        match action:
            case StructurelessCandidate():
                return cls.ERROR
            case SkippedCandidate():
                return cls.SKIPPED
            case AlreadyPresentCandidate():
                return cls.ALREADY_PRESENT
            case AlignCandidateSequence():
                return None
            case AddCandidateStructure():
                return None
            case _:
                assert_never(action)


@enum.unique
class FamilyStatus(enum.Enum):
    """A FamilyStatus models if a family is ready for curation or if it is in
    some other state.

    COMPLETE - All structures have been added to the family
    CURATE - This family has no structures added, but it could
    INCOMPLETE - The family has at least one structure already added
    NO_VALID - All matching structures were rejected from the family
    SKIPPED - The family was skipped by the disallow list
    """

    COMPLETE = "complete"
    CURATE = "curate"
    INCOMPLETE = "incomplete"
    NO_VALID = "no_matches"
    SKIPPED = "skipped"

    @classmethod
    def status_for(cls, action: Action) -> FamilyStatus:
        """Generate a FamilyStatus for a given Action."""
        match action:
            case SkippedFamily():
                return cls.SKIPPED
            case NoAcceptedMatches():
                return cls.NO_VALID
            case CompleteFamily():
                return cls.COMPLETE
            case FamilyUpdate():
                if action.has_aligned():
                    return cls.INCOMPLETE
                return cls.CURATE
            case _:
                assert_never(action)


@frozen
class ChainReport:
    chain_id: str
    status: ChainStatus
    basepairing_status: None | BasepairingStatus
    dot_bracket: str
    has_base_pairs: bool
    has_pseudoknots: bool
    new_basepair_count: None | int
    removed_basepair_count: None | int

    @classmethod
    def from_info(
        cls,
        info: ChainInfo,
        status: ChainStatus,
        basepairing_status: None | BasepairingStatus,
    ) -> ChainReport:
        return cls(
            chain_id=str(info.chain_id),
            status=status,
            basepairing_status=basepairing_status,
            dot_bracket=info.basepairing.psuedoknotted(),
            has_base_pairs=info.basepairing.has_basepairs,
            has_pseudoknots=info.basepairing.has_pk,
            new_basepair_count=None,
            removed_basepair_count=None,
        )


@frozen
class MatchReport:
    sequence_start_position: int
    sequence_stop_position: int
    bit_score: float
    e_value: float
    cm_start_position: int
    cm_end_position: int
    hex_color: str
    is_significant: bool
    match_status: MatchStatus

    @classmethod
    def from_match(cls, match: Match) -> MatchReport:
        return cls(
            sequence_start_position=match.sequence_start_position,
            sequence_stop_position=match.sequence_stop_position,
            bit_score=match.bit_score,
            e_value=match.e_value,
            cm_start_position=match.cm_start_position,
            cm_end_position=match.cm_end_position,
            hex_color=match.hex_color,
            is_significant=match.is_significant,
            match_status=MatchStatus.AUTOMATED_MATCH,
        )


@frozen
class SequenceReport:
    sequence_id: None | str
    sequence_length: int
    sequence_md5: str
    status: SequenceStatus

    @classmethod
    def from_info(
        cls,
        accession: None | SequenceAccession,
        info: ChainInfo,
        status: SequenceStatus,
    ) -> SequenceReport:
        acc = None
        if accession:
            acc = str(accession)
        return cls(
            sequence_id=acc,
            sequence_length=info.sequence_length,
            sequence_md5=info.sequence_md5_hash,
            status=status,
        )


@frozen
class CandidateReport:
    chain_id: PdbChainId
    match_report: MatchReport
    chain: None | ChainReport
    sequence: None | SequenceReport

    @classmethod
    def from_action(
        cls, action: CandidateAction, basepairing_status: None | BasepairingStatus
    ) -> CandidateReport:
        match action:
            case StructurelessCandidate(match=match):
                return CandidateReport(
                    chain_id=match.chain_id,
                    match_report=MatchReport.from_match(match),
                    chain=None,
                    sequence=None,
                )
            case SkippedCandidate(candidate=candidate):
                return CandidateReport(
                    chain_id=candidate.match.chain_id,
                    match_report=MatchReport.from_match(candidate.match),
                    chain=ChainReport.from_info(
                        candidate.chain_info,
                        ChainStatus.SKIPPED,
                        basepairing_status,
                    ),
                    sequence=None,
                )
            case AlreadyPresentCandidate(accession=accession, candidate=candidate):
                return CandidateReport(
                    chain_id=candidate.match.chain_id,
                    match_report=MatchReport.from_match(candidate.match),
                    chain=ChainReport.from_info(
                        candidate.chain_info,
                        ChainStatus.ALREADY_PRESENT,
                        basepairing_status,
                    ),
                    sequence=SequenceReport.from_info(
                        accession, candidate.chain_info, SequenceStatus.KNOWN_SEQUENCE
                    ),
                )
            case AlignCandidateSequence(candidate=candidate, accession=accession):
                return CandidateReport(
                    chain_id=candidate.match.chain_id,
                    match_report=MatchReport.from_match(candidate.match),
                    chain=ChainReport.from_info(
                        candidate.chain_info,
                        ChainStatus.NEW_CHAIN,
                        basepairing_status,
                    ),
                    sequence=SequenceReport.from_info(
                        accession, candidate.chain_info, SequenceStatus.NEW_SEQUENCE
                    ),
                )
            case AddCandidateStructure(candidate=candidate, accession=accession):
                return CandidateReport(
                    chain_id=candidate.match.chain_id,
                    match_report=MatchReport.from_match(candidate.match),
                    chain=ChainReport.from_info(
                        candidate.chain_info,
                        ChainStatus.NEW_CHAIN,
                        basepairing_status,
                    ),
                    sequence=SequenceReport.from_info(
                        accession, candidate.chain_info, SequenceStatus.KNOWN_SEQUENCE
                    ),
                )
            case _:
                assert_never(action)


@frozen
class AlignmentReport:
    num_columns_before: int
    num_columns_after: int


@frozen
class FamilyReport:
    family_id: str
    family_name: str
    rna_type: str
    status: FamilyStatus
    alignment_info: None | AlignmentReport
    candidates: list[CandidateReport]

    @classmethod
    def from_action(cls, family_action: Action) -> FamilyReport:
        candidate_reports = []
        family = family_action.family
        candidate_actions = family_action.actions
        for candidate_action in candidate_actions:
            bp_status = None
            if basepairing := candidate_action.basepairing:
                bp_status = BasepairingStatus.NEW_STRUCTURE
                if family.known_basepairing(basepairing):
                    bp_status = BasepairingStatus.KNOWN_STRUCTURE
            candidate_reports.append(
                CandidateReport.from_action(candidate_action, bp_status)
            )

        return FamilyReport(
            family_id=family_action.family.rfam_accession,
            family_name=family_action.family.info.description,
            rna_type=family_action.family.info.rna_type,
            status=FamilyStatus.status_for(family_action),
            alignment_info=None,
            candidates=candidate_reports,
        )


@frozen
class Report:
    families: list[FamilyReport]

    @classmethod
    def empty(cls) -> Report:
        return cls(families=[])

    def track(self, action: Action):
        self.families.append(FamilyReport.from_action(action))

    def write(self, handle: ty.TextIO):
        logger.trace("Writing report to {}", handle.name)
        converter = Converter()

        def field_converters(cls):
            converters = {}
            for field in fields(cls):
                assert isinstance(field, Attribute)
                converters[field.name] = override(rename=to_camel_case(field.name))
            return converters

        def to_camel_case_unstructure(cls):
            return make_dict_unstructure_fn(cls, converter, **field_converters(cls))

        def to_camel_case_structure(cls):
            return make_dict_structure_fn(cls, converter, **field_converters(cls))

        converter.register_unstructure_hook_factory(has, to_camel_case_unstructure)
        converter.register_structure_hook_factory(has, to_camel_case_structure)
        raw = converter.unstructure(self)
        json.dump(raw["families"], handle)

    def __bool__(self) -> bool:
        return bool(self.families)
