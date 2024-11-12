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
import itertools as it
import json
import operator as op
import typing as ty

from attr import frozen
from attrs import Attribute, fields, has
from cattrs import Converter
from cattrs.gen import make_dict_structure_fn, make_dict_unstructure_fn, override
from loguru import logger

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
class ChainStatus(enum.Enum):
    ERROR = "error"
    SKIPPED = "skipped"
    ALREADY_PRESENT = "already_present"

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
    SKIPPED = "skipped"
    COMPLETE = "complete"
    CURATE = "curate"
    INCOMPLETE = "incomplete"
    NO_VALID = "no_matches"

    @classmethod
    def status_for(cls, action: Action) -> FamilyStatus:
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
    has_base_pairs: bool
    has_pseudoknots: bool
    new_basepair_count: int
    removed_basepair_count: int


@frozen
class StructureReport:
    pdb_id: str
    resolution: None | float
    method: str
    chains: list[ChainReport]


@frozen
class SequenceReport:
    sequence_id: str
    sequence_length: int
    status: SequenceStatus
    structures: list[StructureReport]


@frozen
class FamilyReport:
    family_id: str
    family_name: str
    rna_type: str
    # num_columns: int
    status: FamilyStatus
    sequences: list[SequenceReport]

    @classmethod
    def from_action(cls, family_action: Action) -> FamilyReport:
        sequence_reports = []
        key = lambda a: getattr(a, "accession", "Error")
        pdb_key = op.attrgetter("pdb_id")
        chain_key = op.attrgetter("chain_id")

        def grouped(iterable: ty.Iterable, key) -> ty.Iterable:
            ordered = sorted(iterable, key=key)
            return it.groupby(ordered, key)

        for accession, per_sequence in grouped(family_action.actions, key):
            structure_reports = []
            per_sequence = list(per_sequence)
            for pdb_id, per_structure in grouped(per_sequence, pdb_key):
                chain_reports = []
                per_structure = list(per_structure)
                for chain_id, per_chain in grouped(per_structure, chain_key):
                    per_chain = list(per_chain)
                    pass

            # for _, candidates in it.groupby(structures, chain_key):
            #     candidates = list(candidates)
            #     candidate = None
            #     if len(candidates) == 1:
            #         candidate = candidates[0]
            #     elif len(candidates) == 2:
            #         assert any(
            #             isinstance(c, AlignCandidateSequence) for c in candidates
            #         ) and any(isinstance(c, AddCandidateStructure) for c in candidates)
            #         candidate = next(
            #             c for c in candidates if isinstance(c, AlignCandidateSequence)
            #         )
            #     else:
            #         logger.error("Saw too many candidates {}", len(candidates))
            #         raise ValueError(
            #             "Should never have more than 2 candidates for a chain"
            #         )
            #
            #     has_bp = False
            #     has_pk = False
            #     sequence_id = None
            #     sequence_length = None
            #     if not isinstance(
            #         candidate,
            #         (AlreadyPresentCandidate, StructurelessCandidate, SkippedCandidate),
            #     ):
            #         has_bp = candidate.candidate.chain_info.has_basepairs
            #         has_pk = candidate.candidate.chain_info.has_pk
            #         sequence_id = candidate.accession.raw
            #         sequence_length = candidate.candidate.chain_info.sequence_length
            #
            #     if isinstance(candidate, AlreadyPresentCandidate):
            #         sequence_id = candidate.accession.raw
            #
            #     structure_reports.append(
            #         ChainReport(
            #             chain_id=candidate.chain_id.chain_id,
            #             status=ChainStatus.status_for(candidate),
            #             has_base_pairs=has_bp,
            #             has_pseudoknots=has_pk,
            #             is_unique_sequence=False,
            #             sequence_id=sequence_id,
            #             sequence_length=sequence_length,
            #         )
            #     )

            method = "Failed to find method"
            resolution = None
            if not isinstance(structures[0], StructurelessCandidate):
                resolution = structures[0].candidate.chain_info.info.resolution
                method = structures[0].candidate.chain_info.info.method

            structure_reports.append(
                StructureReport(
                    pdb_id=pdb_id.pdb_id,
                    method=method,
                    resolution=resolution,
                    chains=structure_reports,
                )
            )

        return FamilyReport(
            family_id=family_action.family.rfam_accession,
            family_name=family_action.family.info.description,
            status=FamilyStatus.status_for(family_action),
            rna_type=family_action.family.info.rna_type,
            sequences=sequence_reports,
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
