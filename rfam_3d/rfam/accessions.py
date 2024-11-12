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
import re
import typing as ty

from attrs import frozen


@enum.unique
class AccessionKind(enum.Enum):
    RNACENTRAL = "rnacentral"
    PDB_CHAIN = "pdb_chain"
    HASH = "hash"
    OTHER = "other"


@frozen
class SequenceAccession:
    """This represents an accession within a Rfam alignment. These accession
    are slightly special since they are generally of the form:
    '<seq_id>/<start>-<stop>'. This handles class is meant to do things
    like checking if a given sequence id matches this accession.
    """

    raw: str
    kind: AccessionKind

    @classmethod
    def build(cls, raw: str, kind=None) -> SequenceAccession:
        found = AccessionKind.OTHER
        if kind:
            found = kind
        if raw.startswith("URS"):
            found = AccessionKind.RNACENTRAL
        elif re.match(r"^[1-9]\w{3}", raw):
            found = AccessionKind.PDB_CHAIN
        return cls(raw=raw, kind=found)

    @property
    def sequence_id(self) -> str:
        return self.raw.split("/", 1)[0]

    @property
    def range(self) -> None | ty.Tuple[int, int]:
        parts = self.raw.split("/", 1)
        if "-" in parts[1]:
            start, stop = parts[1].split("-", 1)
            return (int(start), int(stop))
        return None

    def matches_sequence(self, seq_id: str) -> bool:
        return self.sequence_id == seq_id
