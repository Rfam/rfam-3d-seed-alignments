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
    """This models what kind of accession each SequenceAccession is. This can
    be useful to tracking and figuring out what to do with each Accession.
    """

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
        """Build a SequenceAccession from a given string. If given a specific
        kind this will be the default kind, but will be overridden if the name
        matches some known patterns.

        >>> SequenceAccession.build("URS00000001")
        SequenceAccession(raw='URS00000001', kind=<AccessionKind.RNACENTRAL: 'rnacentral'>)
        >>> SequenceAccession.build("URS00000001", kind=AccessionKind.OTHER)
        SequenceAccession(raw='URS00000001', kind=<AccessionKind.RNACENTRAL: 'rnacentral'>)
        >>> SequenceAccession.build("1S72", kind=AccessionKind.OTHER)
        SequenceAccession(raw='1S72', kind=<AccessionKind.PDB_CHAIN: 'pdb_chain'>)
        >>> SequenceAccession.build("NR_1")
        SequenceAccession(raw='NR_1', kind=<AccessionKind.OTHER: 'other'>)
        >>> SequenceAccession.build("NR_1", kind=AccessionKind.HASH)
        SequenceAccession(raw='NR_1', kind=<AccessionKind.HASH: 'hash'>)
        """

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
        """Get the accession part of this SequenceAccession.

        >>> SequenceAccession.build("NR_1/1-10").sequence_id
        'NR_1'
        >>> SequenceAccession.build("1S72").sequence_id
        '1S72'
        """
        return self.raw.split("/", 1)[0]

    @property
    def range(self) -> None | ty.Tuple[int, int]:
        """Compute the start-stop range for this SequenceAccession, if one
        exists. This does not change the coordinate system provided.

        >>> SequenceAccession.build("NR_1/1-10").range
        (1, 10)
        >>> SequenceAccession.build("1S72").range
        >>> SequenceAccession.build("URS000001/100-20").range
        (100, 20)
        """
        if "/" not in self.raw:
            return None
        parts = self.raw.split("/", 1)
        if "-" in parts[1]:
            start, stop = parts[1].split("-", 1)
            return (int(start), int(stop))
        return None

    def matches_sequence(self, accession: str | SequenceAccession) -> bool:
        """Check if this SequenceAccession has the same accession as the given
        accession. If the given accession is a SequenceAccession, then the
        sequence_id are compared

        >>> SequenceAccession.build("NR_1/1-10").matches_sequence("NR_1")
        True
        >>> SequenceAccession.build("1S72").matches_sequence("1S72")
        True
        >>> SequenceAccession.build("URS000001/100-20").matches_sequence("URS000001")
        True
        >>> SequenceAccession.build("URS000001/100-20").matches_sequence("NR_1")
        False
        >>> SequenceAccession.build("NR_1/1-10").matches_sequence(SequenceAccession.build("NR_1/20-30"))
        True
        """

        if isinstance(accession, SequenceAccession):
            return accession.sequence_id == self.sequence_id
        return self.sequence_id == accession

    def __str__(self) -> str:
        return self.raw
