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

"""This module representes Candidate, which are a Match which could be added to an
alignment, alongside their associated metdata."""

from __future__ import annotations

from attrs import frozen
from Bio.Seq import Seq

from rfam_3d.rfam.alignments import Alignment
from rfam_3d.rfam.api import FamilyInfo
from rfam_3d.rfam.matches import Match
from rfam_3d.structures.info import ChainInfo
from rfam_3d.structures.structure_id import PdbChainId, PdbId


@frozen
class Candidate:
    match: Match
    chain_info: ChainInfo

    @property
    def chain_id(self) -> PdbChainId:
        return self.chain_info.chain_id

    @property
    def pdb_id(self) -> PdbId:
        return self.chain_info.chain_id.pdb_id

    # FIXME: Rework this so it strips out the extra gaps in the alignment. This
    #        will have to parse the aligned sequence and remove all columns in
    #        the basepairing are gaps in the sequence. This should produce the
    #        shortened pairing
    # TODO: Consider adding a parameter to extract the sequence pairing only.
    # FIXME: This should use 'pseudoknotted not complete"
    @property
    def complete_basepairing(self) -> str:
        return self.chain_info.basepairing.psuedoknotted()

    def same_sequence(self, candidate: Candidate) -> bool:
        """Check if this candidate and another have the same sequence by
        checking their MD5 hashes.
        """
        return (
            self.chain_info.sequence_md5_hash == candidate.chain_info.sequence_md5_hash
        )

    def full_sequence(self) -> Seq:
        return self.chain_info.seq()

    def truncated_sequence(self) -> Seq:
        assert False


@frozen
class CandidateFamily:
    rfam_accession: str
    info: FamilyInfo
    alignment: Alignment
    candidates: list[Candidate]
