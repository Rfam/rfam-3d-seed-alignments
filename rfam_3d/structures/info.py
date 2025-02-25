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

"""This contains the logic for fetching information about a structure and the
chains in it.
"""

from __future__ import annotations

from attrs import frozen
from Bio.Seq import Seq
from diskcache import Cache
from loguru import logger

from rfam_3d.rfam import alignments
from rfam_3d.rfam.accessions import AccessionKind, SequenceAccession
from rfam_3d.rnacentral import RnacentralApi, RnacentralMapping
from rfam_3d.structures.pdbe import ExperimentalInfo, PdbeApi
from rfam_3d.structures.rna3dhub import Basepairing, Rna3dHubApi
from rfam_3d.structures.structure_id import PdbChainId


@frozen
class ChainInfo:
    chain_id: PdbChainId
    rnacentral_id: str | None
    info: ExperimentalInfo
    sequence: str
    sequence_md5_hash: str
    basepairing: Basepairing

    def seq(self) -> Seq:
        return Seq(self.sequence)

    @property
    def sequence_length(self) -> int:
        return len(self.sequence)

    def unique_accession(self) -> SequenceAccession:
        if self.rnacentral_id:
            return SequenceAccession(
                raw=f"{self.rnacentral_id}/1-{len(self.sequence)}",
                kind=AccessionKind.RNACENTRAL,
            )
        return SequenceAccession(
            raw=self.sequence_md5_hash,
            kind=AccessionKind.HASH,
        )

    @property
    def has_pk(self) -> bool:
        return self.basepairing.has_pk

    @property
    def has_basepairs(self) -> bool:
        return self.basepairing.has_basepairs


@frozen
class StructureLookup:
    pdbe: PdbeApi
    rna3dhub: Rna3dHubApi
    rnacentral: RnacentralApi
    mapping: RnacentralMapping

    @classmethod
    def build(
        cls,
        cache: Cache,
        mapping: RnacentralMapping,
        pdbe_per_second=10,
        rna3dhub_per_second=10,
        rnacentral_per_second=10,
    ) -> StructureLookup:
        pdbe = PdbeApi.build(cache, per_second=pdbe_per_second)
        rna3dhub = Rna3dHubApi.build(cache, per_second=rna3dhub_per_second)
        rnacentral = RnacentralApi.build(cache, per_second=rnacentral_per_second)
        return StructureLookup(
            pdbe=pdbe, rna3dhub=rna3dhub, rnacentral=rnacentral, mapping=mapping
        )

    def fetch_structure_info(self, chain_id: PdbChainId) -> None | ChainInfo:
        logger.info("Fetching structure info for {}", chain_id)
        basepairing = self.rna3dhub.basepairing(chain_id)
        if not basepairing:
            logger.warning("Found no basepairs for {}", chain_id)
            return None

        seq = Seq(basepairing.sequence)
        info = self.pdbe.experiment(chain_id.pdb_id)
        return ChainInfo(
            chain_id=chain_id,
            rnacentral_id=self.mapping.urs_taxid(chain_id),
            info=info,
            sequence=basepairing.sequence,
            sequence_md5_hash=alignments.normalized_hash(seq),
            basepairing=basepairing,
        )
