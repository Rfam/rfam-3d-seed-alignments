# -*- coding: utf-8 -*-

# Copyright [2009-2025] EMBL-European Bioinformatics Institute
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

import csv
import hashlib
import json
import typing as ty

import requests
from attrs import frozen
from Bio.Seq import Seq
from cattr import GenConverter
from diskcache import Cache
from loguru import logger
from requests import Session
from requests_ratelimiter import LimiterAdapter

from rfam_3d.structures.structure_id import PdbChainId

RNACENTRAL_URL = "https://rnacentral.org/api/v1/rna"


@frozen
class RnacentralMapping:
    mapping: dict[PdbChainId, str]

    @classmethod
    def build(cls, mapping: dict[PdbChainId, str]) -> RnacentralMapping:
        """Create a new RnacentralMapping with the given mapping."""
        return cls(mapping=mapping)

    @classmethod
    def parse_tsv(cls, handle: ty.IO) -> RnacentralMapping:
        """Parse the given TSV file from RNAcentral to produce the mapping.
        This is expected to be the mapping file from RNAcentral's API.
        """

        reader = csv.DictReader(
            handle, delimiter="\t", fieldnames=["urs", "db", "key", "taxid", "rna_type"]
        )
        mapping = {}
        for row in reader:
            key = PdbChainId.build(row["key"])
            mapping[key] = f"{row['urs']}_{row['taxid']}"
        return RnacentralMapping.build(mapping)

    @classmethod
    def from_handle(cls, handle: ty.TextIO) -> RnacentralMapping:
        """Load the data in the given TextIO handle to return a mapping. This
        is the counter part to the to_handle method.
        """

        c = GenConverter()
        c.register_structure_hook(PdbChainId, lambda v, _: PdbChainId.build(v))
        raw = json.load(handle)
        return c.structure(raw, RnacentralMapping)

    def urs_taxid(self, chain_id: PdbChainId) -> str | None:
        """Find the URS_taxid for a given PdbChainId, if it exists."""
        return self.mapping.get(chain_id, None)

    def to_handle(self, handle: ty.TextIO):
        """Write this mapping to the given handle. This is stored as a JSON
        object and a new line is appened.
        """

        c = GenConverter()
        c.register_unstructure_hook(PdbChainId, lambda pid: str(pid))
        raw = c.unstructure(self)
        json.dump(raw, handle)
        handle.write("\n")


@frozen
class RnacentralApi:
    """This class wraps talking to the RNAcentral API and caching the results."""

    cache: Cache
    session: Session

    @classmethod
    def with_session(
        cls, cache: Cache, session: Session, per_second=10
    ) -> RnacentralApi:
        rnacentral_limiter = LimiterAdapter(per_second=per_second)
        session.mount("https://rnacentral.org/", rnacentral_limiter)
        return RnacentralApi(cache, session)

    @classmethod
    def build(cls, cache: Cache, per_second=10) -> RnacentralApi:
        """Build a new wrapper from the RNAcentral API."""
        session = Session()
        return RnacentralApi.with_session(cache, session, per_second=per_second)

    def rnacentral_id(self, sequence: Seq) -> str | None:
        """Fetch the URS, if it exists, for the given sequence. This fetches
        the generic URS and uses the given taxid as the taxid component. The
        taxid is not validated, and is just appended. If the sequence is not
        known then None is returned.

        >>> api = RnacentralApi.build(Cache())
        >>> api.rnacentral_id(Seq("UAGAUCUUUGACUCUGGCAGUCUCCAGG"))
        "URS000075E7DE_9606"
        >>> api.rnacentral_id(Seq("TAGATCTTTGACTCTGGCAGTCTCCAGG"))
        "URS000075E7DE_9606"
        >>> api.rnacentral_id(Seq("AAGUAGUUGGUUUGUAUGAGAUGGUU"))
        "URS000075B58F_9606"
        >>> api.rnacentral_id(Seq("NNNNNNNNNNNNNNNNNNNNNNNNNN"))
        None
        """

        logger.debug("Fetching RNAcentral id")
        seq = sequence.replace("U", "T")
        seq = str(seq).encode("ascii")
        m = hashlib.md5()
        m.update(seq)
        md5 = m.hexdigest()
        key = f"rnacentral/{md5}"
        if urs := self.cache.get(key):
            logger.trace("Using existing value")
            assert isinstance(urs, str)
            return urs

        response = requests.get(RNACENTRAL_URL, params={"md5": md5})
        response.raise_for_status()
        data = response.json()
        if not data["results"]:
            logger.warning("No RNAcentral id found")
            return None

        logger.trace("URS found")
        urs = data["results"][0]["rnacentral_id"]
        self.cache.set(key, urs)
        return urs
