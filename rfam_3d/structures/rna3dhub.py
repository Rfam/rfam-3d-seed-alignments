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

import json
import typing as ty

import cattrs
from attrs import frozen
from diskcache import Cache
from loguru import logger
from requests import Session
from requests_ratelimiter import LimiterAdapter

from rfam_3d.structures.structure_id import PdbChainId

RNA3DHUB_URL = "https://rna.bgsu.edu/rna3dhub/rest/getSequenceBasePairs"


@frozen
class Rna3DHubAnnotation:
    seq_id1: int
    nt1: str
    bp: str
    seq_id2: int
    nt2: str
    crossing: int


@frozen
class BasepairAnnotations:
    pdb_id: str
    chain_id: str
    sequence: str
    annotations: ty.List[Rna3DHubAnnotation]

    def dot_bracket(self) -> str:
        structure = list("." * len(self.sequence))
        for annotation in self.annotations:
            # TODO: Check chains are the same for nt1 and nt2
            if annotation.bp == "cWW":
                structure[annotation.seq_id1 - 1] = "("
                structure[annotation.seq_id2 - 1] = ")"
        return "".join(structure)


@frozen
class Basepairing:
    chain_id: PdbChainId
    sequence: str
    _nested: str
    _complete: str

    def psuedoknotted(self) -> str:
        final = []
        assert len(self._nested) == len(
            self._complete
        ), "Basepairs have incorrect lengths"
        for i, character in enumerate(self._complete):
            if character in ["(", ")"] and self._nested[i] == ".":
                if character == "(":
                    final.append("{")
                elif character == ")":
                    final.append("}")
            else:
                final.append(character)
        assert len(final) == len(self._complete)
        return "".join(final)

    @property
    def has_pk(self):
        return "{" in self.psuedoknotted()

    @property
    def has_basepairs(self):
        return "(" in self.psuedoknotted()


@frozen
class Rna3dHubApi:
    cache: Cache
    session: Session

    @classmethod
    def build(cls, cache: Cache, per_second=10) -> Rna3dHubApi:
        """Build a new RNA3dHubApi with the given cache and per_second rate
        limit.
        """

        session = Session()
        bgsu_limiter = LimiterAdapter(per_second=per_second)
        session.mount("http://rna.bgsu.edu/", bgsu_limiter)
        session.mount("https://rna.bgsu.edu/", bgsu_limiter)
        return Rna3dHubApi(cache=cache, session=session)

    def basepair_annotations(
        self, chain_id: PdbChainId, nested: bool
    ) -> None | BasepairAnnotations:
        """Query the getSequenceBasePairs API endpoint for the basepair
        annotations.
        """

        if chain_id.pdb_id.pdb_id.lower() == "6t7t":
            return None

        key = f"rna3dhub/{chain_id.pdb_id}_{chain_id.chain_id}_{nested}"
        if value := self.cache.get(key):
            logger.info("Using existing data for {} (nested: {})", chain_id, nested)
            if isinstance(value, BasepairAnnotations):
                return value
            logger.error("Invalid cached data for {}", key)

        logger.debug(
            "Fetching RNA3DHub basepairs for {} (nested: {})", chain_id, nested
        )
        response = self.session.get(
            RNA3DHUB_URL,
            params={
                "pdb_id": chain_id.pdb_id.pdb_id,
                "chain": chain_id.chain_id,
                "only_nested": str(nested),
            },
        )
        response.raise_for_status()

        try:
            # So it seems that the debug option is on (as of 2024) so that
            # means the response has some HTML elements and comments in it. This
            # always seems to be of the same form (comment, data, comment,
            # rest) so we split the text into lines and take the second line
            # and hope for the best.
            text = response.text
            if text.startswith("<!-- DEBUG-VIEW START"):
                text = response.text.splitlines()[1]
            raw = json.loads(text)
        except json.JSONDecodeError:
            logger.info(
                "Did not get valid JSON data for {} (nested: {}), likely no sequence, skipping",
                chain_id,
                nested,
            )
            return None

        if raw.get("sequence", None) == "No sequence was found for the given id":
            logger.debug(
                "Did not find an RNA3DHub sequence for {} (nested: {})",
                chain_id,
                nested,
            )
            return None

        value = cattrs.structure(raw, BasepairAnnotations)
        self.cache.set(key, value)
        return value

    def basepairing(self, chain_id: PdbChainId) -> None | Basepairing:
        """Fetch all Basepairing information."""
        logger.info("Fetching structure info for {}", chain_id)
        nested_annotations = self.basepair_annotations(chain_id, True)
        complete_annotations = self.basepair_annotations(chain_id, False)

        if not nested_annotations or not complete_annotations:
            logger.warning("Did not find annotation data for {}", chain_id)
            return None

        if nested_annotations.sequence != complete_annotations.sequence:
            logger.warning("RNA3DHub is inconsitent about sequences")
            return None

        if nested_annotations.sequence == "No sequence was found for the given id":
            logger.warning("No sequence found")
            return None

        return Basepairing(
            chain_id=chain_id,
            sequence=nested_annotations.sequence,
            nested=nested_annotations.dot_bracket(),
            complete=complete_annotations.dot_bracket(),
        )
