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

"""This module contains a wrapper around PDBe's API.
"""

from __future__ import annotations

from attrs import field, frozen
from diskcache import Cache
from loguru import logger
from requests import Session
from requests_ratelimiter import LimiterAdapter

from rfam_3d.structures.structure_id import PdbId
from rfam_3d.utils import asupper

PDBE_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/"


@frozen
class ExperimentalInfo:
    """A summary of some experimental information about a structure."""

    pdb_id: PdbId
    resolution: None | float
    method: str = field(converter=asupper)


@frozen
class PdbeApi:
    """The wrapper around PDBe's API. This will cache results from requests."""

    cache: Cache
    session: Session

    @classmethod
    def with_session(cls, session: Session, cache: Cache, per_second=10) -> PdbeApi:
        """Create a new PdbeApi with the given session cache and per_second rate limit.
        The sesion will be modified to have a rate limiter for the PDBe api.
        """

        limiter = LimiterAdapter(per_second=per_second)
        session.mount("https://www.ebi.ac.uk/pdbe/", limiter)
        return PdbeApi(cache=cache, session=session)

    @classmethod
    def build(cls, cache: Cache, per_second=10) -> PdbeApi:
        """Create a new PdbeApi with the given cache and per_second rate limit.
        This will create a new Session for the API.
        """
        return PdbeApi.with_session(Session(), cache, per_second)

    def experiment_info(self, pdb_id: PdbId) -> ExperimentalInfo:
        """Query the PDBe API for the experiment information about the given
        structure.
        """

        pid = pdb_id.pdb_id
        logger.debug("Fetching structure info for {}", pid)
        key = f"pdbe/info-{pid}"
        if value := self.cache.get(key):
            logger.trace("Using existing value")
            assert isinstance(value, ExperimentalInfo)
            return value

        response = self.session.get(PDBE_URL + pid)
        response.raise_for_status()
        data = response.json()
        if pid not in data:
            raise ValueError(f"No structure info found for {pid}")

        raw = data[pid][0]
        value = ExperimentalInfo(
            pdb_id=pdb_id,
            resolution=raw.get("resolution", None),
            method=raw["experimental_method"],
        )
        self.cache.set(key, value)
        return value
