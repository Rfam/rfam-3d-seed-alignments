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

import typing as ty

from attrs import field, frozen
from diskcache import Cache
from loguru import logger
from requests import Session
from requests_ratelimiter import LimiterAdapter

from rfam_3d.structures.structure_id import PdbChainId, PdbId
from rfam_3d.utils import asupper

PDB_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/{name}/{pdb_id}"
# EXPERIMENT_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/"
# MOLECULE_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/"


@frozen
class ExperimentalInfo:
    """A summary of some experimental information about a structure."""

    pdb_id: PdbId
    resolution: None | float
    method: str = field(converter=asupper)


@frozen
class MoleculeInfo:
    """A summary of the molecules found in a structure."""

    chain_id: PdbChainId
    taxid: None | int


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

    def __fetch__(self, name: str, pdb_id: str) -> ty.Any:
        url = PDB_URL.format(name=name, pdb_id=pdb_id)
        logger.trace("Requesting {} info for {}", name, pdb_id)
        response = self.session.get(url)
        response.raise_for_status()
        data = response.json()
        if pdb_id not in data:
            raise ValueError(f"Request for {pdb_id} failed")
        return data[pdb_id]

    def molecules(self, chain_id: PdbChainId) -> MoleculeInfo:
        pid = chain_id.pdb_id.pdb_id
        logger.debug("Fetching experimental info for {}", pid)
        key = f"pdbe/molecule-{chain_id}"
        if value := self.cache.get(key):
            logger.trace("Using existing value for {}", pid)
            assert isinstance(value, MoleculeInfo)
            return value

        logger.debug("Fetching molecule info for {}", pid)
        raw = self.__fetch__("molecules", pid)
        chain_info = [c for c in raw if chain_id.chain_id in c["in_chains"]]
        if len(chain_info) == 0:
            raise ValueError(f"Failed to find information about {chain_id}")
        if len(chain_info) > 1:
            raise ValueError(f"Found duplicate information about {chain_id}")

        chain_info = chain_info[0]
        value = MoleculeInfo(
            chain_id=chain_id,
            taxid=chain_info["source"]["tax_id"],
        )
        self.cache.set(key, value)
        return value

    def experiment(self, pdb_id: PdbId) -> ExperimentalInfo:
        """Query the PDBe API for the experiment information about the given
        structure.
        """

        pid = pdb_id.pdb_id
        logger.debug("Fetching experimental info for {}", pid)
        key = f"pdbe/experiment-{pid}"
        if value := self.cache.get(key):
            logger.trace("Using existing value for {}", pid)
            assert isinstance(value, ExperimentalInfo)
            return value

        raw = self.__fetch__("experiment", pid)[0]
        value = ExperimentalInfo(
            pdb_id=pdb_id,
            resolution=raw.get("resolution", None),
            method=raw["experimental_method"],
        )
        self.cache.set(key, value)
        return value
