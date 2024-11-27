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

"""This module contains RfamApi which is wrapper around Rfam's API.
"""

from __future__ import annotations

from attrs import frozen
from loguru import logger
from requests import Session
from requests_ratelimiter import LimiterAdapter


@frozen
class FamilyInfo:
    """This represents some the information about a family that can be feteched
    from the Rfam API.
    """

    id: str
    accession: str
    rna_type: str
    num_seed: int
    num_full: int
    description: str


@frozen
class RfamApi:
    """A wrapper aound Rfam's API. This does not cache any data."""

    session: Session

    @classmethod
    def with_session(cls, session: Session, per_second=10) -> RfamApi:
        """Build a new RfamApi object with the given session per_second rate
        limit. The session will be modified to have a rate limit for all Rfam
        urls.
        """
        limiter = LimiterAdapter(per_second=per_second)
        session.mount("http://rfam.org", limiter)
        session.mount("https://rfam.org", limiter)
        return cls(session=session)

    @classmethod
    def build(cls, per_second=10) -> RfamApi:
        """Build a new RfamApi object with the given per_second rate limit."""
        return cls.with_session(Session(), per_second=per_second)

    def info(self, accession: str) -> FamilyInfo:
        """Fetches the information about a given family from the Rfam API.

        >>> RfamApi.build().info("RF00008")
        FamilyInfo(id='Hammerhead_3', accession='RF00008', rna_type='Gene; ribozyme;', num_seed=85, num_full=750, description='Hammerhead ribozyme (type III)')
        """

        logger.debug("Fetching family info for {}", accession)
        response = self.session.get(
            "https://rfam.org/family/" + accession,
            params={"content-type": "application/json"},
        )
        response.raise_for_status()
        data = response.json()
        return FamilyInfo(
            id=data["rfam"]["id"],
            accession=data["rfam"]["acc"],
            rna_type=data["rfam"]["curation"]["type"],
            num_seed=int(data["rfam"]["curation"]["num_seed"]),
            num_full=int(data["rfam"]["curation"]["num_full"]),
            description=data["rfam"]["description"],
        )
