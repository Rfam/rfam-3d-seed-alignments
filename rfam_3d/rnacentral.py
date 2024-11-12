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

import hashlib

import requests
from attrs import frozen
from Bio.Seq import Seq
from diskcache import Cache
from loguru import logger
from requests import Session
from requests_ratelimiter import LimiterAdapter

RNACENTRAL_URL = "https://rnacentral.org/api/v1/rna"


@frozen
class RnacentralApi:
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
        """Fetch the URS, if it exists, for the given sequence."""
        logger.debug("Fetching RNAcentral id")
        seq = sequence.replace("U", "T")
        seq = str(seq).encode("ascii")
        m = hashlib.md5()
        m.update(seq)
        md5 = m.hexdigest()
        key = f"rnacentral/{md5}"
        if value := self.cache.get(key):
            logger.trace("Using existing value")
            assert isinstance(value, str)
            return value

        response = requests.get(RNACENTRAL_URL, params={"md5": md5})
        response.raise_for_status()
        data = response.json()
        if not data["results"]:
            logger.warning("No RNAcentral id found")
            return None

        logger.debug("URS found")
        value = data["results"][0]["rnacentral_id"]
        self.cache.set(key, value)
        return value
