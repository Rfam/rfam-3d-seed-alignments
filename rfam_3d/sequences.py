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

"""A module with some utilties for working with sequences.

>>> normalize_sequence(Seq("ac-cg."))
Seq("ACCG")
>>> normalized_hash(Seq("ac-cg."))
"c3e65f1cfdc67810a569e30559d77178"
"""

import hashlib

from Bio.Seq import Seq


def normalize_sequence(sequence: Seq) -> Seq:
    """Clean up and normalize the sequences so they can be compared reilable.
    This strips all gap characters ('.', '-') and converts the sequence to
    upper case before comparing.

    >>> normalize_sequence(Seq("ac-cg."))
    Seq("ACCG")
    """
    return sequence.upper().replace(".", "").replace("-", "")


def normalized_hash(sequence: Seq) -> str:
    """Normalize the sequence, compute its MD5 hash and return the hexdigest.

    >>> normalized_hash(Seq("ac-cg."))
    "c3e65f1cfdc67810a569e30559d77178"
    """
    return hashlib.md5(bytes(normalize_sequence(sequence))).hexdigest()
