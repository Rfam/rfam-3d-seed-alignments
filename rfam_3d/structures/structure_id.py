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

"""A modules which contains some classes meant to sever as identifiers for
structures. These classes are used to ensure case-insensitivity as needed and
to make sure which strings are what does not get confused.
"""

from attrs import field, frozen
from attrs.validators import instance_of

from rfam_3d.utils import aslower


@frozen(order=True)
class PdbId:
    """This models a PDB id. These are case-insensitive identifiers, generally
    4 characters starting with '1'.

    :pdb_id: The PDB id which will always be forced to lower case.
    """

    pdb_id: str = field(converter=aslower)

    def __str__(self):
        return f"{self.pdb_id}"


@frozen(order=True)
class PdbChainId:
    """This models a combination of a PDB id and chain ids. PDB ids are treated
    as PdbId to handle case-insensitivity, while a chain_id is case sensitive.

    :pdb_id: A PdbId
    :chain_id: A chain id in a structure, case sensitive.
    """

    pdb_id: PdbId = field(validator=instance_of(PdbId))
    chain_id: str

    def secondary_structure_id(self) -> str:
        """Create an the id that is used in a GR line of a stockholm file to
        indicate this is a secondary structure annotation. While pdb ids are
        stored in lowercase, they must be written upper case here.

        >>> PdbChainId(PdbId('1S72'), 'A').secondary_structure_id()
        '1S72_A_SS'
        >>> PdbChainId(PdbId('1s72'), 'A').secondary_structure_id()
        '1S72_A_SS'
        """
        pid = str(self.pdb_id).upper()
        return f"{pid}_{self.chain_id}_SS"

    def __str__(self) -> str:
        return f"{self.pdb_id}_{self.chain_id}"
