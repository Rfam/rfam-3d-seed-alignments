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

from __future__ import annotations

import typing as ty

from attr import frozen
from Bio.Align import MultipleSeqAlignment
from loguru import logger

from rfam_3d.structures.info import ChainInfo


class MissingSequences(Exception):
    """This is raised if the FamilyActions is asked to compute RF lines for missing sequences."""


@frozen
class SecondaryStructureAnnotation:
    info: ChainInfo
    line: str


@frozen
class StructureActions:
    structures_to_add: dict[str, list[ChainInfo]]

    @classmethod
    def empty(cls) -> StructureActions:
        return cls(structures_to_add={})

    def annotate_structure(self, accession: str, pdb_info: ChainInfo):
        """Add the structure"""
        current = self.structures_to_add.get(accession, list())
        if pdb_info not in current:
            current.append(pdb_info)
        self.structures_to_add[accession] = current

    def secondary_structure_annotations(
        self, alignment: MultipleSeqAlignment
    ) -> ty.Dict[str, ty.List[SecondaryStructureAnnotation]]:
        lines: ty.Dict[str, ty.List[SecondaryStructureAnnotation]] = {}
        for record in alignment:
            current = lines.get(record.id, [])
            for info in self.structures_to_add.get(record.id, []):
                line = []
                secondary = list(info.basepairing.psuedoknotted())
                for char in str(record.seq):
                    if char in ["-", "."]:
                        line.append(".")
                    else:
                        line.append(secondary.pop(0))
                assert not secondary, "Did not process all 2D"
                assert len(line) == len(record.seq), "Mismatch secondary and 2D"
                annotation = SecondaryStructureAnnotation(info=info, line="".join(line))
                current.append(annotation)
            if current:
                lines[record.id] = current
        missing = set(self.structures_to_add.keys()) - lines.keys()
        if missing:
            logger.error(
                "Missed computing secondary structures for {}", ",".join(missing)
            )
            raise MissingSequences(missing)
        return lines

    def __len__(self) -> int:
        return len(self.structures_to_add)

    def __bool__(self) -> bool:
        return bool(self.structures_to_add)
