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

"""This module contains the logic for parsing an Rfam stockholm alignment into
an Alignment structure. The parsing is very simple and assumes that the
sequence is in the 'pfam' style, which each sequence on one line.
"""

from __future__ import annotations

import json
import re
import typing as ty

import cattrs
from attrs import frozen
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from loguru import logger

from rfam_3d.rfam.accessions import AccessionKind, SequenceAccession
from rfam_3d.sequences import normalized_hash
from rfam_3d.structures.structure_id import PdbChainId, PdbId

NAME_PATTERN = r"([1-9]\w{3})_(\w+)_SS$"


@frozen
class AlignmentSequence:
    """This represents the sequence that is part of an alignment. It doesn't
    allow reconstruction of the whole alignment it is just meant to capture the
    metadata about what is stored. Notably this does not store all annotations,
    just the annotated structures.

    :accession: The accession for the sequence.
    :sequence: The sequence in the alignment.
    :sequence_md5_hash: The MD5 hash of the normalized sequence.
    :annotated_structures: All per letter (GR lines) annotations for this
    sequence which are a secondary structures of a PDB structure.
    """

    accession: SequenceAccession
    sequence: str
    sequence_md5_hash: str
    annotated_structures: ty.List[PdbChainId]

    def has_annotation_for(self, chain: PdbChainId) -> bool:
        return any(c == chain for c in self.annotated_structures)


@frozen
class Alignment:
    """This is meant to represent what is part of an Alignment. This doesn't
    actually store the sequences, or need to just enough to track what is
    stored.

    :rfam_accession: The Rfam accession.
    :sequences: A list of the AlignmentSequences
    """

    rfam_accession: str
    sequences: ty.List[AlignmentSequence]

    def known_accession(self, acc: SequenceAccession) -> bool:
        """Check if any sequence has the given accession."""
        return any(s.accession.raw == acc.raw for s in self.sequences)

    def known_sequence(self, seq: Seq) -> None | SequenceAccession:
        """Check if the given sequence exists within this alignment already.
        This is done on the basis of normalized hashes, so gaps are ignored.
        """
        hash = normalized_hash(seq)
        for sequence in self.sequences:
            if sequence.sequence_md5_hash == hash:
                return sequence.accession
        return None

    def has_structure(self) -> bool:
        """Check if there is any aligned structure in this alignment."""
        return any(s.annotated_structures for s in self.sequences)

    def sequence_annotated_with(self, chain: PdbChainId) -> None | SequenceAccession:
        """Check if any sequence in the alignment has an annotation with the
        given PdbChainId.
        """
        for sequence in self.sequences:
            if sequence.has_annotation_for(chain):
                return sequence.accession
        return None


def parse_file(accession: str, alignment: ty.TextIO, format="stockholm") -> Alignment:
    """Extract all mapped PDB structures and the associated sequence
    accessions, which may be an RNAcentral id. PDB structures are detected as
    being mapped by looking for GR annotations for sequences.

    :accession: The Rfam accession for this alignment.
    :alignment: A file handle to the alignment file to parse.
    :format: The format of the alignment file.
    """

    sequences = []
    align = AlignIO.read(alignment, format)
    for sequence in align:
        assert isinstance(sequence, SeqRecord)
        logger.debug("Checking {} for mapped PDBs", sequence.id)
        structures_mapped: list[PdbChainId] = []
        for name, _ in sequence.letter_annotations.items():
            if not name.startswith("GR:"):
                logger.debug("Annotation {} is not a GR line", name)
                continue

            name = name[3:]
            if match := re.match(NAME_PATTERN, name):
                chain = PdbChainId(
                    pdb_id=PdbId(pdb_id=match.group(1)), chain_id=match.group(2)
                )
                logger.debug(
                    "Accession {} has an associated PDB {}", sequence.id, chain
                )
                structures_mapped.append(chain)
            else:
                logger.debug("Annotation {} is not a PDB", name)

        # Try to find the most specific accession kind we can.
        kind = None
        if len(structures_mapped) == 1 and sequence.id == str(structures_mapped[0]):
            kind = AccessionKind.PDB_CHAIN

        assert sequence.id, "Sequence must have an id"
        sequences.append(
            AlignmentSequence(
                accession=SequenceAccession.build(sequence.id, kind=kind),
                sequence=str(sequence.seq),
                sequence_md5_hash=normalized_hash(sequence.seq),
                annotated_structures=structures_mapped,
            )
        )

    return Alignment(
        rfam_accession=accession,
        sequences=sequences,
    )


def load_jsonl(handle: ty.TextIO) -> list[Alignment]:
    logger.trace("Loading alignment info from {}", handle.name)
    alignments = []
    for line in handle:
        alignment = cattrs.structure(json.loads(line), Alignment)
        alignments.append(alignment)
    logger.trace("Found {} alignments", len(alignments))
    return alignments
