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

import json
import re
import sys
import typing as ty
from pathlib import Path

import cattrs
import click
import yaml
from Bio import AlignIO
from diskcache import Cache
from loguru import logger

from rfam_3d import disallow
from rfam_3d.report import Report
from rfam_3d.rfam import alignments, matches
from rfam_3d.rfam.family import build_families
from rfam_3d.structures.info import StructureLookup
from rfam_3d.update.actions import CompleteFamily, NoAcceptedMatches
from rfam_3d.update.planner import FamilyUpdate, Planner, SkippedFamily
from rfam_3d.update.structure_actions import StructureActions
from rfam_3d.utils import assert_never


@click.group()
@click.option(
    "--cache-path", default="./cache", envvar="RFAM_3D_CACHE", type=click.Path()
)
@click.option(
    "--log-level",
    default="INFO",
    help="Log level to use",
    type=click.Choice(
        ["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"],
        case_sensitive=False,
    ),
)
@click.pass_context
def main(ctx, cache_path="./cache", log_level="INFO"):
    logger.remove()
    logger.add(sys.stderr, level=log_level)

    ctx.ensure_object(dict)
    path = Path(cache_path)
    path.mkdir(parents=True, exist_ok=True)
    ctx.obj["cache_path"] = path


@main.command("compute-actions")
@click.option("--disallow-file", type=click.Path())
@click.argument("matches-file", type=click.File("r"))
@click.argument("alignments-file", type=click.File("r"))
@click.argument("report-file", type=click.File("w"))
@click.argument("sequences", type=click.Path())
@click.argument("info", type=click.Path())
@click.pass_context
def compute_actions_cmd(
    ctx,
    matches_file: ty.TextIO,
    alignments_file: ty.TextIO,
    sequences: str | Path,
    report_file: ty.TextIO,
    info: str | Path,
    disallow_file: None | str,
):
    """This computes the actions needed to take to update the given family.

    This reads the mapping file, which should contain all PDBs which map to the given
    accession, and determines what actions must be taken. It will write out a
    possibly empty file of sequences to align and a JSON file of PDB info. That
    info will include the secondary structures to add. If this file is not
    present that means there are no actions to take. There should never be a
    sequences file without an info file. If the PDB file is not present
    sequence file must not exist. No sequences, indicates the alignment does
    not need any new sequences aligned to it.

    Arguments\b
    ---------\b

    matches-file:\b
        A TSV file listing

    alignment-file:\b
        A JSONL file of that shows comes from using extract-mapped on all SEED
        alignments.
    """

    info = Path(info)
    info.mkdir(parents=True, exist_ok=True)

    sequences = Path(sequences)
    sequences.mkdir(parents=True, exist_ok=True)

    disallowed = disallow.Disallowed.empty()
    if disallow_file:
        logger.trace("Loading disallow file {}", disallow_file)
        with open(disallow_file, "r") as raw:
            disallowed = cattrs.structure(
                yaml.load(raw, yaml.Loader), disallow.Disallowed
            )

    all_alignments = alignments.load_jsonl(alignments_file)
    all_matches = matches.load_all(matches_file)
    families = build_families(all_alignments, all_matches)
    report = Report.empty()
    cache_dir: Path = ctx.obj["cache_path"]
    with Cache(str(cache_dir)) as cache:
        fetcher = StructureLookup.build(cache)
        planner = Planner.build(fetcher, disallowed)

        for family in families:
            logger.info("Computing actions for {}", family.rfam_accession)
            todo = planner.actions(family)
            report.track(todo)
            match todo:
                case SkippedFamily():
                    logger.info("Skipping {}", family.rfam_accession)
                case CompleteFamily():
                    logger.info("Family {} is up-to-date", family.rfam_accession)
                case NoAcceptedMatches():
                    logger.info("No accepted matches for {}", family.rfam_accession)
                case FamilyUpdate():
                    logger.info("Will update {}", family.rfam_accession)
                    seq_out = sequences / f"{todo.family.rfam_accession}.fa"
                    if not todo.has_sequences():
                        logger.info(
                            "No sequences to align for {}", todo.family.rfam_accession
                        )
                        seq_out.touch()
                    info_out = info / f"{todo.family.rfam_accession}.json"
                    with seq_out.open("w") as seq, info_out.open("w") as out:
                        todo.write_data(seq, out)
                case _:
                    assert_never(todo)
    report.write(report_file)


@main.command("parse-alignment")
@click.argument("accession")
@click.argument("alignment", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def extract_mapped_cmd(accession: str, alignment: ty.TextIO, output: ty.TextIO):
    """Extract all mapped PDB structures and their associated RNAcentral
    sequences, if any from the given alignment.

    This assumes that the given file is a single alignment.

    Arguments\b
    ---------\b
    accession:\b
      The Rfam accession for the given alignment.\b
    desc:\b
      The DESC file for the family.\b
    alignment:\b
      A Rfam SEED alignment in stockholm format.
    """

    mapped = alignments.parse_file(accession, alignment)
    raw = cattrs.unstructure(mapped)
    json.dump(raw, output)
    output.write("\n")


@main.command("add-structure-information")
@click.argument("alignment", type=click.Path())
@click.argument("pdb-info", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def add_structure_cmd(alignment: str, pdb_info: ty.TextIO, output: ty.TextIO):
    """Update the alignment with the secondary structures from the pdb-info
    file.

    This assumes that the sequences for each structure have already been
    aligned.

    Arguments\b
    ---------\b
    alignment:\b
      An alignment in stockholm format
    pdb-info:\b
      The information on the PDB structures to align.
    output:\b
      The filename to write to.
    """

    align = AlignIO.read(alignment, "stockholm")
    structure_actions = cattrs.structure(json.load(pdb_info), StructureActions)
    annotations = structure_actions.secondary_structure_annotations(align)
    with open(alignment, "r") as raw:
        for line in raw:
            output.write(line)
            maybe_id = line.split(" ", 1)[0]
            for secondary in annotations.get(maybe_id, []):
                ss_id = secondary.info.chain_id.secondary_structure_id()
                output.write(f"#=GR {maybe_id} {ss_id} {secondary.line}\n")


@main.command("clear-cache")
@click.argument("pattern", default="")
@click.pass_context
def clear_cache_cmd(ctx, pattern: str):
    pdb_cache: Path = ctx.obj["cache_path"]
    with Cache(str(pdb_cache)) as cache:
        assert isinstance(cache, (ty.Sized))  # For type hinting
        logger.info("Current cache size: {}", len(cache))
        if not pattern:
            logger.info("Removing all entries")
            cache.clear()
        else:
            logger.info("Removing all entries that match {}", pattern)
            pat = re.compile(pattern)
            to_remove = []
            for key in cache.iterkeys():
                if re.search(pat, str(key), re.IGNORECASE):
                    to_remove.append(key)
            logger.info("Will remove {} keys", len(to_remove))
            for key in to_remove:
                logger.debug("Removing {}", key)
                cache.delete(key)
        logger.info("Final cache size: {}", len(cache))


if __name__ == "__main__":
    main()
