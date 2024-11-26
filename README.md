# Rfam 3D Seed Alignments

The goal of this project is to automate the incorporation of the 3D structural information into the [Rfam](https://rfam.org) seed alignments using the following workflow:

- The Rfam-PDB [mapping file](./pdb_full_region.txt) is used to find out which PDB files need to be added to which seed alignments.

- The 3D structural annotations are downloaded from the [RNA 3D Hub](http://rna.bgsu.edu/rna3dhub) database which regularly annotates all RNA 3D structures using [FR3D](http://rna.bgsu.edu/FR3D).

- The PDB sequences and secondary structures in dot-bracket notation are iteratively added to the Rfam seed alignments using the _cmalign_ [Infernal](http://eddylab.org/Infernal) program.

- The PDB accessions are replaced with [RNAcentral](https://rnacentral.org) identifiers in the final alignments, if needed.

## Usage

- To update one or more Rfam families:

  ```
  add_3d.py RF00162
  add_3d.py RF00162 RF00507
  ```

  Use `--nocache` to force recomputing the output and download the latest PDB-Rfam and PDB-RNAcentral mapping files.

- To update all families:

  ```
  add_3d.py all --nocache
  ```

- To get FR3D secondary structure for a PDB id:

  ```
  fr3d_2d.py 2QUS_B
  >2QUS_B
  GGGAGCCCUGUCACCGGAUGUGCUUUCCGGUCUGAUGAGUCCGUGAGGACAAAACAGGGCUCCCGAAUU
  .((((((((((.((((((.....{.))))))(....).((((...}))))...))))))))))......
  ```

The updated seed alignments with the added 3D structures will be in the `output` folder (see [precomputed results](./data/output)).

## Manually curated Rfam-PDB mapping file

It is possible to manually add mapping between Rfam accessions and PDB ids to [pdb_full_region_curated.txt](./pdb_full_region_curated.txt). This step is needed in order to analyse PDB sequences that do not match Rfam covariance models automatically. This can happen when a PDB sequence gets a bit score below the Rfam threshold because it is much shorter than the corresponding Rfam model.

## Workflow

Some things to clarify:

- This is meant to update a SEED alignment of a given family. This does not
  edit the CMs, just suggested a new SEED.

- A PDB structure is made of 1 or more chains. The chains in the structure may
  be of the same or different molecules. The 3D structure of a molecule does not
  have to be consitent between different experiments.

The chain observed in a structure, here called the 'chain' or the 'observed
sequence', may be incomplete relative to the sequence that was used in the
experiment, here called the 'experimental sequence'. This is because
experiments may not have enough electron density to see all nucleotides.

An experimental sequence may match one or more Rfam families. In this case
matching means the sequence has a hit above a families bit score threshold.

Each match within an experimental may be complete relative to the model or not.

Basepairs from 3D structures are annotated with FR3D. We do not run FR3D
ourselves, but instead fetch the pairing information from RNA BGSU.

We are interested in aligning all matches to experimental sequences from all 3D
structures into Rfam SEED alignments.

The aligned matches are used to infer the alignment of basepairs from 3D
structures. The basepairing information may or may not be used to udpate the
SEED alignment with new basepairs. Some factors that will effect the addition
of basepairs:

- Resolution. Likely no reason to include basepairs from models 3.5Å structures

- Unusual structure. Not all solved structures are a native conformation, nor
  are they in the state that Rfam models.

In general this pipeline fetches all novel experimental sequences, aligns them
into the approiate models and then infers the new base pair annotations.

Basepair annotations are per letter (GR lines) annotations to their
experimental sequence.

The sequence is meant to be an RNAcentral identifier if possible. If it is not
possible a hash is used to force all structures with the same sequence into the
same set of GR annotations instead of spread across several sequences.

While the match to an experimental sequence may not include the entire
sequence, this will align the entire sequence into the alignment. This is
because it is hard to know if the alignment is too short or not ahead of time.
Sometimes they are, sometimes aligning the whole thing in is a mistake.

TODO In order to allow for maximum flexibility there is a script to align in
only the matching region.

This fetches the current SEED and CMs from Rfam's public SVN repository ([https://svn.rfam.org/]).

Each SEED alignment is parsed to find which, if any, structures and associated sequences, are already a member of the seed alignment.

The mapping file is then read to determine which new pdbs

## Feedback

Please feel free to [raise an issue](https://github.com/Rfam/rfam-3d-seed-alignments/issues) to report any problems with the code or the data.

## Acknowledgements

We would like to thank [Sri Devan Appasamy](http://sridevan.me) and [Craig Zirbel](https://www.bgsu.edu/arts-and-sciences/mathematics-and-statistics/faculty-and-staff/craig-zirbel.html) for developing an RNA 3D Hub API to provide FR3D annotations for RNA 3D structures.
