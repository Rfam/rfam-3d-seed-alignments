# Rfam 3D Seed Alignments

The goal of this project is to automate the incorporation of the 3D structural information into the [Rfam](https://rfam.org) seed alignments using the following workflow:

- The Rfam-PDB mapping file is used to find out which PDB files need to be added to which seed alignments.

- The 3D structural annotations are extracted from the [RNA 3D Hub](http://rna.bgsu.edu) database which regularly annotates all RNA 3D structures using [FR3D](http://rna.bgsu.edu/FR3D).

- The PDB sequences and secondary structures in dot-bracket notation are iteratively added to the Rfam seed alignments using the _cmalign_ [Infernal](http://eddylab.org/Infernal) program.

- The PDB accessions are replaced with RNAcentral identifiers in the final alignments.

## Installation

- Download the repository or use `git clone`.

The following steps are only required to update 3D annotations from RNA 3D Hub:

```
virtualenv env
pip install -r requirements.txt
source env/bin/activate
```

## Usage

```
python add-3d.py
```

The updated seed alignments with the added 3D structures will be in the `output` folder (see [precomputed results](./output)).

## Updating the data

- Delete `pdb_full_region.txt` to download the latest Rfam-PDB mapping from the Rfam FTP archive.

- Delete `pdb.tsv` to download the latest RNAcentral-PDB mapping from the RNAcentral FTP archive.

- Run `python export-basepairs.py` to update 3D annotations (last updated on June 16, 2020) - _requires RNA 3D Hub credentials_.

## Feedback

Please feel free to raise issues to report any problems with the code or the data.
