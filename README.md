# Rfam 3D Seed Alignments

The goal of this project is to automate the incorporation of the 3D structural information into the [Rfam](https://rfam.org) seed alignments using the following workflow:

- The Rfam-PDB mapping file is used to find out which PDB files need to be added to which seed alignments.

- The 3D structural annotations are downloaded from the [RNA 3D Hub](http://rna.bgsu.edu/rna3dhub) database which regularly annotates all RNA 3D structures using [FR3D](http://rna.bgsu.edu/FR3D).

- The PDB sequences and secondary structures in dot-bracket notation are iteratively added to the Rfam seed alignments using the _cmalign_ [Infernal](http://eddylab.org/Infernal) program.

- The PDB accessions are replaced with [RNAcentral](https://rnacentral.org) identifiers in the final alignments.

## Installation

- Download the repository or use `git clone`.

- Start an interactive session using Docker:

    ```
    docker run rfam
    ```

## Usage

To update a specific Rfam family:

```
add-3d.py RF00162
```

Alternatively, run `add-3d.py` to update all families.

To get FR3D secondary structure for a PDB id:

```
fr3d_2d.py 2QUS_B
>2QUS_B
GGGAGCCCUGUCACCGGAUGUGCUUUCCGGUCUGAUGAGUCCGUGAGGACAAAACAGGGCUCCCGAAUU
.((((((((((.((((((.....{.))))))(....).((((...}))))...))))))))))......
```

The updated seed alignments with the added 3D structures will be in the `output` folder (see [precomputed results](./data/output)).

## Updating the data

- Delete `pdb_full_region.txt` to download the latest Rfam-PDB mapping from the Rfam FTP archive.

- Delete `pdb.tsv` to download the latest RNAcentral-PDB mapping from the RNAcentral FTP archive.

## Feedback

Please feel free to raise an issue to report any problems with the code or the data.

## Acknowledgements

We would like to thank [Sri Devan Appasamy](http://sridevan.me) and [Craig Zirbel](https://www.bgsu.edu/arts-and-sciences/mathematics-and-statistics/faculty-and-staff/craig-zirbel.html) for developing an RNA 3D Hub API to provide FR3D annotations for RNA 3D structures.
