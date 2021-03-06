# Rfam 3D Seed Alignments

The goal of this project is to automate the incorporation of the 3D structural information into the [Rfam](https://rfam.org) seed alignments using the following workflow:

- The Rfam-PDB [mapping file](./pdb_full_region.txt) is used to find out which PDB files need to be added to which seed alignments.

- The 3D structural annotations are downloaded from the [RNA 3D Hub](http://rna.bgsu.edu/rna3dhub) database which regularly annotates all RNA 3D structures using [FR3D](http://rna.bgsu.edu/FR3D).

- The PDB sequences and secondary structures in dot-bracket notation are iteratively added to the Rfam seed alignments using the _cmalign_ [Infernal](http://eddylab.org/Infernal) program.

- The PDB accessions are replaced with [RNAcentral](https://rnacentral.org) identifiers in the final alignments.

## Installation

- Download the repository or use `git clone`.

- Start an interactive session using Docker:

    ```
    docker-compose run rfam
    ```

Alternatively, follow instructions in the [Dockerfile](./Dockerfile) to install locally.

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

## Feedback

Please feel free to [raise an issue](https://github.com/Rfam/rfam-3d-seed-alignments/issues) to report any problems with the code or the data.

## Acknowledgements

We would like to thank [Sri Devan Appasamy](http://sridevan.me) and [Craig Zirbel](https://www.bgsu.edu/arts-and-sciences/mathematics-and-statistics/faculty-and-staff/craig-zirbel.html) for developing an RNA 3D Hub API to provide FR3D annotations for RNA 3D structures.
