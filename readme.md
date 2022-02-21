# lvar_anno

Add rRNA and mitochondrial annotations to the Lvar 3.0 genome.

## Running the pipeline.

1. Ensure snakemake is installed, and that `bedtools`, `bedops`, and `gffutils`.
Python package are accessible from path. If [agat](https://github.com/NBISweden/AGAT#installation) is available,
the final two rules can be run to produce error-checked `.gff` and `.gtf` annotation files.

2. Ensure files names in `files/config.yaml` point to the correct locations.

3. Run the pipeline issueing the command `snakemake` from the parent directory.

## References
1. Lytechinus variegatus genome provided by PL Davidson et al. from the Wray lab at Duke
    - Assembly: https://db.cngb.org/search/assembly/CNA0013733/
    - Paper: https://academic.oup.com/gbe/article/12/7/1080/5841217

2. Mitochondrial genome retrieved from Bronstein and Krol.
    - Data:  https://www.ncbi.nlm.nih.gov/nuccore/MG676468
    - Paper: https://www.sciencedirect.com/science/article/pii/S088875431730174X#s0070

3. rRNA annotations provided by PL Davidson 

