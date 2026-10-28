# Getting started

To run the pipeline you will need:

* The source code
* A `nextflow.config` file
* An installation of [Nextflow](https://www.nextflow.io)
* An installation of [Singularity](https://docs.sylabs.io/guides/latest/user-guide)
* [Annotation files](annotation_files.md)

## Source code

Clone the source code from [GitHub](https://github.com/Clinical-Genomics-Lund/nextflow_wgs).

```
git clone https://github.com/Clinical-Genomics-Lund/nextflow_wgs.git
```

## The `nextflow.config` file

This file contains various settings and annotation files required to run the pipeline. Template config files (used if running on Hopper or Trannel in Lund) are found in the `configs` folder. Copy the file next to your `main.nf` file and adjust it to your needs. More information on the annotation files are found in the [annotation files section](annotation_files.md).

## Dependencies

### Nextflow

[Nextflow](https://www.nextflow.io/) is a programming language designed to build workflows. At the moment, the pipeline is implemented using Nextflow's DSL1 syntax. More recent versions of Nextflow only supports the new DSL2 syntax. Thus an older version of Nextflow (latest 21) is required to run the workflow.

As mentioned on the [Nextflow page](https://www.nextflow.io/docs/latest/getstarted.html), it is possible to run an older version as such:

```
NXF_VER=20.04.0 nextflow run hello
```

To run Nextflow, you will also need to have Java installed.

### Singularity

[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) is used to manage dependencies for individual steps of the pipeline, called "processes". This means that you will not need to install the required dependencies on the computer that you will run on - the processes can execute directly inside these containers. These containers are further discussed in the [container section](input_containers.md) of this README.

### Loading dependencies on a cluster

If you are running this on a computational cluster providing the `module` command, the required dependencies can be loaded as such:

```
module load Java
module load nextflow
module load singularity
```

## Running the pipeline

### Stub-runs

Stub-runs allow you to test-run the pipeline without actually performing the full analysis. The stub rub creates dummy files for each processing steps, and completes in a matter of minutes. This is a useful tool for testing and debugging the pipeline.

```
nextflow run main.nf -stub-run
```

### Running real data

To run a full dataset, you need to provide a profile (`-profile` argument) and input CSV (`--csv` argument):

```
nextflow run main.nf \
    -profile wgs \
    --csv path/to/input.csv
```

Note the difference between single dash arguments (`-profile`) and double dash arguments (`--csv`). Single-dashed arguments are provided directly to Nextflow while double-dashed arguments are provided as params to the workflow itself. These overrides any params specified in the `nextflow.config` file.

Additional useful arguments:

* `-resume` If possible, continue a previous run and only rerun required steps.
* `-w` Specify the location of the work folder (default is in the same folder as `main.nf`)

### Refresh SV results after the AnnotSV join repair

The TSV join now verifies original variant IDs, alleles and positions instead
of requiring AnnotSV's derived deletion endpoint to equal the VCF REF span.
It also supports exact allele matches after IDs change, skips split-mode rows,
and rejects conflicting annotations. Missing classes remain missing. This does
not change the scoring model, variant alleles, genotypes or calling settings.

Deploy `main.nf`, `bin/add_annotsv.py` and `bin/annotsv_join.py` together. The join
module is an explicit task input so its changes invalidate the cached join.
No container rebuild is required: it uses the existing Python/pysam/cmdvcf
dependencies and the scripts from the pipeline checkout.

Resume the original production command with its original configuration, inputs,
launch directory, work directory and `-resume` (or the original run/session ID).
Nextflow can reuse alignment, callers, VEP, database queries and AnnotSV TSVs.
The required refresh starts at `add_annotsv_to_svvcf`, then runs its downstream
OMIM, caller-penalty, genetic-model and `score_sv` stages, compression/indexing,
and dependent final exports. Existing AnnotSV TSVs already contain the classes;
they do not need regenerating solely for this join repair.

The joined VCF contains `##OLWGS_ANNOTSV_JOIN=original_allele_v1`. Keep baseline
outputs separately if a published output directory would otherwise overwrite
them. Update downstream comparison/report inputs to the refreshed scored VCFs;
already imported cases require the site's normal refresh/reimport to show new
scores. Re-running reports against old scored VCFs does not apply the fix.

Validate the join helper with:

```bash
python -m pytest -q tests/test_annotsv_join.py
```

### Running as a batch script

Note that if you are running jobs on Hopper in CMD, Lund, production jobs should be started using the `start_nextflow_analysis.pl` script.

Otherwise, if working on a computational cluster, the jobs will typically be executed using SLURM. A minimal example of a SLURM run script is shown below. 

```
# SBATCH --job-name=job_name
# SBATCH --output=slurm_log_%j.log
# SBATCH --ntasks=2
# SBATCH --mem=4gb
# SBATCH --time=2-00:00:00

module load Java
module load nextflow
module load singularity

nextflow run main.nf \
    -profile wgs \
    --csv path/to/input.csv
```

Assuming it is named `jobfile.run`, then it can be queued by running:

```
sbatch jobfile.run
```

