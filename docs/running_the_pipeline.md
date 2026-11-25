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

### Annotating an existing VCF

For a single-sample SNV/indel VCF, add the `snv_vcf` profile after the WGS
profiles. Put the indexed, normalized VCF and its index in the samplesheet's
`read1` and `read2` columns. The VCF sample name must match `id` in the
samplesheet. This follows the existing SNV annotation and scoring path.

```bash
nextflow -c configs/nextflow.hopper.config run main.nf \
    -profile wgs,hg38_masked,snv_vcf \
    --csv /path/to/sample.csv \
    --outdir /path/to/evaluation/snv \
    --crondir /path/to/evaluation/snv/cron
```

For a sorted, unannotated GRCh38 SV VCF, run the standalone SV entry point
with the Hopper VEP and gnomAD resources. It runs VEP and gnomAD annotation and
writes an indexed VCF to `<outdir>/vcf/<sample_id>.gnomad.vcf.gz`. It ends before
the production SVDB merge and scoring steps.

```bash
nextflow -c configs/nextflow.hopper.config run annotate_sv.nf \
    --vcf /path/to/unannotated.sv.vcf.gz \
    --sample_id HG002 \
    --outdir /path/to/evaluation/sv
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

