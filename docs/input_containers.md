# Containers

Several containers are used within the wgs pipeline. There is a "default" container used for many processes, and some more narrow containers used for specific pipelines. For processes which use other than the default, containers are specified in the processes.

To build a Singularity container, you run the following command.

```
sudo singularity build name_of_your_container.sif Singularity
```

To retain the environment variables from the user, add the `-E` flag to the `sudo` command.

```
sudo -E singularity build name_of_your_container.sif Singularity
```

## The base container

The base container can be built by navigating into the `container` directory and running the following command:

```
sudo -E singularity build name_of_your_container.sif Singularity
```

The container can be directly downloaded from [here](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/master/container/Singularity).

## Madeline2

`madeline2` is installed separately. The Singularity recipe for building it is found [here](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/master/container/Singularity_madeline2).

## VEP

The VEP container is built from a publicly available docker container. Currently, the pipeline uses version 103.

```
singularity pull docker://ensemblorg/ensembl-vep:release_103
```

When running this container, additional reference files needs to be provided. These are outlined in the [[4.2 Inputs - Annotation files]] page.

## CADD

The CADD container can be retrieved from the [Galaxy Project](https://depot.galaxyproject.org/singularity/). 


Alternatively (if working at CMD in Lund), it can be clone from GitHub as such:

```
git clone https://github.com/Clinical-Genomics-Lund/CADD-container/tree/v1.6
```

In addition, reference files needs to be downloaded. Note that the file is really large - 200GB zipped. Use `wget` with the `-c` flag such to allow continuing the download if it is interrupted.

```
wget -c https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/annotationsGRCh38_v1.6.tar.gz
```

Make sure that the lines in the `Singularity` recipe matches to where you have saved the annotations folder, i.e. the path in the following line:

```
ln -sr /fs1/resources/ref/hg38/annotation_dbs/CADD_v1.6/data/annotations/GRCh38_v1.6 GRCh38_v1.6
```

Note: At the moment this container assumes that you have access to the `/fs1/viktor/misc-scripts` directory, i.e. you are running at the servers in Lund.

Finally you build the container:

```
sudo -E singularity build cadd_v1.6.sif Singularity
```

## GATK

This container is used for CNV-calling. It can retrieved from a Docker container provided by Broad Institute.

```
singularity pull docker://broadinstitute/gatk
```

## Sentieon

Normally, Sentieon containers are retrieved from the [Galaxy Project](https://depot.galaxyproject.org/singularity/)

The container for the 2023 version seems not to work correctly (i.e. it crashes when running `sentieon driver`). An alternative container was be built from the Docker container found [here](https://github.com/Sentieon/sentieon-docker) for the 2023 version.

## Stranger

A Singularity recipe for Stranger is found [here](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/blob/master/container/stranger/Singularity). 

Alternatively, the Stranger container can be retrieved from the [Galaxy Project](https://depot.galaxyproject.org/singularity/).

## Twist myeloid

Different versions of the twist-myeloid container are used for the `cnv-kit` and `freebayes` processes.

**FIXME**: While we still use this container, would we include the recipe?

## AnnotSV

A repository containing the recipe to build the AnnotSV container can be retrieved as such:

```
git clone https://github.com/Clinical-Genomics-Lund/annotsv_container.git
```

To build the container, `cd` into the repository and run the following command:

```
sudo singularity build annotsv.v2.3.sif Singularity
```

## POD

A POD container can be built using the recipe found [here](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/tree/master/container/pod).

## Genmod

A Genmod container can be retrieved from the [Galaxy Project](https://depot.galaxyproject.org/singularity/).

## REViewer

A REViewer container can be built using the recipe found [here](https://github.com/Clinical-Genomics-Lund/nextflow_wgs/tree/master/container/reviewer).

## Melt

A Melt container can be retrieved from the [Galaxy Project](https://depot.galaxyproject.org/singularity/).
