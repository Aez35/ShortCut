# ShortCut
small RNA trimming, library quality control, and alignment all in one go! Currently in development.


### Create environment and install dependencies
```
conda create -n shortcut shortstack r-base r-tidyverse r-ggplot2
conda install -c conda-forge rpy2
```

# Usage

## Plotrim
```
python ShortCut.py [-help] -fastq FASTQFILES.fq/fastq -kingdom plant/animal [-trimkey KEY] [-threads THREADS] [-out outputDirectory, default= shortcut_output/] [-annotate] [-genome GENOME.fa] [-known_mirnas MIRNAS.fa] [-dn_mirna] [-ssout ssOutputDirectory]
```

## Options

|Option     |Description                                                     |
|:---------:|:-----------------------------------------------------------:   |
|help       | prints help message                                            |
|fastq      | List of small RNA-seq libraries in FASTQ format (zipped ok)    |
|kingdom    | Specify either 'plant' or 'animal'                             |
|out        | Output directory, default is plotrim_output                    |
|threads    | Specify number of threads to use during trimming step          |
|annotate   | Switch to trigger sRNA annotations with ShortStack             |
|known_mirnas| FASTA file containing known miRNAs for ShortStack annotating  |
|genome     | Genome used for ShortStack alignment and annotation of sRNAs   |
|trimkey    | miRNA sequence for detecting adapter sequences for trimming    |
|ssout      | ShortStack output directory. ShortStack will create default    |
|dn_mirna   | De novo miRNA search option in ShortStack                      |

