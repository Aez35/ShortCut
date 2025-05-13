# ShortCut
small RNA trimming, library quality control, and alignment all in one go! Currently in development.


### Create environment and install dependencies
```
conda create -n shortcut shortstack r-base r-tidyverse r-ggplot2 cutadapt
conda install -c conda-forge rpy2
```

# Usage

## ShortCut
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



# Results
## Cutadapt summary
A summary of libraries that were trimmed using Cutadapt. This provides the user an idea of how many total reads were in their unprocessed libraries, and how successful adapter trimming was at trimming these reads.

![Cutadapt_summary](https://github.com/user-attachments/assets/61996985-6490-4810-8daf-146b8dcffaef)

## Library Summary 
Barplot of all library sizes (including libraries that did not require trimming)
![Library_summary](https://github.com/user-attachments/assets/031fa326-e5be-4465-a818-88236d138316)

## Read size distribution PDF
PDF containing the sRNA read size distribution for each library. sRNA libraries generally expect a peak at 21- and 24-nt in accordance with miRNA and siRNA peaks.
<img width="1439" alt="image" src="https://github.com/user-attachments/assets/5882e7e0-ae73-42f0-8bc5-8bd153b9c4d3" />


