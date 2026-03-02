# ShortCut
Adapter trimming and length distributions from small RNA-seq data

# Introduction
Small RNA sequencing usually produces reads that are longer than the biological RNAs. The extra sequence derives from adapters and is present at the 3' end of each read. Before analysis, these residual adapter sequences need to be identified and removed. After trimming, the size distribution of the trimmed reads is informative to assess sample quality. High rates of RNA degradation in the input are often reflected in the size profile of the trimmed reads. 

ShortCut automates adapter discovery, trimming, and provides an easy to use report on the results. The trimmed reads are output in a condensed FASTA format (see below) and are ready for further analysis.

## Input
One or more FASTQ files (can be .gz compressed) containing untrimmed small RNA seq data.

## Output
- Three condensed FASTA files per input FASTQ file. The three files contained properly trimmed sequences, sequences that were too short to be acceptable after trimming, and sequences where no adapter was identified, respectively.
- A csv-formatted file containing summary data on the output.
- An html report file giving an overview of the results.

# Install
1. Install `pixi` on your system <https://pixi.prefix.dev/latest/installation/>
2. Configure `pixi` to use BioConda per <https://bioconda.github.io/>
3. ```
   mkdir myShortCut
   cd myShortCut
   pixi init
   pixi add pandas cutadapt plotly
   pixi shell
   ```
4. Download the `ShortCut` file, `chmod +x ShortCut`, and run using `./ShortCut`. Or install to your path as you see fit.

(Mike to do - place this version of ShortCut on BioConda to simplify installation)

# Usage
With `pixi shell` activated as described above:
```
usage: ShortCut [-h] --fastq FASTQ [FASTQ ...] [--min_read_size MIN_READ_SIZE] [--out_dir OUT_DIR] [--threads THREADS] [--dicermin DICERMIN] [--dicermax DICERMAX]
                [--trim_key TRIM_KEY | --adapter ADAPTER]

options:
  -h, --help            show this help message and exit
  --fastq FASTQ [FASTQ ...]
                        One or more FASTQ alignment files
  --min_read_size MIN_READ_SIZE
                        Minimum length of trimmed reads to retain
  --out_dir OUT_DIR     output directory
  --threads THREADS     Number of threads to use. Default=1
  --dicermin DICERMIN   Minimum trimmed read length for true regulatory small RNA. Default=21
  --dicermax DICERMAX   Minimum trimmed read length for true regulatory small RNA. Default=24
  --trim_key TRIM_KEY   Trim key to use to infer adapter sequence (default: UCGGACCAGGCUUCAUUCCCC)
  --adapter ADAPTER     Adapter sequence to use for trimming
```

# Test / Example
Get some example small RNA-seq files from SRA. (Assumes installation of `fasterq-dump` from `sra-tools`).
`fasterq-dump SRR3222443 SRR3222444 SRR3222445`

Then, using four threads:
`ShortCut --fastq *.fastq --threads 4`

## Output details
### Condensed FASTA files
Each input FASTQ file results in three condensed FASTA files: 
- The 'Short' directory contains FASTA files beginning with 's_'. These result from reads that were successfully adapter trimmed, but after trimming were shorter than the length set by option `--min_read_size` (default = 12).
- The 'trimmed_ok' directory contains FASTA files beginning with 't_'. These result from reads that were successfully adapter trimmed, and were >= `--min_read_size` after trimming.
- The 'untrimmed' directory contains FASTA files beginning with 'u_'. These result from reads were the adapter was not identified, and thus were not trimmed.

Small RNA-seq data often contain multiple reads with the same sequence. These are often microRNAs or siRNAs that were highly abundant in the sample. "Read condensation" refers to the process of condensing the input for alignment such that each unique sequence is only represented once in the condensed FASTA file. Suppose the sequence TCGGACCAGGCTTCATTCCCC was the sequence for 509,066 reads. Instead of repetitively writing, and aligning, the same sequence over and over we can instead write a single FASTA entry. Example of condensed format:
```
>t_SRR3222443_Cd1_509066
TCGGACCAGGCTTCATTCCCC
```

- `t_SRR3222443` : Prefix specific for the input read file. The `t_` indicates trimmed.
- `Cd1` : `Cd` means "condensed". The integer after is a unique integer within that file.
- `_509066` : The trailing integer that represents the total number of reads that have this sequence. In this case, 509,066 reads with this sequence.

### Reads vs. Sequences
A **read** is a single FASTQ output from the sequencer. A **sequence** is a unique sequence of DNA/RNA which may have been observed one or more times in an experiment.

### ShortCut_tallys.csv
This is a .csv file containing summary information on read counts and sequence counts by RNA length, Category (either trimmed_ok, short, or untrimmed), and Source file.  [View example csv](example/ShortCut_tallys.csv)


### ShortCut_Results.html
This is a report that is viewable with a web browser. It shows several interactive plots that are informative for assessing library trimming results. [View example report](example/ShortCut_Results.html)





