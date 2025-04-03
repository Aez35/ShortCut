#!/usr/bin/env python
version = '0.1'

import re, gzip
import argparse, shutil, subprocess
from subprocess import Popen, PIPE, STDOUT
import os, glob
from pathlib import Path
import csv, sys, time
from csv import writer
import rpy2.robjects as robjects
robjects.r('library(ggplot2)')
robjects.r('library(tidyverse)')


print("Checking required modules...")
req_packs=["bowtie","samtools","RNAfold","cutadapt"]
for pack in req_packs:
        path = shutil.which(pack)
        if path is not None:
            print('Requires package {0} : {1}'.format(pack,path))
        else:
            msg = 'Requires package {0} : Not found!'.format(pack)
            sys.exit(msg)


ap = argparse.ArgumentParser()

ap.add_argument('-fastq',nargs='*',required=True,help='One or more fastq alignment files')
ap.add_argument("-n", help="Results file name")
ap.add_argument("-threads", help="Specify number of threads for samtools", default = 1)
ap.add_argument("-kingdom", required=True,choices=["plant","animal"],help="Specify animal or plant")
ap.add_argument("-annotate",help="Specify whether you want ShortStack to annotate",nargs='?', const='')
ap.add_argument("-out", help="output directory",default="plotrim_output")
ap.add_argument('-genome',help='Genome for annotation')
ap.add_argument('-known_mirnas',help='FASTA-formatted file of known mature miRNA sequences')
ap.add_argument("-ssout", help="output directory",default="Shortstack_output")
ap.add_argument("-trimkey", help="Abundant miRNA used to find adapters for trimming with option -autotrim")
ap.add_argument("-dn_mirna", help="De novo miRNA search in ShortStack")
args = ap.parse_args()

#print log
print("--------")
print("sRNATrim version " + version)
print("Options:")
print("     'Threads' " + str(args.threads))
print("     'Kingdom' "+ args.kingdom)
print("     'Output directory' "+ args.out+"/")
if args.fastq is not None:
    print("     'Fastqs' " +str(args.fastq))
if args.annotate is not None:
    print("     'Annotate sRNAs' yes")
    print("     'Genome' " +str(args.genome))
    print("     'known_mirnas' " +str(args.known_mirnas))
    if args.dn_mirna == True:
        print("     'denovo miRNAs' " +str(args.known_mirnas))
    print("     'ssout' " +str(args.ssout))

print("--------")

if args.annotate is True:
    req_packs=["ShortStack"]
    for pack in req_packs:
        path = shutil.which(pack)
        if path is not None:
            print('Requires package {0} : {1}'.format(pack,path))
        else:
            msg = 'Requires package {0} : Not found!'.format(pack)
            sys.exit(msg)

path=args.out
isExist = os.path.exists(path)
if isExist==True:
    msg = ("Output directory '" + args.out + "/' already exists. Please assign a new output directory using option '-out'.")
    sys.exit(msg)

#_____________ Functions ______________
def run(cmd) :
#Run subprocess
    proc = subprocess.call(cmd, shell=True,stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)

def process_fastqs(fastqs, keydna, args):
    for fastq in fastqs:
        print("Detecting adapter for " + fastq + " with key "  + keydna)
        if not os.path.exists(fastq):
            sys.exit(f"Error: File {fastq} does not exist. Please check the path.")

        # Adapter detection command
        if fastq.endswith(".gz"):
            cmd = (
                f"gzip -dc {fastq} | awk -v target='{keydna}' '{{idx = index($0, target); if (idx) print substr($0, idx + length(target),20)}}' | "
                "uniq -c | sort -nr | head -n1 | awk '{print $2}'"
            )
        else:
            cmd = (
                f"awk -v target='{keydna}' '{{idx = index($0, target); if (idx) print substr($0, idx + length(target),20)}}' {fastq} | "
                "uniq -c | sort -nr | head -n1 | awk '{print $2}'"
            )

        try:
            output = subprocess.check_output(cmd, shell=True, text=True).strip()
        except subprocess.CalledProcessError as e:
            sys.exit(f"Error while processing {fastq}: {e.output}")

        if "G" not in output:
            print("No adapter detected. Copying untrimmed file to trimmedLibraries instead.")
            
            # Create destination filename
            head, tail = os.path.split(fastq)
            tfile = os.path.join('trimmedLibraries', f"t_{tail}")
            
            # Copy file as-is
            shutil.copy(fastq, tfile)
            
            continue  # move on to the next file
        else:
            print(f"Adapter detected: {output}")

        head, tail = os.path.split(fastq)
        tfile = os.path.join('trimmedLibraries', f"t_{tail}")
        report_file = os.path.join('trimmedLibraries', f"{tail}_cutadapt_report.txt")

        # Construct the Cutadapt command
        if fastq.endswith(".gz"):
            cmd = (
                f"gzip -dc {fastq} | cutadapt -j {args.threads} -a {output} -o {tfile} -m 12 -"
            )
        else:
            cmd = f"cutadapt -j {args.threads} -a {output} -o {tfile} -m 12 {fastq}"

        # Run Cutadapt and capture its output
        try:
            result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
            with open(report_file, "w") as rf:
                rf.write(result.stdout)
        except Exception as e:
            sys.exit(f"Error trimming {fastq}: {e}")

def summarize_cutadapt_reports(trimmed_dir='trimmedLibraries', output_csv='trimming_summary.csv'):
    summary = {}

    for filename in os.listdir(trimmed_dir):
        if filename.endswith("_cutadapt_report.txt"):
            sample = filename.replace("_cutadapt_report.txt", "")
            filepath = os.path.join(trimmed_dir, filename)
            
            with open(filepath) as f:
                content = f.read()

            total_reads=re.search(r"Total reads processed:\s+([\d,]+)", content)
            if total_reads:
                total_count = total_reads.group(1).replace(",", "")


            # Extract "Reads with adapters"
            match_adapters = re.search(r"Reads with adapters:\s+([\d,]+) \(([\d.]+)%\)", content)
            if match_adapters:
                adapters_count = match_adapters.group(1).replace(",", "")
                adapters_percent = match_adapters.group(2) + "%"
            else:
                adapters_count = adapters_percent = "N/A"

            # Extract "Reads written"
            match_written = re.search(r"Reads written \(passing filters\):\s+([\d,]+) \(([\d.]+)%\)", content)
            if match_written:
                written_count = match_written.group(1).replace(",", "")
                written_percent = match_written.group(2) + "%"
            else:
                written_count = written_percent = "N/A"

            summary[sample] = {
                "Reads with adapters": adapters_count,
                "% Reads with adapters": adapters_percent,
                "Reads written (passing filters)": written_count,
                "% Reads written": written_percent,
                "Total reads": total_count
            }

    # Transpose into a dataframe-like CSV structure
    metrics = ["Reads with adapters", "% Reads with adapters", 
               "Reads written (passing filters)", "% Reads written", "Total reads"]

    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        header = ["Metric"] + list(summary.keys())
        writer.writerow(header)
        
        for metric in metrics:
            row = [metric]
            for sample in summary:
                row.append(summary[sample].get(metric, "N/A"))
            writer.writerow(row)
    
    print(f"✅ Cutadapt summary written to {output_csv}")

def DNAcheck(dna):
    """
    Check whether a DNA/RNA sequence contains only valid nucleotides (A, T, G, C, U).

    Parameters:
        sequence (str): DNA/RNA sequence to validate.

    Returns:
        int: 2 if the sequence is valid, 1 otherwise.
    """
#Check that all characters in provided sequences are ATGCU
    y = dna.upper()
    if re.match("^[ATGCU]+$", y):
        return(2)
    else:
        return(1)

def get_distribution_from_trimmedLibraries(trimmedLibraries):
    # Get a list of all FASTQ and gzipped FASTQ files in the directory
    fastq_files = [
        os.path.join(trimmedLibraries, f)
        for f in os.listdir(trimmedLibraries)
        if f.endswith('.fastq') or f.endswith('.fastq.gz')
    ]
    
    if not fastq_files:
        print("No FASTQ files found in the trimmedLibraries.")
        return

    # Dictionary to store read length distributions
    distributions = {}

    for fastq_file in fastq_files:
        length_count = {}

        # Open normally if .fastq, else use gzip for .fastq.gz
        open_func = gzip.open if fastq_file.endswith('.gz') else open

        with open_func(fastq_file, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 1:  # Sequence lines (2nd line in every FASTQ entry)
                    read_length = len(line.strip())
                    length_count[read_length] = length_count.get(read_length, 0) + 1

        distributions[fastq_file] = length_count

    # Find the union of all read lengths across all files
    all_lengths = sorted({length for dist in distributions.values() for length in dist})

    # Prepare CSV data
    csv_data = [['Read Length'] + [os.path.basename(file) for file in fastq_files]]  # Header row

    for length in all_lengths:
        row = [length]  # First column is the read length
        for file in fastq_files:
            row.append(distributions[file].get(length, 0))
        csv_data.append(row)

    # Write to CSV file
    output_csv = 'read_length_distribution.csv'
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(csv_data)

    print(f"CSV file '{output_csv}' created successfully.")

def run_ss(t_fastq_list, genome, known_mirnas):
    print("Annotating sRNAs in " + str(args.known_mirnas))
    if args.dn_mirna == True:
        print("Conducting de novo annotation of sRNAs...")

        command=(f"ShortStack --genomefile ../{genome} --readfile {t_fastq_list}* --dn_mirna --known_miRNAs ../{known_mirnas}")
        run(command)
    else:
        command = f"ShortStack --genomefile ../{genome} --readfile {t_fastq_list}* --known_miRNAs ../{known_mirnas}"
        # Open process and stream output in real-time
    # Use Popen for real-time output streaming
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1, universal_newlines=True)

    # Print output line by line in real time
    for line in process.stdout:
        print(line, end="")  # Print and flush each line
        sys.stdout.flush()   # Ensure immediate printing

    # Print error output if any
    for line in process.stderr:
        print(line, end="", file=sys.stderr)
        sys.stderr.flush()

    process.wait()  # Ensure the process completes



#_________ run ____________
def main():
            subprocess.call(["mkdir",args.out])
            os.chdir(args.out)

            fastqs = ['../' + item for item in args.fastq]

            print("_____________________________")
            print(" ")
            print("Trimming adapters...")
            #If no key provided, set key according to kingdom.
            if args.trimkey is None:
                subprocess.call(["mkdir","trimmedLibraries"])
                if args.kingdom == 'plant':
                    #ath-miR166a
                    print("No key provided for trimming.")
                    key="UCGGACCAGGCUUCAUUCCCC"
                    print("     trimkey used for discovering adapter sequence: "+key)
                else:
                    #hsa-let-7a
                    print("No key provided for trimming.")
                    key="UGAGGUAGUAGGUUGUAUAGUU"
                    print("     trimkey used for discovering adapter sequence: " + key)
            else:
                subprocess.call(["mkdir","trimmedLibraries"])
                key=args.trimkey
            #Change key to uppercase and to DNA
            key=key.upper()
            keydna=key.translate(str.maketrans("uU", "tT"))

            key_check=DNAcheck(keydna)
            if key_check ==1:
                sys.exit("Error! Key for trimming adapters has characters besides A, T, G, C, or U!")
            if (len(key) < 20) or (len(key) > 30):
                sys.exit("Trim key must be between 20 and 30 letters")

            #Trim adapters
            process_fastqs(fastqs, keydna, args)
            summarize_cutadapt_reports()

            print("Adapter trimming complete.")

            get_distribution_from_trimmedLibraries("trimmedLibraries/")

            # library size graph:
            #R lib distribution
            r_script1 = """

            # Load data
            df <- read_csv("trimming_summary.csv")

            # Pivot into long format, including Total reads
            df_long <- df %>%
            filter(Metric %in% c( "Reads with adapters", "Reads written (passing filters)","Total reads")) %>%
            pivot_longer(-Metric, names_to = "Library", values_to = "Count") %>%
            mutate(
                Count = as.numeric(Count),
                Metric = factor(
                Metric,
                levels = c("Reads with adapters", "Reads written (passing filters)","Total reads")
                )
            )
            
            # Pivot wider to get totals for each library
            totals <- df_long %>%
            filter(Metric == "Total reads") %>%
            select(Library, Total = Count)

            # Join with long data and compute percentage
            df_labeled <- df_long %>%
            left_join(totals, by = "Library") %>%
            mutate(Percent = round(100 * Count / Total, 1),
                    Label = paste0(Percent, "%"))

            # Extract % Reads written to use as labels
            percent_labels <- df %>%
            filter(Metric == "% Reads written") %>%
            pivot_longer(-Metric, names_to = "Library", values_to = "Percent")

            # Plot function
            lib_plot <- function(data) {
            ggplot(df_labeled, aes(x = Library, y = Count, fill = Metric)) +
                geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
                        color = "black", width = 0.6, alpha = 0.85) +
                geom_text(
                    aes(label = Label),
                    position = position_dodge(width = 0.7),
                    vjust = -0.5,
                    size = 3,
                    fontface = "bold"
                ) +
                scale_fill_manual(values = c(
                    "Total reads" = "grey80",
                    "Reads written (passing filters)" = "steelblue",
                    "Reads with adapters" = "lightblue"
                )) +
                scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.1))) +
                labs(
                    title = "Cutadapt Summary",
                    x = "Library",
                    y = "Read Count",
                    fill = "Read Type"
                ) +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            }

            # Generate and save the plot
            plot_obj <- lib_plot(df_labeled)
            ggsave("Library_summary.png", plot = plot_obj, width = 10, height = 7, dpi = 300, bg = "white")

            """

            robjects.r(r_script1)
            
            #Python code:
            file_path ="read_length_distribution.csv"
            r_script="""

            file_path <- "read_length_distribution.csv"

            # Load and filter data
            data <- read.csv(file_path)
            data_filt <- data[data$Read.Length < 51, ]

            # Get list of library names (excluding Read.Length column)
            libraries <- colnames(data_filt)[-1]

            # Start PDF device
            pdf("Read_distribution.pdf", width = 10, height = 6, bg = "white")

            # Loop through libraries and plot each one
            for (lib in libraries) {
            df_single <- data_filt %>%
                select(Read.Length, all_of(lib)) %>%
                rename(Count = all_of(lib))

            plot_obj <- ggplot(df_single, aes(x = Read.Length, y = Count)) +
                geom_line(color = "steelblue", linewidth = 1) +
                theme_minimal() +
                labs(
                title = paste("Read Length Distribution:", lib),
                x = "Read Length",
                y = "Count"
                ) +
                scale_x_continuous(breaks = seq(min(df_single$Read.Length), max(df_single$Read.Length), by = 1)) +
                theme(legend.position = "none")

            print(plot_obj)  # Each plot gets its own page
            }

            # Close PDF device
            dev.off()
            """

            robjects.r(r_script)

            t_fastq_files = "trimmedLibraries/"
            if args.annotate is not None:
                run_ss(t_fastq_files,args.genome,args.known_mirnas)
            

            #cmd = f"rm -rf trimmedLibraries/*.txt"
            #os.system(cmd)
            cmd = f"rm trimming_summary.csv read_length_distribution.csv"
            os.system(cmd)
if __name__ == "__main__":
    main()

