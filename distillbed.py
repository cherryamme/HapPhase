from Bio import SeqIO
import re

def find_repeats(sequence, min_length=4):
    # Regular expression pattern for finding repeated bases
    pattern = re.compile(r"(A{4,}|T{4,}|C{4,}|G{4,}|N{4,})")
    matches = pattern.finditer(str(sequence))
    return [(match.start(), match.end()) for match in matches]

def write_bed(chrom, repeats, output_file):
    with open(output_file, "a") as bed_file:  # use append mode
        for start, end in repeats:
            bed_file.write(f"{chrom}\t{start}\t{end}\n")

def main(fasta_file, bed_output_file):
    # Parse the FASTA file using SeqIO
    for record in SeqIO.parse(fasta_file, "fasta"):
        chrom = record.id
        sequence = record.seq
        if chrom == "chr5":
            print("Process in chr5")
            repeats = find_repeats(sequence)
            write_bed(chrom, repeats, bed_output_file)

# Usage example
fasta_file = "/home/jiangchen/00_software/reference/GRCh38/hg38.fa"
bed_output_file = "repeats.bed"
main(fasta_file, bed_output_file)
