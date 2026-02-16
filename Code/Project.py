
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re
import pandas as pd
import matplotlib.pyplot as plt
import sys

def main():

    inp = int(input(f"1.Enter FASTA file: \n2.Enter DNA sequence:\n3.Enter RNA sequence:\n4.Enter complementary DNA sequence\nChoose from the following:").strip())

    if inp == 1:
        sequence = check_Fasta(input("Enter your Fasta file: ")).strip()
    elif inp == 2:
        sequence = input("Enter DNA sequence: ").strip().upper()
    elif inp == 3:
        mrna = input("Enter the RNA sequence: ").strip().upper()
        sequence = Seq(mrna).back_transcribe()
    elif inp == 4:
        sequence = Seq(input("Enter the complementary DNA sequence:").strip().upper()).reverse_complement()
    else:
        sys.exit("Choose a valid method of input")
    df = find_crispr_targets(sequence)

    if df.empty:
        print("No CRISPR target sites found.")
    else:
        print("\nCRISPR Target Sites Found:")
        print(df)

        df.to_csv("crispr_targets.csv", index=False)
        print("\nSaved results to crispr_targets.csv")

        plot_targets(df, len(sequence))


def find_crispr_targets(sequence):

    results = []
    seq_str = str(sequence).upper()

    for match in re.finditer(r'([ACGT]{20})(GG)', seq_str):
        gRNA = match.group(1)
        pam = match.group(2)
        start = match.start()   
        end = match.end()
        gc_content = gc_fraction(gRNA) * 100

        results.append({
        "gRNA": gRNA,
        "PAM": pam,
        "Start": start,
        "End": end,
        "GC_Content(%)": gc_content
        })

    return pd.DataFrame(results)


def plot_targets(df, seq_length):

    if df.empty:
        sys.exit("No sites to visualize.")

    plt.figure(figsize=(12, 3))
    plt.scatter(df["Start"], [1] * len(df), color="red", marker="|", s=200, label="PAM sites")
    plt.title("CRISPR Target Sites (NGG PAM)")
    plt.xlabel("Sequence Position")
    plt.yticks([])
    plt.xlim(0, seq_length)
    plt.legend()
    plt.tight_layout()
    plt.savefig("crispr_targets_plot.png")
    plt.show()
    print("Visualization saved as crispr_targets_plot.png")

def check_Fasta(fasta_file):
    try:
        record = SeqIO.read(fasta_file, "fasta")
    except FileNotFoundError:
        sys.exit("File not found")

    else:
        sequence = record.seq
        return sequence


if __name__ == "__main__":
    main()
