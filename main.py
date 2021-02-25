from Bio import SeqIO
from io import StringIO
import numpy as np
import matplotlib.pyplot as plt

def read_covid_dna_fasta_file(_file):
    dna_seq_record_iterator = SeqIO.parse(_file, "fasta")
    dna_seq = next(dna_seq_record_iterator)

    dna = {
        "ORF1ab": [dna_seq.seq[265:13483], (115, 115), 7],
        "S": [dna_seq.seq[21562:25384], (62, 62), 22],
        "ORF3a": [dna_seq.seq[25392:26220], (28, 30), 12],
        "E": [dna_seq.seq[26244:26472], (15, 16), 12],
        "M": [dna_seq.seq[26522:27191], (26, 27), 33],
        "ORF6": [dna_seq.seq[27201:27387], (14, 14), 10],
        "ORF7a": [dna_seq.seq[27393:27759], (19, 20), 14],
        "ORF7b": [dna_seq.seq[27755:27887], (12, 12), 12],
        "ORF8": [dna_seq.seq[27893:28259], (19, 20), 14],
        "N": [dna_seq.seq[28273:29533], (36, 36), 36],
        "ORF10": [dna_seq.seq[29557:29674], (11, 11), 4]
    }

    return dna

def add_n(genome):
    genome_keys = list(genome.keys())
    for k in genome_keys:
        n = genome[k][2]
        for i in range(0, n):
            genome[k][0] += 'N'

def numpfy(gene):
    arr = ""
    for i in gene:
        if i == "G":
            arr += "0 "
        if i == "C":
            arr += "255 "
        if i == "U":
            arr += "100 "
        if i == "A":
            arr += "200 "
        if i == "N":
            arr += "75 "

    arr_np = np.fromstring(arr, dtype=np.uint8, sep=" ")
    return arr_np

if __name__ == "__main__":
    wuhan = read_covid_dna_fasta_file("data/wuhan_jan2020.fasta")
    saudi = read_covid_dna_fasta_file("data/saudi_feb2020.fasta")

    add_n(wuhan)
    add_n(saudi)

    dna = { 'gene=ORF1ab':[(115,115),7],
            'gene=S':[(62,62),22],
            'gene=ORF3a':[(28,30),12],
            'gene=E':[(15,16),12],
            'gene=M':[(26,27),33],
            'gene=ORF6':[(14,14),10],
            'gene=ORF7a':[(19,20),14],
            'gene=ORF7b':[(12,12),12],
            'gene=ORF8':[(19,20),14],
            'gene=N':[(36,36),36],
            'gene=ORF10':[(11,11),4]}

    gene_name = list(dna.keys())
    f, ax = plt.subplots(11, 3, figsize=(25, 30))
    row = 0
    col = 0
    mut_dict={}
    for i in gene_name:
        G = i[5:]
        wuhan_gene = wuhan[G][0].transcribe()
        wuhan_gene = numpfy(wuhan_gene)
        wuhan_gene = wuhan_gene.reshape(wuhan[G][1])
        saudi_gene = saudi[G][0].transcribe()
        saudi_gene = numpfy(saudi_gene)
        saudi_gene = saudi_gene.reshape(saudi[G][1])
        mut = wuhan_gene - saudi_gene

        # sub-plot the numpy array with matplotlib pcolor method.
        ax[row][col].pcolor(wuhan_gene)
        ax[row][col].set_title(G+' Gene - Wuhan')
        col+=1

        ax[row][col].pcolor(saudi_gene)
        ax[row][col].set_title(G+' Gene - Saudi')
        col+=1

        # title to subplot
        ax[row][col].pcolor(mut)
        ax[row][col].set_title(G+' Gene - Mutataion')
        row+= 1
        col=0

    f.tight_layout()
    f.savefig("wuhan.jpg")
