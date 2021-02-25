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
    usa = read_covid_dna_fasta_file("data/usa_apr2020.fasta")

    add_n(wuhan)
    add_n(usa)

    # Gene ORF1ab
    wuhan_orf1ab = wuhan["ORF1ab"][0].transcribe()
    wuhan_orf1ab = numpfy(wuhan_orf1ab)
    wuhan_orf1ab = wuhan_orf1ab.reshape(wuhan["ORF1ab"][1])
    usa_orf1ab = usa["ORF1ab"][0].transcribe()
    usa_orf1ab = numpfy(usa_orf1ab)
    usa_orf1ab = usa_orf1ab.reshape(usa["ORF1ab"][1])
    mut_orf1ab = wuhan_orf1ab - usa_orf1ab

    # Gene S
    wuhan_s = wuhan["S"][0].transcribe()
    wuhan_s = numpfy(wuhan_s)
    wuhan_s = wuhan_s.reshape(wuhan["S"][1])
    usa_s = usa["S"][0].transcribe()
    usa_s = numpfy(usa_s)
    usa_s = usa_s.reshape(usa["S"][1])
    mut_s = wuhan_s - usa_s

    # Gene ORF3a
    wuhan_orf3a = wuhan["ORF3a"][0].transcribe()
    wuhan_orf3a = numpfy(wuhan_orf3a)
    wuhan_orf3a = wuhan_orf3a.reshape(wuhan["ORF3a"][1])
    usa_orf3a = usa["ORF3a"][0].transcribe()
    usa_orf3a = numpfy(usa_orf3a)
    usa_orf3a = usa_orf3a.reshape(usa["ORF3a"][1])
    mut_orf3a = wuhan_orf3a - usa_orf3a

    # Gene E
    wuhan_e = wuhan["E"][0].transcribe()
    wuhan_e = numpfy(wuhan_e)
    wuhan_e = wuhan_e.reshape(wuhan["E"][1])
    usa_e = usa["E"][0].transcribe()
    usa_e = numpfy(usa_e)
    usa_e = usa_e.reshape(usa["E"][1])
    mut_e = wuhan_e - usa_e

    # Gene M
    wuhan_m = wuhan["M"][0].transcribe()
    wuhan_m = numpfy(wuhan_m)
    wuhan_m = wuhan_m.reshape(wuhan["M"][1])
    usa_m = usa["M"][0].transcribe()
    usa_m = numpfy(usa_m)
    usa_m = usa_m.reshape(usa["M"][1])
    mut_m = wuhan_m - usa_m

    # Gene ORF6
    wuhan_orf6 = wuhan["ORF6"][0].transcribe()
    wuhan_orf6 = numpfy(wuhan_orf6)
    wuhan_orf6 = wuhan_orf6.reshape(wuhan["ORF6"][1])
    usa_orf6 = usa["ORF6"][0].transcribe()
    usa_orf6 = numpfy(usa_orf6)
    usa_orf6 = usa_orf6.reshape(usa["ORF6"][1])
    mut_orf6 = wuhan_orf6 - usa_orf6

    # Gene ORF7a
    wuhan_orf7a = wuhan["ORF7a"][0].transcribe()
    wuhan_orf7a = numpfy(wuhan_orf7a)
    wuhan_orf7a = wuhan_orf7a.reshape(wuhan["ORF7a"][1])
    usa_orf7a = usa["ORF7a"][0].transcribe()
    usa_orf7a = numpfy(usa_orf7a)
    usa_orf7a = usa_orf7a.reshape(usa["ORF7a"][1])
    mut_orf7a = wuhan_orf7a - usa_orf7a

    # Gene ORF7b
    wuhan_orf7b = wuhan["ORF7b"][0].transcribe()
    wuhan_orf7b = numpfy(wuhan_orf7b)
    wuhan_orf7b = wuhan_orf7b.reshape(wuhan["ORF7b"][1])
    usa_orf7b = usa["ORF7b"][0].transcribe()
    usa_orf7b = numpfy(usa_orf7b)
    usa_orf7b = usa_orf7b.reshape(usa["ORF7b"][1])
    mut_orf7b = wuhan_orf7b - usa_orf7b

    # Gene ORF8
    wuhan_orf8 = wuhan["ORF8"][0].transcribe()
    wuhan_orf8 = numpfy(wuhan_orf8)
    wuhan_orf8 = wuhan_orf8.reshape(wuhan["ORF8"][1])
    usa_orf8 = usa["ORF8"][0].transcribe()
    usa_orf8 = numpfy(usa_orf8)
    usa_orf8 = usa_orf8.reshape(usa["ORF8"][1])
    mut_orf8 = wuhan_orf8 - usa_orf8

    # Gene N
    wuhan_n = wuhan["N"][0].transcribe()
    wuhan_n = numpfy(wuhan_n)
    wuhan_n = wuhan_n.reshape(wuhan["N"][1])
    usa_n = usa["N"][0].transcribe()
    usa_n = numpfy(usa_n)
    usa_n = usa_n.reshape(usa["N"][1])
    mut_n = wuhan_n - usa_n

    # Gene ORF10
    wuhan_orf10 = wuhan["ORF10"][0].transcribe()
    wuhan_orf10 = numpfy(wuhan_orf10)
    wuhan_orf10 = wuhan_orf10.reshape(wuhan["ORF10"][1])
    usa_orf10 = usa["ORF10"][0].transcribe()
    usa_orf10 = numpfy(usa_orf10)
    usa_orf10 = usa_orf10.reshape(usa["ORF10"][1])
    mut_orf10 = wuhan_orf10 - usa_orf10

    f, axs = plt.subplots(11, 3, figsize=(25, 30))

    ax = axs[0, 0]
    ax.pcolor(wuhan_orf1ab)
    ax.set_title("Wuhan - Gene ORF1ab")

    ax = axs[0, 1]
    ax.pcolor(usa_orf1ab)
    ax.set_title("USA - Gene ORF1ab")

    ax = axs[0, 2]
    ax.pcolor(mut_orf1ab)
    ax.set_title("Mutation - Gene ORF1ab")

    ax = axs[1, 0]
    ax.pcolor(wuhan_s)
    ax.set_title("Wuhan - Gene S")

    ax = axs[1, 1]
    ax.pcolor(usa_s)
    ax.set_title("USA - Gene S")

    ax = axs[1, 2]
    ax.pcolor(mut_s)
    ax.set_title("Mutation - Gene S")

    ax = axs[2, 0]
    ax.pcolor(wuhan_orf3a)
    ax.set_title("Wuhan - Gene ORF3a")

    ax = axs[2, 1]
    ax.pcolor(usa_orf3a)
    ax.set_title("USA - Gene ORF3a")

    ax = axs[2, 2]
    ax.pcolor(mut_orf3a)
    ax.set_title("Mutation - Gene ORF3a")

    ax = axs[3, 0]
    ax.pcolor(wuhan_e)
    ax.set_title("Wuhan - Gene E")

    ax = axs[3, 1]
    ax.pcolor(usa_e)
    ax.set_title("USA - Gene E")

    ax = axs[3, 2]
    ax.pcolor(mut_e)
    ax.set_title("Mutation - Gene E")

    ax = axs[4, 0]
    ax.pcolor(wuhan_m)
    ax.set_title("Wuhan - Gene M")

    ax = axs[4, 1]
    ax.pcolor(usa_m)
    ax.set_title("USA - Gene M")

    ax = axs[4, 2]
    ax.pcolor(mut_m)
    ax.set_title("Mutation - Gene M")

    ax = axs[5, 0]
    ax.pcolor(wuhan_orf6)
    ax.set_title("Wuhan - Gene ORF6")

    ax = axs[5, 1]
    ax.pcolor(usa_orf6)
    ax.set_title("USA - Gene ORF6")

    ax = axs[5, 2]
    ax.pcolor(mut_orf6)
    ax.set_title("Mutation - Gene ORF6")

    ax = axs[6, 0]
    ax.pcolor(wuhan_orf7a)
    ax.set_title("Wuhan - Gene ORF7a")

    ax = axs[6, 1]
    ax.pcolor(usa_orf7a)
    ax.set_title("USA - Gene ORF7a")

    ax = axs[6, 2]
    ax.pcolor(mut_orf7a)
    ax.set_title("Mutation - Gene ORF7a")

    ax = axs[7, 0]
    ax.pcolor(wuhan_orf7b)
    ax.set_title("Wuhan - Gene ORF7b")

    ax = axs[7, 1]
    ax.pcolor(usa_orf7b)
    ax.set_title("USA - Gene ORF7b")

    ax = axs[7, 2]
    ax.pcolor(mut_orf7b)
    ax.set_title("Mutation - Gene ORF7b")

    ax = axs[8, 0]
    ax.pcolor(wuhan_orf8)
    ax.set_title("Wuhan - Gene ORF8")

    ax = axs[8, 1]
    ax.pcolor(usa_orf8)
    ax.set_title("USA - Gene ORF8")

    ax = axs[8, 2]
    ax.pcolor(mut_orf8)
    ax.set_title("Mutation - Gene ORF8")

    ax = axs[9, 0]
    ax.pcolor(wuhan_n)
    ax.set_title("Wuhan - Gene N")

    ax = axs[9, 1]
    ax.pcolor(usa_n)
    ax.set_title("USA - Gene N")

    ax = axs[9, 2]
    ax.pcolor(mut_n)
    ax.set_title("Mutation - Gene N")

    ax = axs[10, 0]
    ax.pcolor(wuhan_orf10)
    ax.set_title("Wuhan - Gene ORF10")

    ax = axs[10, 1]
    ax.pcolor(usa_orf10)
    ax.set_title("USA - Gene ORF10")

    ax = axs[10, 2]
    ax.pcolor(mut_orf10)
    ax.set_title("Mutation - Gene ORF10")

    f.tight_layout()
    f.savefig("wuhan.jpg")
