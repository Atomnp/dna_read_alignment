import os
import requests


def download_fasta(accession, filename, out_dir):
    url = (
        f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi"
        f"?id={accession}&db=nuccore&report=fasta&retmode=txt"
    )
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, filename)
    print(f"Downloading {out_path} ...")
    r = requests.get(url)
    r.raise_for_status()

    with open(out_path, "w") as f:
        f.write(r.text)

    print(f"Saved: {out_path}\n")


genomes = [
    # (Accession, Filename)
    #     ("NC_001422.1", "phix174.fasta"),                 # PhiX174
    #     ("MN908947.3", "sars_cov2.fasta"),               # SARS-CoV-2
    #     ("NC_001416.1", "lambda_phage.fasta"),           # Lambda phage
    #     ("NC_000908.2", "mycoplasma_genitalium.fasta"),  # Mycoplasma genitalium
    #     ("NC_001133.9", "yeast_chrI.fasta"),             # Yeast chr I (chromosomal FASTA)
    #     ("U00096.3", "ecoli_k12.fasta"),                 # E. coli K-12 MG1655
    #     ("NC_003070.9", "arabidopsis.fasta"),            # Arabidopsis thaliana (Chr 1)
    #     ("NC_004354.4", "drosophila_chr2L.fasta"),       # Drosophila melanogaster (Chr 2L)
    ("NC_003070.9", "arabidopsis_chr1.fasta"),  # Arabidopsis chromosome 1 (~30 Mbp)
    ("NC_004354.4", "drosophila_chr2L.fasta"),  # Drosophila chromosome 2L (~23 Mbp)
    ("NC_006088.5", "chicken_chr1.fasta"),
]


out_dir = os.path.join(os.path.dirname(__file__), "..", "genome_data")
out_dir = os.path.abspath(out_dir)

for accession, filename in genomes:
    download_fasta(accession, filename, out_dir)
