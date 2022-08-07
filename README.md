# Polymorpheus
Integrated toolkit with web GUI for analysing genetic polymorphisms and their potential effects.

The software takes a VCF file with the genetic variants, annotates it and allows browsing the results in a user-friendly manner.

It uses [SnpEff](https://github.com/pcingola/SnpEff) for annotating the VCFs.

## Setup

Download the repository:

```bash
$ git clone https://github.com/mzdravkov/polymorpheus.git
```

Alternatively, you can download a zip archive by visiting https://github.com/mzdravkov/polymorpheus then clicking on `Code->Download ZIP`. Then decompress it.


Go to the root directory of the repository:

```bash
$ cd gene_variants
```

Download [SnpEff](https://pcingola.github.io/SnpEff/) (If you already have SnpEff, you can use your existing installation. See the note below.):
```bash
# Download the latest version
$ wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip -O SnpEff.zip

# Unzip the files
$ unzip SnpEff.zip
```


> **Note**: If the installation directory of SnpEff is not `polymorpheus/snpEff` you'll have to change the property `snpEff_path` in `config.yaml` to the correct directory

Create a virtual environment and install dependencies:

```bash
$ python -m venv env
$ source env/bin/activate
$ pip install -r requirements.txt
```

## Running

To run the web app:

```bash
$ python src/run_app.py
```

Then visit the application on http://127.0.0.1:5000/.

## Usage

### Gene sets
Polymorpheus uses the concept of gene sets to specify which genes you're interested in analysing. You can define a gene set by uploading a file with HGNC gene names.
Subsequently, when you upload a VCF file you specify which gene set you wish to analyse it against. Polymorpheus will extract polymorphisms from the VCF only for the genes in the gene set.

### Browsing parsed VCFs
When you upload a VCF file, Polymorpheus will use SnpEff to annotate it. 

> **Note**: The first time you process a VCF that uses a reference genome that hasn't been used before, SnpEff will automatically download it. This may take a few minutes.

Once the VCF is marked as processed then Polymorpheus has already annotated, filtered, parsed and stored the data in a local database for fast access to the data. You can browse the information for the polymorphisms in the VCF, this includes:
- Listing affected genes with information on their function and processes they participate in.
- Listing genetic variants with their predicted effect for a given gene. You can filter them by various properties.
- Embedded genomic browser, showing the variants, transcripts and genes.
- Finding how a polymorphism modifies the protein sequence.
