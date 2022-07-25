# gene_variants
Integrated toolkit with web GUI for analysing gene variants and their potential effects.

It uses [SnpEff](https://github.com/pcingola/SnpEff) for annotating the VCFs.

## Setup

Download the repository:

```bash
$ git clone https://github.com/mzdravkov/gene_variants.git
```

Alternatively, you can download a zip archive by visiting https://github.com/mzdravkov/gene_variants then clicking on `Code->Download ZIP`. Then decompress it.


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


> **Note**: If the installation directory of SnpEff is not `gene_variants/snpEff` you'll have to change the property `snpEff_path` in `config.yaml` to the correct directory

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

> **Note**: The first time you process a VCF that uses a reference genome that hasn't been used before, SnpEff will automatically download it. This may take a few minutes.