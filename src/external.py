import requests
from collections import defaultdict

HGNC_FETCH_URL = 'http://rest.genenames.org/fetch/symbol/{}'
ENSEMBL_FETCH_PROTEIN_FROM_ENSEMBL_ID = "https://rest.ensembl.org/sequence/id/{}?type=protein;species=homo_sapiens;db_type=core"
NEXTPROT_FETCH_PROTEIN_FUNCTION = "https://api.nextprot.org/entry/{}/function"


def get_hgnc_info(gene_hgnc):
    req = requests.get(HGNC_FETCH_URL.format(gene_hgnc), headers={"Accept":"application/json"})

    if not req.ok:
        req.raise_for_status()

    return req.json()['response']['docs'][0]


def get_protein_seq_from_transcript_id(transcript_id):
    req = requests.get(
        ENSEMBL_FETCH_PROTEIN_FROM_ENSEMBL_ID.format(transcript_id),
        headers={"Accept":"application/json"})

    if not req.ok:
        req.raise_for_status()

    return req.json()['seq']


def get_protein_annotation_from_nextprot(uniprot_id):
    nextprot_id = "NX_" + uniprot_id

    req = requests.get(
        NEXTPROT_FETCH_PROTEIN_FUNCTION.format(nextprot_id),
        headers={"Accept":"application/json"})

    if not req.ok:
        req.raise_for_status()

    return req.json()['entry']['annotationsByCategory']