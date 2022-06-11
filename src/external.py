import requests

HGNC_FETCH_URL = 'http://rest.genenames.org/fetch/symbol/{}'


def get_hgnc_info(gene_hgnc):
    req = requests.get(HGNC_FETCH_URL.format(gene_hgnc), headers={"Accept":"application/json"})
    return req.json()['response']['docs'][0]