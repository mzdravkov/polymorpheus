import duckdb

import numpy as np

from utils import sha256sum

from rwmutex import RWLock


__lock = RWLock()


with __lock.write:
	db = duckdb.connect(database='db.duckdb', read_only=False)
	db.execute(
	"""
	CREATE TABLE IF NOT EXISTS files (
		hash VARCHAR(40) PRIMARY KEY, 
		name VARCHAR(1000) NOT NULL,
		path VARCHAR(1000) NOT NULL,
		created_at TIMESTAMP NOT NULL,
		status VARCHAR(32) NOT NULL
	);

	CREATE TABLE IF NOT EXISTS tasks (
		id UINTEGER PRIMARY KEY,
		created_at TIMESTAMP, 
		file_hash VARCHAR(40) NOT NULL,

		FOREIGN KEY(file_hash) REFERENCES files(hash)
	);

	CREATE TABLE IF NOT EXISTS genes (
		file_hash VARCHAR(40) NOT NULL,
		gene_hgnc VARCHAR NOT NULL,

		PRIMARY KEY (file_hash, gene_hgnc),
		FOREIGN KEY(file_hash) REFERENCES files(hash)
	);

	CREATE TABLE IF NOT EXISTS variants (
		file_hash VARCHAR(40) NOT NULL,
		gene_hgnc VARCHAR NOT NULL,
		gene_variation UINTEGER NOT NULL,

		-- start of standard VCF-fields
		chrom VARCHAR,
		pos LONG,
		id VARCHAR,
		ref VARCHAR,
		alt VARCHAR[],
		qual DOUBLE,
		filter VARCHAR[],
		info JSON,
		format VARCHAR,

		-- start of additional PyVCF fields
		start_pos UBIGINT,
		end_pos UBIGINT,
		alleles VARCHAR[],
		affected_start UBIGINT,
		affected_end UBIGINT,

		PRIMARY KEY (file_hash, gene_hgnc, gene_variation),
		FOREIGN KEY(file_hash, gene_hgnc) REFERENCES genes(file_hash, gene_hgnc),
	);

	CREATE TABLE IF NOT EXISTS annotations (
		file_hash VARCHAR(40) NOT NULL,
		gene_hgnc VARCHAR NOT NULL,
		gene_variation UINTEGER NOT NULL,
		variation_annotation UINTEGER NOT NULL,
		
		-- start of snpEff annotation fields
		alt VARCHAR,
		effect VARCHAR,
		impact VARCHAR,
		gene VARCHAR,
		gene_id VARCHAR,
		feature_type VARCHAR,
		feature_id VARCHAR,
		transcript_biotype VARCHAR,
		rank_to_total VARCHAR,
		hgvs_dna VARCHAR,
		hgvs_protein VARCHAR,
		cdna_pos_to_cdna_len VARCHAR,
		cds_pos_to_cds_len VARCHAR,
		prot_pos_to_prot_len VARCHAR,
		distance_to_feature VARCHAR,
		note VARCHAR,

		PRIMARY KEY (file_hash, gene_hgnc, gene_variation, variation_annotation),
		FOREIGN KEY(file_hash, gene_hgnc, gene_variation) REFERENCES variants(file_hash, gene_hgnc, gene_variation),
	);
	"""
	)
	db.close()


INSERT_VARIANTS_QUERY = """
INSERT INTO variants
SELECT
	'{}' as file_hash,
	'{}' as gene_hgnc,
	gene_variation,
	chrom,
	pos,
	id,
	ref,
	str_split(alt, ',') as alt,
	qual, str_split(filter, ',') as filter,
	info,
	format,
	start_pos,
	end_pos,
	str_split(alleles, ',') as alleles,
	affected_start,
	affected_end
FROM {}
"""

INSERT_ANNOTATIONS_QUERY = """
INSERT INTO annotations
SELECT
	'{}' as file_hash,
	'{}' as gene_hgnc,
	gene_variation,
	variation_annotation,
	alt,
	effect,
	impact,
	gene,
	gene_id,
	feature_type,
	feature_id,
	transcript_biotype,
	rank_to_total,
	hgvs_dna,
	hgvs_protein,
	cdna_pos_to_cdna_len,
	cds_pos_to_cds_len,
	prot_pos_to_prot_len,
	distance_to_feature,
	note
FROM {}
"""


def save_gene_data(file_hash, gene, variants, annotations):
	with __lock.write:
		db = duckdb.connect(database='db.duckdb', read_only=False)
		db.execute('INSERT INTO genes (file_hash, gene_hgnc) VALUES (?, ?)', (file_hash, gene))

		# add index inplace as a new column and rename it gene_variation
		variants.reset_index(inplace=True)
		variants.rename({'index': 'gene_variation'}, axis=1, inplace=True)

		# join array values as csv, because duckdb doesn't parse them properly
		# so we need to split them manually
		variants['alt'] = variants['alt'].apply(lambda alt: ','.join(str(a) for a in alt))
		variants['filter'] = variants['filter'].apply(lambda filter: ','.join(filter))
		variants['alleles'] = variants['alleles'].apply(lambda alleles: ','.join(str(a) for a in alleles))

		db.register('variants_df', variants)
		db.register('annotations_df', annotations)

		insert_variants_query = INSERT_VARIANTS_QUERY.format(file_hash, gene, 'variants_df')
		db.execute(insert_variants_query)

		insert_annotations_query = INSERT_ANNOTATIONS_QUERY.format(file_hash, gene, 'annotations_df')
		db.execute(insert_annotations_query)

		# When we register the dataframes, duckdb would keep references to them.
		# We unregister them so that the memory can be freed.
		db.unregister('variants_df')
		db.unregister('annotations_df')
		db.close()


def get_file_by_sha(sha):
	with __lock.read:
		db = duckdb.connect(database='db.duckdb', read_only=True)
		file = db.execute('SELECT * FROM files WHERE hash = ?', (sha,)).fetch_df().to_dict('records')[0]
		db.close()
		return file


def save_file(filename, sha, path, created_at, status='unprocessed'):
	with __lock.write:
		db = duckdb.connect(database='db.duckdb', read_only=False)
		db.execute(
			'INSERT INTO files (hash, name, path, created_at, status) VALUES (?, ?, ?, ?, ?)',
			(sha, filename, path, created_at, status))
		db.close()


def update_file_status(sha, status):
	with __lock.write:
		db = duckdb.connect(database='db.duckdb', read_only=False)
		db.execute('UPDATE files SET status=? WHERE hash = ?', (status, sha))
		db.close()


def get_files():
	with __lock.read:
		db = duckdb.connect(database='db.duckdb', read_only=True)
		files =  db.execute('SELECT * FROM files').fetch_df().to_dict('records')
		db.close()
		return files


def get_chromosome_for_gene(gene_hgnc):
	with __lock.read:
		db = duckdb.connect(database='db.duckdb', read_only=True)
		chrom =  db.execute('SELECT chrom FROM variants WHERE gene_hgnc = ? LIMIT 1', (gene_hgnc,)).fetchone()[0]
		db.close()
		return chrom


def read_query(query, params):
	with __lock.read:
		db = duckdb.connect(database='db.duckdb', read_only=True)
		df = db.execute(query, params).fetch_df()
		db.close()
		return df