import duckdb
import json

import numpy as np

from datetime import datetime
from utils import sha256sum

from rwmutex import RWLock


__lock = RWLock()

DATABASE = 'db.duckdb'


with __lock.write:
	db = duckdb.connect(database=DATABASE, read_only=False)
	db.execute(
	"""
	CREATE TABLE IF NOT EXISTS files (
		hash VARCHAR(40) PRIMARY KEY, 
		name VARCHAR(1000) NOT NULL,
		path VARCHAR(1000) NOT NULL,
		genome_ref VARCHAR(64) NOT NULL,
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
		var_type VARCHAR,
		var_subtype VARCHAR,

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

	CREATE SEQUENCE IF NOT EXISTS gene_sets_id_seq START 1;

	CREATE TABLE IF NOT EXISTS gene_sets (
		id UINTEGER PRIMARY KEY,
		name VARCHAR(1000) NOT NULL,
		description VARCHAR(1000) NOT NULL,
		created_at TIMESTAMP NOT NULL
	);

	CREATE SEQUENCE IF NOT EXISTS gene_set_members_id_seq START 1;

	CREATE TABLE IF NOT EXISTS gene_set_members (
		id UINTEGER PRIMARY KEY,
		gene_set_id UINTEGER NOT NULL,
		name VARCHAR(1000) NOT NULL,

		FOREIGN KEY(gene_set_id) REFERENCES gene_sets(id)
	);
	"""
	)
	db.close()


INSERT_VARIANTS_QUERY = """
INSERT INTO variants
SELECT
	'{}' AS file_hash,
	'{}' AS gene_hgnc,
	gene_variation,
	chrom,
	pos,
	id,
	ref,
	str_split(alt, ',') AS alt,
	qual,
	str_split(filter, ',') AS filter,
	info,
	format,
	start_pos,
	end_pos,
	str_split(alleles, ',') AS alleles,
	affected_start,
	affected_end,
	var_type,
	var_subtype
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
		db = duckdb.connect(database=DATABASE, read_only=False)
		db.execute('INSERT INTO genes (file_hash, gene_hgnc) VALUES (?, ?)', (file_hash, gene))

		# add index inplace as a new column and rename it gene_variation
		variants.reset_index(inplace=True)
		variants.rename({'index': 'gene_variation'}, axis=1, inplace=True)

		# join array values as csv, because duckdb doesn't parse them properly
		# so we need to split them manually
		variants['alt'] = variants['alt'].apply(lambda alt: ','.join(str(a) for a in alt))
		variants['filter'] = variants['filter'].apply(lambda filter: ','.join(filter))
		variants['alleles'] = variants['alleles'].apply(lambda alleles: ','.join(str(a) for a in alleles))
		variants['info'] = variants['info'].apply(lambda info: json.dumps(info))

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
		db = duckdb.connect(database=DATABASE, read_only=True)
		files = db.execute('SELECT * FROM files WHERE hash = ?', (sha,)).fetch_df().to_dict('records')
		file = files[0] if files else None
		db.close()
		return file


def save_file(filename, sha, path, genome_ref, created_at, status='unprocessed'):
	with __lock.write:
		db = duckdb.connect(database=DATABASE, read_only=False)
		db.execute(
			'INSERT INTO files (hash, name, path, genome_ref, created_at, status) VALUES (?, ?, ?, ?, ?, ?)',
			(sha, filename, path, genome_ref, created_at, status))
		db.close()


def update_file_status(sha, status):
	with __lock.write:
		db = duckdb.connect(database=DATABASE, read_only=False)
		db.execute('UPDATE files SET status=? WHERE hash = ?', (status, sha))
		db.close()


def get_files():
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		files =  db.execute('SELECT * FROM files').fetch_df().to_dict('records')
		db.close()
		return files


def get_chromosome_for_gene(gene_hgnc):
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		chroms =  db.execute('SELECT chrom FROM variants WHERE gene_hgnc = ? LIMIT 1', (gene_hgnc,)).fetchone()
		chrom = chroms[0] if chroms else None
		db.close()
		return chrom


def __in_filter(column, values):
	return '  AND {} IN ({})'.format(column, ','.join(['?']*len(values)))


def get_variants(sha, gene_hgnc, effects=None, impacts=None, biotypes=None, feature_types=None):
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		variants_df = None
		query = """
		SELECT DISTINCT v.gene_variation, start_pos, end_pos, ref, a.alt, var_type, var_subtype
		FROM variants v
		JOIN annotations a ON v.file_hash = a.file_hash AND v.gene_hgnc = a.gene_hgnc AND v.gene_variation = a.gene_variation
		WHERE v.file_hash = ?
			AND v.gene_hgnc = ?
        """

		if effects:
			query += __in_filter('effect', effects)
		if impacts:
			query += __in_filter('impact', impacts)
		if biotypes:
			query += __in_filter('transcript_biotype', biotypes)
		if feature_types:
			query += __in_filter('feature_type', feature_types)

		params = [sha, gene_hgnc] + [v for p in [effects, impacts, biotypes, feature_types] for v in p if p]

		variants_df =  db.execute(query, params).fetch_df()
		# convert 0-based index to 1-based and half-open interval, i.e [) to closed, i.e. []
		variants_df['start_pos'] += 1
		db.close()
		return variants_df


def delete_file(sha):
	with __lock.write:
		db = duckdb.connect(database=DATABASE, read_only=False)
		db.execute('DELETE FROM annotations WHERE file_hash = ?', (sha,))
		db.execute('DELETE FROM variants WHERE file_hash = ?', (sha,))
		db.execute('DELETE FROM genes WHERE file_hash = ?', (sha,))
		db.execute('DELETE FROM tasks WHERE file_hash = ?', (sha,))
		db.execute('DELETE FROM files WHERE hash = ?', (sha,))
		db.close()


def get_gene_sets():
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		gene_sets =  db.execute('SELECT * FROM gene_sets').fetch_df().to_dict('records')
		db.close()
		return gene_sets


def save_gene_set(name, description, genes):
	with __lock.write:
		db = duckdb.connect(database=DATABASE, read_only=False)
		cursor = db.cursor()
		insert_gene_set_query = """
		INSERT INTO gene_sets (id, name, description, created_at)
		VALUES (nextval('gene_sets_id_seq'), ?, ?, ?)
		"""
		cursor.execute(insert_gene_set_query, (name, description, datetime.now()))
		insert_gene_query = """
		INSERT INTO gene_set_members (id, name, gene_set_id)
		VALUES (nextval('gene_set_members_id_seq'), ?, currval('gene_sets_id_seq'))
		"""
		db.executemany(insert_gene_query, [[gene] for gene in genes])
		db.close()


def delete_gene_set(id):
	with __lock.write:
		db = duckdb.connect(database=DATABASE, read_only=False)
		db.execute('DELETE FROM gene_set_members WHERE gene_set_id = ?', (id,))
		db.execute('DELETE FROM gene_sets WHERE id = ?', (id,))
		db.close()


def get_gene_set_by_id(id):
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		query = 'SELECT * FROM gene_sets WHERE id = ?'
		gene_sets = db.execute(query, (id,)).fetch_df().to_dict('records')
		gene_set = gene_sets[0] if gene_sets else None
		db.close()
		return gene_set


def get_genes_for_gene_set(id):
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		query = 'SELECT * FROM gene_set_members WHERE gene_set_id = ?'
		genes = db.execute(query, (id,)).fetch_df().to_dict('records')
		db.close()
		return genes


def save_gene_set_member(name, gene_set_id):
	with __lock.write:
		db = duckdb.connect(database=DATABASE, read_only=False)
		query = """
		INSERT INTO gene_set_members (id, name, gene_set_id)
		VALUES (nextval('gene_set_members_id_seq'), ?, ?)
		"""
		db.execute(query, (name, gene_set_id))
		db.close()


def delete_gene_set_member(id):
	with __lock.write:
		db = duckdb.connect(database=DATABASE, read_only=False)
		db.execute('DELETE FROM gene_set_members WHERE id = ?', (id,))
		db.close()


def get_variant(file_hash, gene_hgnc, variant_id):
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		query = """
		SELECT *
		FROM variants
		WHERE file_hash = ?
		  AND gene_hgnc = ?
		  AND gene_variation = ?
		LIMIT 1
		"""
		variants = db.execute(query, (file_hash, gene_hgnc, variant_id)).fetch_df().to_dict('records')
		variant = None
		if variants:
			variant = variants[0]
			# convert 0-based index to 1-based and half-open interval, i.e [) to closed, i.e. []
			variant['start_pos'] += 1
		db.close()
		return variant


def get_variant_annotations(file_hash, gene_hgnc, variant_id):
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		query = """
		SELECT *
		FROM annotations
		WHERE file_hash = ?
		  AND gene_hgnc = ?
		  AND gene_variation = ?
		"""
		annotations = db.execute(query, (file_hash, gene_hgnc, variant_id)).fetch_df().to_dict('records')
		db.close()
		return annotations


def get_variant_annotation(file_hash, gene_hgnc, variant_id, annotation_id):
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		query = """
		SELECT *
		FROM annotations
		WHERE file_hash = ?
		  AND gene_hgnc = ?
		  AND gene_variation = ?
		  AND variation_annotation = ?
		LIMIT 1
		"""
		annotations = db.execute(query, (file_hash, gene_hgnc, variant_id, annotation_id)).fetch_df().to_dict('records')
		annotation = annotations[0] if annotations else None
		db.close()
		return annotation


def get_transcripts_for_variant(file_hash, gene_hgnc, variation_id):
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		query = """
		SELECT feature_id
		FROM annotations
		WHERE file_hash = ?
		  AND gene_hgnc = ?
		  AND gene_variation = ?
		  AND feature_type = 'transcript'
		"""
		genes = db.execute(query, (id,)).fetch_df().to_dict('records')
		db.close()
		return genes


def read_query(query, params):
	with __lock.read:
		db = duckdb.connect(database=DATABASE, read_only=True)
		df = db.execute(query, params).fetch_df()
		db.close()
		return df