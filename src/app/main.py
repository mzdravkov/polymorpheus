import os
import threading
import json
import functools
from datetime import datetime
from flask import Flask
from flask import Blueprint
from flask import render_template
from flask import flash
from flask import request
from flask import redirect
from flask import url_for
from flask import Markup
from flask import send_from_directory
from werkzeug.utils import secure_filename
import pandas as pd

import db
import analysis
import utils
import proteins
from tasks import parse
from external import get_hgnc_info
from external import get_protein_seq_from_transcript_id
from external import get_protein_annotation_from_nextprot
from vcf_processing import VCFParsingException
from vcf_processing import get_header_lines
from vcf_processing import validate_vcf_version
from vcf_processing import validate_and_get_genome_reference

main = Blueprint('main', __name__)

GENES_FILE = 'data/genes.csv'

UPLOAD_FOLDER = 'uploads'
VCF_EXTENSIONS = {'vcf', 'vcf.gz'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@main.route("/")
def files():
    files = db.get_files()
    return render_template('files.html', files=files)


def allowed_file(filename, extensios):
    for ext in extensios:
        if filename.endswith('.' + ext):
            return True
    return False


def is_vcf(filename):
    return allowed_file(filename, VCF_EXTENSIONS)


def validate_file_upload(fn):
  @functools.wraps(fn)
  def decorated_function(*args, **kwargs):
    # check if the post request has the file part
    if 'file' not in request.files:
        flash('No file part.', category='danger')
        return redirect(request.url)
    file = request.files['file']

    # If the user does not select a file, the browser submits an
    # empty file without a filename.
    if file.filename == '':
        flash('No selected file.', category='danger')
        return redirect(request.url)
    return fn(*args, **kwargs)
  return decorated_function


@main.route('/files/new', methods=['GET'])
def upload_vcf_page():
    gene_sets = db.get_gene_sets()
    return render_template('upload_vcf.html', gene_sets=gene_sets)


@main.route('/files/new', methods=['POST'])
@validate_file_upload
def upload_vcf():
    file = request.files['file']

    if not is_vcf(file.filename):
        flash('Please upload a VCF file.', category='danger')
        return redirect(request.url)

    if file:
        filename = secure_filename(file.filename)
        path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(path)

        vcf_sha = utils.sha256sum(path)
        gene_set_id = request.form['gene_set']

        existing_row = db.get_file(vcf_sha, gene_set_id)
        if existing_row:
            flash(Markup('File has already been uploaded and processed. <a href="/files/{}">Link to existing file</a>'.format(vcf_sha)), category='info')
            return redirect(url_for('main.files'))

        try:
            header = get_header_lines(path)
            validate_vcf_version(header)
            reference_genome = validate_and_get_genome_reference(header)
        except VCFParsingException as e:
            flash('Cannot process the VCF file: ' + str(e), category='danger')
            return redirect(url_for('main.files'))

        db.save_file(filename, vcf_sha, path, reference_genome, gene_set_id, datetime.now())

        genes = db.get_genes_for_gene_set(gene_set_id)
        utils.save_genes_to_file(genes, GENES_FILE)

        th = threading.Thread(target=lambda: parse(path, GENES_FILE, gene_set_id))
        th.start()
        flash('File uploaded. Will start processing it in the background now. This may take a couple of minutes depending on the size of the file.', category='success')
        return redirect(url_for('main.files'))


@main.route('/files/<sha>/<gene_set_id>')    
def file_summary(sha, gene_set_id):
    selected_chromosomes = request.args.getlist('chromosomes')
    file = db.get_file(sha, gene_set_id)
    gene_set = db.get_gene_set_by_id(file['gene_set_id'])
    file_summary = analysis.file_summary(sha).to_dict('records')[0]
    impact_summary = analysis.impact_summary(sha).to_dict('records')
    chromosomes = list({row['chrom']: None for row in impact_summary})

    if selected_chromosomes:
        impact_summary = [x for x in impact_summary if x['chrom'] in selected_chromosomes]

    return render_template(
        'file.html',
        file=file,
        gene_set=gene_set,
        file_summary=file_summary,
        impact_summary=impact_summary,
        chromosomes=chromosomes,
        selected_chromosomes=selected_chromosomes)


@main.route('/files/<sha>/<gene_set_id>/delete')
def delete_file(sha, gene_set_id):
    db.delete_file(sha, gene_set_id)
    flash('The file and all its information was deleted.', category='success')
    return redirect(url_for('main.files'))


@main.route('/files/<sha>/<gene_set_id>/gene/<gene_hgnc>')
def get_gene(sha, gene_set_id, gene_hgnc):
    file = db.get_file(sha, gene_set_id)
    gene_set = db.get_gene_set_by_id(gene_set_id)
    chromosome = db.get_chromosome_for_gene(gene_hgnc)
    hgnc_info = get_hgnc_info(gene_hgnc)
    selected_biotypes = request.args.getlist('biotypes')
    effects_summary = analysis.effects_by_impact_summary_for_gene(
        sha,
        gene_hgnc,
        biotypes=selected_biotypes
        ).to_dict('records')
    ordering = {
        'high': 4,
        'moderate': 3,
        'low': 2,
        'modifier': 1,
    }
    for row in effects_summary:
        row['order'] = ordering[row['impact'].lower()]
    effects_summary.sort(reverse=True, key=lambda row: row['order'])

    transcript_biotypes = analysis.get_transcript_biotypes(sha, gene_hgnc)

    protein_annotation = get_protein_annotation_from_nextprot(hgnc_info['uniprot_ids'][0])

    return render_template(
        'gene.html',
        file=file,
        gene_set=gene_set,
        gene_hgnc=gene_hgnc,
        hgnc_info=hgnc_info,
        chromosome=chromosome,
        transcript_biotypes=transcript_biotypes,
        selected_biotypes=selected_biotypes,
        effects_summary=effects_summary,
        protein_annotation=protein_annotation)


@main.route('/files/<file_hash>/<gene_set_id>/<gene_hgnc>/variants')
def get_gene_variants(file_hash, gene_set_id, gene_hgnc):
    file = db.get_file(file_hash, gene_set_id)
    gene_set = db.get_gene_set_by_id(gene_set_id)
    hgnc_info = get_hgnc_info(gene_hgnc)

    selected_biotypes = request.args.getlist('biotypes')
    selected_effects = request.args.getlist('effects')
    selected_impacts = request.args.getlist('impacts')
    selected_feature_types = request.args.getlist('feature_types')

    variants_df = db.get_variants(
        file_hash,
        gene_set_id,
        gene_hgnc,
        biotypes=selected_biotypes,
        effects=selected_effects,
        impacts=selected_impacts,
        feature_types=selected_feature_types
    )

    transcript_biotypes = analysis.get_transcript_biotypes(file_hash, gene_hgnc)
    effects = analysis.get_effects(file_hash, gene_hgnc)
    impacts = analysis.get_impacts(file_hash, gene_hgnc)
    feature_types = analysis.get_feature_types(file_hash, gene_hgnc)

    chromosome = db.get_chromosome_for_gene(gene_hgnc)
    min_variant_pos = variants_df['start_pos'].min()
    max_variant_pos = variants_df['end_pos'].max()
    distance = max_variant_pos - min_variant_pos
    start_pos = max(min_variant_pos - 0.1*distance, 0)
    end_pos = max_variant_pos + 0.1*distance

    variants = variants_df.to_dict('records')

    # info = pd.json_normalize(variants['info'].apply(lambda i: json.loads(i)), max_level=1)
    # variants = pd.concat([variants.drop(['info'], axis=1), info], axis=1)
    return render_template(
        'variants.html',
        file=file,
        gene_set=gene_set,
        hgnc_info=hgnc_info,
        gene_hgnc=gene_hgnc,
        variants=variants,
        transcript_biotypes=transcript_biotypes,
        selected_biotypes=selected_biotypes,
        effects=effects,
        selected_effects=selected_effects,
        impacts=impacts,
        selected_impacts=selected_impacts,
        feature_types=feature_types,
        selected_feature_types=selected_feature_types,
        chromosome=chromosome,
        start_pos=start_pos,
        end_pos=end_pos)


@main.route('/files/<file_hash>/<gene_set_id>/<gene_hgnc>/variants/<variant_id>')
def show_variant(file_hash, gene_set_id, gene_hgnc, variant_id):
    file = db.get_file(file_hash, gene_set_id)
    gene_set = db.get_gene_set_by_id(gene_set_id)
    variant = db.get_variant(file_hash, gene_set_id, gene_hgnc, variant_id)
    annotations = db.get_variant_annotations(file_hash, gene_set_id, gene_hgnc, variant_id)

    return render_template(
        'variant.html',
        file=file,
        gene_set=gene_set,
        gene_hgnc=gene_hgnc,
        variant=variant,
        annotations=annotations)


@main.route('/files/<file_hash>/<gene_set_id>/<gene_hgnc>/variants/<variant_id>/effects/<annotation_id>')
def show_effect(file_hash, gene_set_id, gene_hgnc, variant_id, annotation_id):
    file = db.get_file(file_hash, gene_set_id)
    gene_set = db.get_gene_set_by_id(gene_set_id)
    variant = db.get_variant(file_hash, gene_set_id, gene_hgnc, variant_id)
    annotation = db.get_variant_annotation(file_hash, gene_set_id, gene_hgnc, variant_id, annotation_id)
    transcript_id = annotation['feature_id'][0:15]
    ref_protein = get_protein_seq_from_transcript_id(transcript_id)
    hgvs = annotation['hgvs_protein']

    # If there's no HGVS, then there's no expected change in the protein sequence.
    alt_protein = ref_protein
    protein_change_range = (-1, -1)
    if hgvs:
        alt_protein = proteins.get_protein_variant(ref_protein, hgvs)
        protein_change_range = proteins.get_range_from_hgvs(hgvs)

    return render_template(
        'effect.html',
        file=file,
        gene_set=gene_set,
        gene_hgnc=gene_hgnc,
        variant=variant,
        annotation=annotation,
        reference_protein=ref_protein,
        alternative_protein=alt_protein,
        protein_change_range=protein_change_range)


@main.route('/files/<sha>/<gene_set_id>/<gene_hgnc>/vcf')
def get_gene_vcf(sha, gene_set_id, gene_hgnc):
    file = db.get_file(sha, gene_set_id)
    directory = os.path.join(main.root_path, '..', '..', 'data', 'intermediary', file['name'])

    include_modifiers = request.args.get('include_modifiers', default=False, type=bool)

    file_name = gene_hgnc + '_filtered.vcf.gz'
    if include_modifiers:
        file_name = gene_hgnc + '.vcf.gz'

    return send_from_directory(directory, file_name, as_attachment=True, attachment_filename=file_name)


@main.route('/files/<sha>/<gene_set_id>/<gene_hgnc>/index')
def get_gene_index(sha, gene_set_id, gene_hgnc):
    file = db.get_file(sha, gene_set_id)
    directory = os.path.join(main.root_path, '..', '..', 'data', 'intermediary', file['name'])

    include_modifiers = request.args.get('include_modifiers', default=False, type=bool)

    file_name = gene_hgnc + '_filtered.vcf.gz.tbi'
    if include_modifiers:
        file_name = gene_hgnc + '.vcf.gz.tbi'

    return send_from_directory(directory, file_name, as_attachment=True, attachment_filename=file_name)


@main.route('/gene_sets')
def list_gene_sets():
    gene_sets = db.get_gene_sets()
    return render_template('gene_sets.html', gene_sets=gene_sets)


@main.route('/gene_sets/new', methods=['GET'])
def upload_gene_set_page():
    return render_template('upload_gene_set.html')


@main.route('/gene_sets/new', methods=['POST'])
@validate_file_upload
def upload_gene_set():
    file = request.files['file']

    name = request.form['name']
    description = request.form['description']
    genes = utils.get_genes_from_file(file)

    db.save_gene_set(name, description, genes)

    flash('The new gene set was created.', category='success')
    return redirect(url_for('main.list_gene_sets'))


@main.route('/gene_sets/<id>/delete')
def delete_gene_set(id):
    db.delete_gene_set(id)
    flash('The gene set was deleted.', category='success')
    return redirect(url_for('main.list_gene_sets'))


@main.route('/gene_sets/<id>')    
def show_gene_set(id):
    gene_set = db.get_gene_set_by_id(id)
    genes = db.get_genes_for_gene_set(id)

    return render_template(
        'gene_set.html',
        gene_set=gene_set,
        genes=genes)


@main.route('/gene_sets/<id>/add_gene')    
def add_gene_to_gene_set_page(id):
    gene_set = db.get_gene_set_by_id(id)

    return render_template('gene_set_add_gene.html', gene_set=gene_set)


@main.route('/gene_sets/<id>/add_gene', methods=["POST"])    
def add_gene_to_gene_set(id):
    gene_set = db.get_gene_set_by_id(id)
    name = request.form['name']

    db.save_gene_set_member(name, gene_set['id'])

    flash('The gene was added to the dataset', category='success')
    return redirect(url_for('main.show_gene_set', id=gene_set['id']))


@main.route('/gene_sets/<gene_set_id>/remove_gene/<member_id>')
def delete_gene_set_member(gene_set_id, member_id):
    db.delete_gene_set_member(member_id)
    flash('The gene was removed from the gene set.', category='success')
    return redirect(url_for('main.show_gene_set', id=gene_set_id))


@main.route('/gencode40')
def get_gencode40():
    file_name = 'gencode.v40.annotation.sorted.gtf.gz'
    directory = os.path.join(main.root_path, '..', '..', 'data')
    return send_from_directory(directory, file_name, as_attachment=True, attachment_filename=file_name)


@main.route('/gencode40_index')
def get_gencode40_index():
    file_name = 'gencode.v40.annotation.sorted.gtf.gz.tbi'
    directory = os.path.join(main.root_path, '..', '..', 'data')
    return send_from_directory(directory, file_name, as_attachment=True, attachment_filename=file_name)