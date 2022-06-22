import os
import threading
import json
from datetime import datetime
from flask import Flask
from flask import Blueprint
from flask import render_template
from flask import flash
from flask import request
from flask import redirect
from flask import url_for
from flask import Markup
from werkzeug.utils import secure_filename
import pandas as pd

import db
import analysis
from utils import sha256sum
from tasks import parse
from external import get_hgnc_info
from vcf_processing import validate_vcf, VCFParsingException

main = Blueprint('main', __name__)


UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'vcf', 'vcf.gz'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@main.route("/")
def files():
    files = db.get_files()
    return render_template('files.html', files=files)


def allowed_file(filename):
    for ext in ALLOWED_EXTENSIONS:
        if filename.endswith('.' + ext):
            return True
    return False


@main.route('/files/new', methods=['GET'])
def upload_file_page():
    return render_template('upload_file.html')


@main.route('/files/new', methods=['POST'])
def upload_file():
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

    if not allowed_file(file.filename):
        flash('Please upload a VCF file.', category='danger')
        return redirect(request.url)

    if file:
        filename = secure_filename(file.filename)
        path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(path)
        vcf_sha = sha256sum(path)
        existing_row = db.get_file_by_sha(vcf_sha)
        if existing_row:
            flash(Markup('File has already been uploaded and processed. <a href="/files/{}">Link to existing file</a>'.format(vcf_sha)), category='info')
            return redirect(url_for('main.files'))

        try:
            validate_vcf(path)
        except VCFParsingException as e:
            print(e)
            flash('Cannot process the VCF file: ' + str(e), category='danger')
            return redirect(url_for('main.files'))

        db.save_file(filename, vcf_sha, path, datetime.now())
        # TODO: support changing the genes files
        th = threading.Thread(target=lambda: parse(path, 'data/genes.csv'))
        th.start()
        flash('File uploaded. Will start processing it in the background now. This may take a couple of minutes depending on the size of the file.', category='success')
        return redirect(url_for('main.files'))


@main.route('/files/<sha>')    
def file_summary(sha):
    selected_chromosomes = request.args.getlist('chromosomes')
    file = db.get_file_by_sha(sha)
    file_summary = analysis.file_summary(sha).to_dict('records')[0]
    impact_summary = analysis.impact_summary(sha).to_dict('records')
    chromosomes = list({row['chrom']: None for row in impact_summary})

    if selected_chromosomes:
        impact_summary = [x for x in impact_summary if x['chrom'] in selected_chromosomes]

    return render_template(
        'file.html',
        file=file,
        file_summary=file_summary,
        impact_summary=impact_summary,
        chromosomes=chromosomes,
        selected_chromosomes=selected_chromosomes)


@main.route('/files/<sha>/delete')
def delete_file(sha):
    db.delete_file(sha)
    flash('The file and all its information was deleted.', category='success')
    return redirect(url_for('main.files'))


@main.route('/files/<sha>/gene/<gene_hgnc>')
def get_gene(sha, gene_hgnc):
    file = db.get_file_by_sha(sha)
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

    return render_template(
        'gene.html',
        file=file,
        gene_hgnc=gene_hgnc,
        hgnc_info=hgnc_info,
        chromosome=chromosome,
        transcript_biotypes=transcript_biotypes,
        selected_biotypes=selected_biotypes,
        effects_summary=effects_summary)


@main.route('/files/<sha>/<gene_hgnc>/variants')
def get_gene_variants(sha, gene_hgnc):
    file = db.get_file_by_sha(sha)
    hgnc_info = get_hgnc_info(gene_hgnc)

    effect = request.args.get('effect')
    impact = request.args.get('impact')

    variants = db.get_variants(sha, gene_hgnc, effect=effect, impact=impact).to_dict('records')
    # info = pd.json_normalize(variants['info'].apply(lambda i: json.loads(i)), max_level=1)
    # variants = pd.concat([variants.drop(['info'], axis=1), info], axis=1)
    return render_template(
        'variants.html',
        file=file,
        hgnc_info=hgnc_info,
        gene_hgnc=gene_hgnc,
        variants=variants)