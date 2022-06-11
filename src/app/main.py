import os
import threading
from datetime import datetime
from flask import Flask
from flask import Blueprint
from flask import render_template
from flask import flash
from flask import request
from flask import redirect
from flask import url_for
from werkzeug.utils import secure_filename

import db
from utils import sha256sum
from tasks import parse

main = Blueprint('main', __name__)


UPLOAD_FOLDER = 'uploads'
ALLOWED_EXTENSIONS = {'vcf'}

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER


@main.route("/")
def files():
    files = db.get_files()
    return render_template('files.html', files=files)


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@main.route('/files/new', methods=['GET'])
def upload_file_page():
    return render_template('upload_file.html')


@main.route('/files/new', methods=['POST'])
def upload_file():
    # check if the post request has the file part
    if 'file' not in request.files:
        flash('No file part.')
        return redirect(request.url)
    file = request.files['file']
    # If the user does not select a file, the browser submits an
    # empty file without a filename.
    if file.filename == '':
        flash('No selected file.')
        return redirect(request.url)
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(path)
        vcf_sha = sha256sum(path)
        db.save_file(filename, vcf_sha, path, datetime.now())
        # TODO: support changing the genes files
        th = threading.Thread(target=lambda: parse(path, 'data/genes.csv'))
        th.start()
        flash('File uploaded. Will start processing it in the background now. This may take a couple of minutes depending on the size of the file.')
        return redirect(url_for('main.files'))

    