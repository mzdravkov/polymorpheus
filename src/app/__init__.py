from flask import Flask
from babel import dates

from .main import main as main_blueprint

from config import CONFIG


def create_app():
    app = Flask(__name__)

    app.config['SECRET_KEY'] = 'very_secret-TODO-change-it'
    app.config['SERVER_NAME'] = CONFIG['hostname']

    app.register_blueprint(main_blueprint)

    @app.template_filter()
    def format_datetime(value, format='medium'):
        format="HH:mm dd.MM.y"
        return dates.format_datetime(value, format)

    @app.template_filter()
    def normalize_chromosome(chrom):
        if chrom.startswith('chr'):
            return chrom
        return 'chr' + chrom

    @app.template_filter()
    def chunkstring(string, length):
        return (string[i:length+i] for i in range(0, len(string), length))

    @app.template_filter()
    def snakecase_to_title(string):
        return string.replace('-', ' ').upper()

    return app