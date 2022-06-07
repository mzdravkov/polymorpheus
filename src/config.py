import yaml

with open("config.yml", "r") as ymlfile:
    CONFIG = yaml.safe_load(ymlfile)