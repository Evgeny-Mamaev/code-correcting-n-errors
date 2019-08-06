import os

import yaml


def remove_files():
    config = yaml.safe_load(open("config.yml"))
    if os.path.isfile(config['coder-generator']):
        os.remove(config['coder-generator'])
    if os.path.isfile(config['decoder-parity-check']):
        os.remove(config['decoder-parity-check'])
    if os.path.isfile(config['decoder-syndrome-decoding']):
        os.remove(config['decoder-syndrome-decoding'])
