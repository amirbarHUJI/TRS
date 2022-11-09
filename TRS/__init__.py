import yaml
from collections import defaultdict

from TRS.Lib import Lib


def process_config(config_file):
    with open(config_file, "rt") as fl:
        config_dct = yaml.safe_load(fl)

    lib_list = [
        Lib(lib_dct["name"], lib_dct["group"], lib_dct["count_file"])
        for lib_dct in config_dct["libs"]
    ]

    group_dct = defaultdict(list)

    for lib in lib_list:
        group_dct[lib.Group].append(lib)

    return lib_list, group_dct
