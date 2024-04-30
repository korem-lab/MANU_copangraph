import glob
import os
import re
from itertools import product

CLUSTER_BASE_PATH = '/ifs/scratch/tk2829_gp/shared/Projects/COMET'
CLUSTER_BASE_PATH_DATA = os.path.join(CLUSTER_BASE_PATH, 'data')
CLUSTER_ANLS_PATHS = {
    'analysis_1': os.path.join(CLUSTER_BASE_PATH, 'analysis_1'),
    'analysis_2': os.path.join(CLUSTER_BASE_PATH, 'analysis_2'),
    'analysis_3': os.path.join(CLUSTER_BASE_PATH, 'analysis_3'),
    'analysis_4': os.path.join(CLUSTER_BASE_PATH, 'analysis_4')
}

AWS_BASE_PATH = '/home/ubuntu/COMET/'
AWS_BASE_PATH_DATA = os.path.join(AWS_BASE_PATH, 'data')
AWS_ANLS_PATHS = {
    'analysis_1': os.path.join(AWS_BASE_PATH, 'analysis_1'),
    'analysis_2': os.path.join(AWS_BASE_PATH, 'analysis_2'),
    'analysis_3': os.path.join(AWS_BASE_PATH, 'analysis_3'),
    'analysis_4': os.path.join(AWS_BASE_PATH, 'analysis_4')
}

LOCAL_BASE_PATH = '/home/izaak/data/CoPan'
LOCAL_BASE_PATH_DATA = os.path.join(LOCAL_BASE_PATH, 'data')
LOCAL_ANLS_PATHS = {
    'analysis_1': os.path.join(LOCAL_BASE_PATH, 'analysis_1'),
    'analysis_2': os.path.join(LOCAL_BASE_PATH, 'analysis_2'),
    'analysis_3': os.path.join(LOCAL_BASE_PATH, 'analysis_3'),
    'analysis_4': os.path.join(LOCAL_BASE_PATH, 'analysis_4')
}

DATASETS = [
    '5G0', '5G1', '5G2', '10G0', '10G1', '10G2'
]
DEPTHS = [
    '0.1', '0.5', '1', '2', '4', '8', '16', '32'
]
TOOLS = [
    'ms', 'mh', 'METACARVELms', 'METACARVELmh', 'COMETms', 'COMETmh', 'METACARVELms+contigs', 'METACARVELmh+contigs'
]

depth_map = lambda x: int(float(x) * 10**6)


def parse_key(key):
    if re.match('(.+)-(.+)M_(.+)', key):
        x, y, z = re.findall('(.+)-(.+)M_(.+)', key).pop()
        assert x in DATASETS and y in DEPTHS and z in TOOLS
        return x, y, z
    elif re.match('(.+)-(.+)M', key):
        x, y = re.findall('(.+)-(.+)M', key).pop()
        assert x in DATASETS and y in DEPTHS
        return x, y
    else:
        return None


def make_key(dataset, depth, tool=None):
    if tool is not None:
        return f'{dataset}-{depth}M_{tool}'
    else:
        return f'{dataset}-{depth}M'


def data_path(key, platform='cluster', raw=False):
    assert platform in ['cluster', 'aws']
    base = CLUSTER_BASE_PATH_DATA if platform == 'cluster' else AWS_BASE_PATH_DATA
    base = os.path.join(base, 'raw_data') if raw else base
    path = os.path.join(base, key)
    assert os.path.isdir(path)
    return path


def analysis_path(key, analysis, platform='cluster'):
    assert platform in ['cluster', 'aws', 'local']
    assert analysis in ['analysis_1', 'analysis_2', 'analysis_3', 'analysis_4']
    if platform == 'cluster':
        path = os.path.join(CLUSTER_ANLS_PATHS[analysis], key)
    elif platform == 'aws':
        path = os.path.join(AWS_ANLS_PATHS[analysis], key)
    else:
        path = os.path.join(LOCAL_ANLS_PATHS[analysis], key)
    print(path)
    assert os.path.isdir(path)
    return path


def generate_keys(*args, tools=True):
    """
    Generate all keys, excluding those which contain an argument in args as a substring.
    """

    # Make sure each argument is a valid substring
    for a in args:
        assert a in DATASETS or a in DEPTHS or a in TOOLS

    for dataset, depth, tool in product(DATASETS, DEPTHS, TOOLS if tools else [None]):
        if dataset in args or depth in args or tool in args:
            continue
        yield make_key(dataset, depth, tool)


def get_gffs(key, platform='cluster'):
    assert platform in ['cluster', 'aws']

    dataset, depth, tool = parse_key(key)
    base = CLUSTER_BASE_PATH_DATA if platform == 'cluster' else AWS_BASE_PATH_DATA
    anno_pattern = os.path.join(base, f'/annotations/{dataset}', '*.gff')
    return glob.glob(anno_pattern)


def get_blasts(key, analysis='analysis_2', platform='cluster'):
    assert analysis in ['analysis_2', 'analysis_3', 'analysis_4']
    dataset, depth, tool = parse_key(key)
    base = analysis_path(key, analysis, platform)
