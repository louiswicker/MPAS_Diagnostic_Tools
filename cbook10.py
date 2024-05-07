import yaml

#=======================================================================================================
#
def read_yaml(filename):

    with open(f'{filename}', 'r') as f:
        config = yaml.safe_load(f)

    return config

#=======================================================================================================
#
def list2dict(list1, list2):

    return dict(zip(list1, list2))

