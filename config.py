from os import path

# paths
PATH = dict()
PATH['repo'] = path.dirname(path.abspath(__file__))
PATH['docs'] = path.join(PATH['repo'], 'docs')
PATH['package'] = path.join(PATH['repo'], 'k_combinat_for_sage')
PATH['dist'] = path.join(PATH['repo'], 'dist')
PATH['version'] = path.join(PATH['repo'], 'VERSION')
PATH['setup'] = path.join(PATH['repo'], 'setup.py')
PATH['namespace_deletions'] = path.join(PATH['docs'], 'namespace_deletions.py')
