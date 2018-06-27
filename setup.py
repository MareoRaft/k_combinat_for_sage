from os import path
import subprocess

from distutils.core import setup

# get paths
PATH_REPO = path.dirname(path.abspath(__file__))
PATH_VERSION = path.join(PATH_REPO, 'VERSION')

# config
name = 'k_combinat_for_sage'
version = open(PATH_VERSION).read().strip()
url = 'https://github.com/mareoraft/{name}'.format(name=name)
download_url = '{url}/archive/{tag}.tar.gz'.format(url=url, tag=version)
setup(
	name=name,
	packages=[name], # must be same as name above
	version=version,
	description='k-Schur combinatorics for SageMath',
	author='Matthew Lancellotti',
	author_email='mvlancellotti@gmail.com',
	url=url, # URL to github repo
	download_url=download_url,
	keywords=['morse', 'k-boundary', 'k-rim', 'k-shape', 'skew-linked diagram', 'k-irreducible', 'root ideals', 'k-schur', 'combinatorics', 'sage', 'sagemath'],
	classifiers=[],
)
