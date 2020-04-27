# PyPI distribution configuration.
# To actually BUILD/PACKAGE/SHIP TO PYPI a new version, run `build_version.py`.
from distutils.core import setup
import sys

# config
name = 'k_combinat_for_sage'
# version = WILL BE AUTOMATICALLY POPULATED IN LINE BELOW.  DO NOT EDIT LINE BELOW.
version = '0.1.5'
url = 'https://github.com/mareoraft/{name}'.format(name=name)
download_url = '{url}/archive/{tag}.tar.gz'.format(url=url, tag=version)
setup(
	name=name,
	packages=[name], # must be same as name above
	version=version,
	description='k-Schur combinatorics for SageMath',
	author='Matthew Lancellotti and George H. Seelinger',
	author_email='mvlancellotti@gmail.com',
	url=url, # URL to github repo
	download_url=download_url,
	keywords=['morse', 'k-boundary', 'k-rim', 'k-shape', 'skew-linked diagram', 'k-irreducible', 'root ideals', 'k-schur', 'combinatorics', 'sage', 'sagemath', 'raising root operators', 'catalan functions'],
	classifiers=[],
)
