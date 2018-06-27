#!/usr/bin/env python3
# To build a new VERSION of this library for PyPI, do the following:
from os import path
import subprocess

# get paths
PATH_REPO = path.dirname(path.dirname(path.abspath(__file__)))
PATH_DOCS = path.join(PATH_REPO, 'docs')
PATH_VERSION = path.join(PATH_REPO, 'VERSION')
PATH_SETUP = path.join(PATH_REPO, 'setup.py')
PATH_DIST_ALL = path.join(PATH_REPO, 'dist/*')

# Warning message to the script user
proceed = input('This script will (among other things) perform a git commit, pull, and push.  You should be on the master branch and confident.  Do you wish to proceed?\n').strip()
if proceed != 'yes':
	raise SystemExit('User aborted version builder.  No action was taken.')

# Change version number in `VERSION` file.
version_old = open(PATH_VERSION).read().strip()
version_new = input('The current version number is {}.  What should the new version number be?\n'.format(version_old)).strip()
with open(PATH_VERSION, 'w') as file:
	file.write(version_new)

# Commit all changes and push to master.
print('git committing and pushing to GitHub...')
cmd = ['git', 'add', PATH_VERSION]
subprocess.run(cmd, check=True)
cmd = ['git', 'commit', '-m', version_new]
subprocess.run(cmd, check=True)
cmd = ['git', 'pull']
subprocess.run(cmd, check=True)
cmd = ['git', 'push']
subprocess.run(cmd, check=True)

# git tag with version number and push to github
cmd = ['git', 'tag', version_new]
subprocess.run(cmd, check=True)
cmd = ['git', 'push', '--tags', 'origin', 'master']
subprocess.run(cmd, check=True)

# create a distribution for PyPI and upload to PyPI
# detailed instructions here: https://packaging.python.org/tutorials/packaging-projects/
print('uploading version to PyPI...')
cmd = ['python', PATH_SETUP, 'sdist']
subprocess.run(cmd, check=True)
cmd = ['twine', 'upload', PATH_DIST_ALL]
subprocess.run(cmd, check=True)
print('version building script completed.')
