#!/usr/bin/env python3
# To build a new VERSION of this library for PyPI, do the following:
import os
from os import path
import subprocess

# get paths
PATH_CWD = os.getcwd()
PATH_REPO = path.dirname(path.abspath(__file__))
PATH_VERSION = path.join(PATH_REPO, 'VERSION')
PATH_SETUP = path.join(PATH_REPO, 'setup.py')
PATH_DIST = path.join(PATH_REPO, 'dist')

# helpers
# input = raw_input # if using python2
def run(cmd_list):
	assert isinstance(cmd_list, list)
	# use 'check_output' instead of 'run' if using python2
	subprocess.run(cmd_list, check=True)

# Warning message to the script user
if PATH_CWD != PATH_REPO:
	raise SystemExit('This script MUST be executed from the root of the repository.\nrepo root: {}\nyour cwd : {}'.format(PATH_REPO, PATH_CWD))
proceed = input('This script will (among other things) perform a git commit, pull, and push.  You should be on the master branch and confident.  Do you wish to proceed?\n').strip()
if proceed != 'yes':
	raise SystemExit('User aborted version builder.  No action was taken.')

# Change version number in `VERSION` file.
version_old = open(PATH_VERSION).read().strip()
version_new = input('The current version number is {}.  What should the new version number be?\n'.format(version_old)).strip()
if version_new == version_old:
	raise SystemExit('That is the same as the previous version number.  You must increment the version number.')
with open(PATH_VERSION, 'w') as file:
	file.write(version_new)

# Commit all changes and push to master.
print('git committing and pushing to GitHub...')
run(['git', 'add', PATH_VERSION])
run(['git', 'commit', '-m', version_new])
run(['git', 'pull'])
run(['git', 'push'])

# git tag with version number and push to github
run(['git', 'tag', version_new])
run(['git', 'push', '--tags', 'origin', 'master'])

# create a distribution for PyPI and upload to PyPI
# detailed instructions here: https://packaging.python.org/tutorials/packaging-projects/
print('uploading version to PyPI...')
run(['python', PATH_SETUP, 'sdist'])
PATH_DIST_VERSION = path.join(PATH_DIST, 'k_combinat_for_sage-{}.tar.gz'.format(version_new))
run(['twine', 'upload', PATH_DIST_VERSION])
print('Version building script completed successfully!')
