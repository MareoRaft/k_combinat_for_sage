#!/usr/bin/env python2.7
from os import path
import subprocess

# get paths
PATH_REPO = path.dirname(path.dirname(path.abspath(__file__)))
PATH_DOCS = path.join(PATH_REPO, 'docs')
PATH_PACKAGE = path.join(PATH_REPO, 'k_combinat_for_sage')

# get content to insert into code
PATH_INSERTION_CONTENT = path.join(PATH_DOCS, 'namespace_deletions.py')
with open(PATH_INSERTION_CONTENT, 'r') as file:
	insertion_content = file.read()

FILE_NAMES = [
	'core.py',
	'partition.py',
	'k_shape.py',
	'skew_partition.py',
	'root_ideal.py',
	'all.py',
	'shorthands.py',
]
file_name_to_code = dict()
for file_name in FILE_NAMES:
	PATH_CODE = path.join(PATH_PACKAGE, file_name)
	if path.exists(PATH_CODE):
		# get content of code
		with open(PATH_CODE, 'r') as file:
			file_name_to_code[file_name] = file.read()

		# split code into pieces and recombine
		pieces = file_name_to_code[file_name].split('# ^*^ sphinx insert ^*^')
		(head, tail) = (pieces[0], pieces[1])
		code_new = head + insertion_content + tail

		# replace old code with new code
		with open(PATH_CODE, 'w') as file:
			file.write(code_new)

# run sphinx to build the docs
cmd = ['sage', '-sh', '-c', 'make html']
subprocess.check_output(cmd)

# restore code to original state
for file_name in FILE_NAMES:
	PATH_CODE = path.join(PATH_PACKAGE, file_name)
	if path.exists(PATH_CODE):
		with open(PATH_CODE, 'w') as file:
			file.write(file_name_to_code[file_name])
