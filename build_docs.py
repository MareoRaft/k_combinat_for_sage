#!/usr/bin/env python2.7
from __future__ import print_function
from os import path
import subprocess
import sys
import time

path_repo = path.dirname(path.dirname(path.abspath(__file__)))
sys.path.append(path_repo)
from config import PATH

start_time = time.time()
print('started building docs...')

# get content to insert into code
with open(PATH['namespace_deletions'], 'r') as file:
	insertion_content = file.read()

FILE_NAMES = [
	'core.py',
	'partition.py',
	'k_shape.py',
	'skew_partition.py',
	'root_ideal.py',
	'strong_marked_tableau.py',
	'all.py',
	'shorthands.py',
]
file_name_to_code = dict()
print('modifying code files...')
for file_name in FILE_NAMES:
	PATH_CODE = path.join(PATH['package'], file_name)
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
print('running sphinx...')
cmd = ['sage', '-sh', '-c', 'make html']
subprocess.check_output(cmd)

# restore code to original state
print('reverting code files...')
for file_name in FILE_NAMES:
	PATH_CODE = path.join(PATH['package'], file_name)
	if path.exists(PATH_CODE):
		with open(PATH_CODE, 'w') as file:
			file.write(file_name_to_code[file_name])

# all done!
print('finished building docs.', end='')
end_time = time.time()
elapsed_time = end_time - start_time
print(' Elapsed time = {}'.format(elapsed_time))
