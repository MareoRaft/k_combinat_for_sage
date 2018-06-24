#!/usr/bin/env python2.7
""" To build the docs, one should *not* run `make html` or `sphinx-build` directly, but rather run this file with `./build_docs.py`. """
from os import path
import subprocess

# get local path
PATH_REPO = path.dirname(path.dirname(path.abspath(__file__)))
PATH_DOCS = path.join(PATH_REPO, 'docs')

# get content to insert into code
PATH_INSERTION_CONTENT = path.join(PATH_DOCS, 'namespace_deletions.py')
with open(PATH_INSERTION_CONTENT, 'r') as file:
	insertion_content = file.read()

# get content of code
PATH_PACKAGE = path.join(PATH_REPO, 'k_combinat_for_sage')
PATH_CODE = path.join(PATH_PACKAGE, 'partition.py')
with open(PATH_CODE, 'r') as file:
	code = file.read()

# split code into pieces and recombine
pieces = code.split('# ^*^ sphinx insert ^*^')
(head, tail) = (pieces[0], pieces[1])
code_new = head + insertion_content + tail

# replace old code with new code
with open(PATH_CODE, 'w') as file:
	file.write(code_new)

# run sphinx to build the docs
cmd = ['sage', '-sh', '-c', 'make html']
subprocess.check_output(cmd)

# restore code to original state
with open(PATH_CODE, 'w') as file:
	file.write(code)
