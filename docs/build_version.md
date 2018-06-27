To build a new VERSION of this library for PyPI, do the following:

  1. Change version number in `VERSION` file.
  2. Commit all changes and push to master.
  3. git tag <version number>
  4. git push --tags origin master
  5. python setup.py sdist
  6. twine upload dist/*

detailed instructions here: https://packaging.python.org/tutorials/packaging-projects/
