Requested Features for the tutorial
------------------------------------

  * On the web
  * pretty
  * mini `code quotes`
  * large ```code blocks```
  * inline $\\latex$
  * runs code (sage kernel)



It seems that an ipynb will do all this except for the sage kernel, which AFAIK is only available on Cocalc.

I think $$ works for latex insertions in nbs.

I think

~~~
this
~~~

is one way to 'fence' code blocks in a ipynb.


notebook metadata example:
 "metadata": {
  "kernelspec": {
   "display_name": "Ruby 2.2.1",
   "language": "ruby",
   "name": "ruby"
  },
  "language_info": {
   "file_extension": ".rb",
   "mimetype": "application/x-ruby",
   "name": "ruby",
   "version": "2.2.1"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 0

example of a nb w latex and code quotes:
https://nbviewer.jupyter.org/url/jdj.mit.edu/~stevenj/IJulia%20Preview.ipynb

GITHUB will automatically make your notebook viewable online (https://blog.github.com/2015-05-07-github-jupyter-notebooks-3/).  ALL YOU NEED TO DO is put a ipynb file into your repo.

The following PROGRAM (https://github.com/aaren/notedown) allows you to write your IPython Notebooks in markdown!
md to ipynb:
notedown input.md > output.ipynb
ipynb to md:
notedown input.ipynb --to markdown --strip > output.md


You can run an ipynb with any installed kernel (including sage) with a command such as:
sage -n jupyter test_sage_nb.ipynb

We MAY be able to use "%load_ext sage" to get sage to work in PYTHON KERNEL notebooks, which means we could then use GITHUB!
https://stackoverflow.com/questions/23384070/taking-advantage-of-sage-and-ipython-notebook-in-the-same-page-or-rather-combi

