.. image:: https://raw.githubusercontent.com/jlvdb/balltree/main/docs/source/_static/logo.png
    :width: 742
    :alt: balltree

|

.. image:: https://img.shields.io/pypi/v/balltree?logo=pypi&logoColor=blue
    :target: https://pypi.org/project/balltree/
.. image:: https://github.com/jlvdb/balltree/actions/workflows/python-ci.yml/badge.svg
    :target: https://github.com/jlvdb/balltree/actions/workflows/python-ci.yml
.. image:: https://readthedocs.org/projects/balltree/badge/?version=latest
    :target: https://balltree.readthedocs.io/en/latest/?badge=latest

|

A fast ball tree implementation for three dimensional (weighted) data with an
Euclidean distance norm. The base implementation is in `C` and there is a
wrapper for `Python`. The tree is optimised towards spatial correlation function
calculations since it provides fast counting routinge, e.g. by implementing a
dualtree query algorithm.

- Code: https://github.com/jlvdb/balltree.git
- Docs: https://balltree.readthedocs.io/
- PyPI: https://pypi.org/project/balltree/

.. toc

Installation
------------

A `C` library can be built with the provided make file, the python wrapper is
automatically compiled and installed with ``pip install balltree``.

The installation does not require any external `C` libraries, the python wrapper
requires the ``Python.h`` header (which should be included in a default python
installation) and `numpy` (including ``numpy/arrayobject.h``).


Maintainers
-----------

- Jan Luca van den Busch
  (*author*, Ruhr-Universit??t Bochum, Astronomisches Institut)
