.. _developer_manual:

Building the developer manual
*****************************

The developer manual is not available online. Rather, it must be built locally with the tool `FORD (FORtran Documenter) <https://github.com/Fortran-FOSS-Programmers/ford>`__ as follows\:

* Installing FORD with `pip <https://pypi.org/project/pip/>`__: ``pip install ford``.
* Creating/updating the manual: ``ford ford.md`` (in the main directory of the SICOPOLIS installation).
* Output is in HTML -> ``docs/ford/index.html``.

If required, see the `FORD documentation <https://forddocs.readthedocs.io>`__ for further information.
