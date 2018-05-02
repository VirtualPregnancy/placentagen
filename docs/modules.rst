=======
Modules
=======

Modules should be collections of functions that are linked to one another in a logical way
make sure that before you add some code to a module, it actually fits in that module.

Its pretty simple to make a new module, just make a file in the `source/placentagen` directory called
`module_name.py`, add a corresponding `test_module_name.py` to the `tests` directory, and then add some
documentation (list on this page, and create a `.rst` file to the `docs/Modules` directory that auto-creates
documentation from doc strings in the module itself.

List of modules
---------------

.. toctree::
   :maxdepth: 1

   Modules/analyse_tree
   Modules/generate_shapes
   Modules/grow_tree
   Modules/imports_and_exports