===========
How to test
===========

Running tests-  :ref:`Running tests`

Creating tests-  :ref:`Creating tests`

Running tests
=============

The codebase on the repository is tested using travis continuous integration for python 2.7 and 3.0

To test locally on your machine, within the root directory you can run tests, first install nosetests:

.. code-block:: console

    >>> pip install nose

Then install coverage:

.. code-block:: console

    >>> pip install coverage

You can then check how your tests are doing using the following command:

.. code-block:: console

    >>> nosetests --with-coverage --cover-package=placentagen

The terminal output should look something like

.. code-block:: console

    ...............................
    Name                                 Stmts   Miss  Cover
    --------------------------------------------------------
    placentagen/__init__.py                  4      0   100%
    placentagen/analyse_tree.py            381    144    62%
    placentagen/generate_shapes.py         102      0   100%
    placentagen/grow_tree.py               673     55    92%
    placentagen/imports_and_exports.py     243    188    23%
    placentagen/pg_utilities.py             81      7    91%
    --------------------------------------------------------
    TOTAL                                 1484    394    73%
    ----------------------------------------------------------------------
    Ran 31 tests in 0.163s

    OK

If tests fail, and you are a developer, try to understand why they have failed and fix them! If tests fail and you have
never looked inside this python libraries, send any error messages to a developer.

We have continuous integration going, using Travis. So, if you want to contribute to this library, the tests have to pass!

If you don't want to use nosetests, an alternative is (from the root directory of the libraries):
    >>> python setup.py test


Creating tests
==============

We aim to have this code 100% covered by unit tests, and if you are developing we expect pull requests to include unit
tests, commented code, and documentation for new code. We try to practise what we preach, but we're almost exclusively
ex-spaghetti code developers so please bear with us if we don't have this right every time.

If you want to learn more about unit testing and test driven development, how about looking at
`this software carpentry course on testing <https://v4.software-carpentry.org/test/index.html>`_.

In brief, before you write your next bit of code, think first about how you will test it. What are the simplest steps you
could break your problem down into? and how would you make sure each of these steps work (for each simple step, is there something
you know the code should do?)? Then write your test, then write your code. Changes to the codebase should not (without good reason)
cause a test to fail.

In this code, each module (python file in `source/placentagen/modulename.py`) has a corresponding test in the `tests` directory (`tests/test_modulename.py`). We use the python
`unittest libraries <https://docs.python.org/2/library/unittest.html>`_.

Each testing file (`test_modulename.py`) starts with the following

.. code-block:: python

    from unittest import TestCase #Imports the unittest libraries we need

    import placentagen #Imports our libraries

    import anyotherlibraries #Imports any other libraries you need (keep this to a minimum - dont import libraries for the sake of it!)

We then set up `classes` of related tests. These are basically groups of tests that should be run together, and test the same or related functions.

.. code-block:: python

    class Test_group_of_tests(TestCase): #Group of tests

        def test_1(self): #The first thing we want to do here
            array_test = placentagen.function(function_inputs) #Function call
            self.assertTrue(array_test.all) #check something is true (in this case the function must output an array of trues and falses)

        def test_2(self): #Let's test something else
            epsilon = placentagen.function2(function_inputs) #function call
            self.assertTrue(epsilon == 0.1) #Something we know about the function output in a simple case

For any given python module, we will have a number of test classes, and test cases. The output of :ref:`Running tests` should
help you to understand how good a job you are doing on covering your modules with tests and whether they work. Talk to others about your
tests, because we don't want anything overcomplicated. Think about small simple tests and steps to follow, not your big problem (even
if you want to generate a million points in the end and check they are in a cylinder, this is not what your test should do, if you can test
your code with one point shouldn't it work on a million?).




