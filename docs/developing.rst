==========
Developing
==========

Pre-requisites-  :ref:`Pre-requisites`

Virtual environments-  :ref:`Virtual environments`

Clone the libraries from GitHub- :ref:`Clone the libraries from GitHub`

Installing libraries via pip- :ref:`Installing libraries via pip`

Pre-requisites
==============

You will need

- Python
- Sphinx
- Virtual environments (optional)


Windows
'''''''
For Sphinx you will first need to have Python installed.
Python is availble from the Python `download page <PythonDownload>`_.
Once Python is installed Sphinx can be installed with pip.

.. code-block:: console

  pip install sphinx

It is also useful (but not a requirement to make use of virtual environments.
this allows you to have easy access to multiple versions of the library (e.g. working and development)

.. code-block:: console

  pip install virtualenv

OS X
''''

To install fun python things  you will need pip, if you don't already have
it you can use the `get-pip.py script <https://bootstrap.pypa.io/get-pip.py>`_.

.. code-block:: console

  curl https://bootstrap.pypa.io/get-pip.py > get-pip.py
  sudo -H python get-pip.py

Once pip is installed you can get Sphinx using

.. code-block:: console

    sudo -H pip install sphinx
    sudo -H pip install sphinx-fortran

It is also useful (but not a requirement to make use of virtual environments.
this allows you to have easy access to multiple versions of the library (e.g. working and development)

.. code-block:: console

    sudo -H pip install virtualenv

GNU/Linux
'''''''''
Use the command line to get all these things:

.. code-block:: console

    sudo apt-get install pythonX.Y-dev # Where X and Y are the major and minor version numbers of the Python you want to install, any version above 2.6 will work
    sudo apt-get install python-sphinx
    sudo apt-get install python-virtualenv


Virtual environments
====================

Possibly the best way to develop with this library is  through virtual environments.
This allows a library configuration to be matched to a  virtual environment that can be
easily moved between i.e. for a release version of the code you use to run models,
and a develop one which might have stuff you are doing that won't work for everyone

To create a virtual environement

Create a home for all virtual environments
''''''''''''''''''''''''''''''''''''''''''

The first task is to create a directory to hold the virtual environment installations

.. code-block:: console

  mkdir virtual_environments

This directory can be created anywhere.

Create a virtual environment
''''''''''''''''''''''''''''

The second task is to create a Python virtual environment to install the libraries

.. code-block:: console

  cd virtual_environments # change directory to where the virtual environment should be created
  virtualenv --system-site-packages nameofyourenvironment

The *--system-site-packages* flag allows the virtual environment to access all the packages that the system python has installed.
This is useful for big packages which may be required, for example; numpy or scipy.
The name of the virtual environment (in this case *nameofyourenvironment*) is determined by you, whatever makes sense, maybe
the name of the branch of these libraries you are on.

Activate virtual environment
''''''''''''''''''''''''''''

The third task is to activate the Python environment.

.. code-block:: console

  source /path/to/env/bin/activate

or in windows:

.. code-block:: console

  \path\to\env\Scripts\activate


The activate script may alter the command prompt to indicate the active virtual environment.
This script will also make changes to your path variables.
To undo these changes execute the *deactivate* script


.. code-block:: console

  deactivate


Clone the libraries from GitHub
===============================

`Set up a development triangle <https://gist.github.com/anjohnson/8994c95ab2a06f7d2339>`_.

If you are at University of Auckland, you can see some ABI documentation on development triangles here, including
a `step-by-step-guide <https://drive.google.com/open?id=1520jNTRMWXjvd47kS8UWOn02tsLZnppY>`_.

At the end of this process you should have a copy of this codebase local to your machine.
You can then install this local copy via pip, within a virtual environment, so you can update and run
as you go.


Installing libraries via pip
============================
You can now install the locally cloned package via pip, with the *-e* option enabling you to make changes
to the libraries and to use these changes in your virtual environment immediately

.. code-block:: console

    pip install -e /path/to/placentagen


