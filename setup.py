from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='placentagen',
    version='0.1.0',
    packages=find_packages('source', exclude=['tests', 'tests.*', 'docs']),
    package_dir={'': 'source'},
    url='https://github.com/alysclark/placentagen.git',
    license=license,
    author='Alys Clark',
    author_email='alys.clark@auckland.ac.nz',
    description=''
)
