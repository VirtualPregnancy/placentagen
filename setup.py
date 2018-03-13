from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='genericplacentas',
    version='0.1.0',
    packages=find_packages(exclude=('tests', 'docs')),
    url='https://github.com/alysclark/generic_placentas.git',
    license=license,
    author='Alys Clark',
    author_email='alys.clark@auckland.ac.nz',
    entry_points={
        'console_scripts': ['genericplacentas=generic_placentas.application:main']
    },
    description=''
)
