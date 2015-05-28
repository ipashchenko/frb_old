try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup


setup(
    name='frb',
    version='0.0',
    author='Ilya Pashchenko',
    author_email='in4pashchenko@gmail.com',
    packages=['frb', 'tests'],
    scripts=['bin/search_file.py'],
    url='https://github.com/ipashchenko/frb',
    license='LICENSE',
    description='Play with FRB search on RA data',
    long_description=open('README.rst').read(),
    install_requires=[
        "numpy >= 1.7.2"
    ],)
