
try:
    import setuptools
    #from setuptools import setup, find_packages
except ImportError:
    sys.exit("Please install setuptools.")
import os
import urllib
from setuptools.command.install import install as _install

# try to download the bitbucket counter file to count downloads
# this has been added since PyPI has turned off the download stats
# this will be removed when PyPI Warehouse is production as it
# will have download stats
COUNTER_URL="http://bitbucket.org/biobakery/ppanini/downloads/counter.txt"
counter_file="counter.txt"
if not os.path.isfile(counter_file):
    print("Downloading counter file to track ppanini downloads"+
        " since the global PyPI download stats are currently turned off.")
    try:
        file, headers = urllib.urlretrieve(COUNTER_URL,counter_file)
    except EnvironmentError:
        print("Unable to download counter")

setuptools.setup(
    name="ppanini",
    version="0.6.2",
    license="MIT",
    description="PPANINI: Prioritization and Prediction of functional Annotations for Novel and Important genes via automated data Network Integration.",
    author="Gholamali Rahnavard, Afrah Shafquat, Eric A. Franzosa, and Curtis Huttenhower",
    author_email="rahnavar@hsph.harvard.edu",
    maintainer="Gholamali  Rahnavard",
    maintainer_email="gholamali.rahnavard@gmail.com",
    url="http://huttenhower.sph.harvard.edu/ppanini",
    keywords=["important", "gene", "prioritization", "annotation" ],
    classifiers=[
        "Programming Language :: Python",
        "Development Status :: 1 - Development",
        "Environment :: Console",
        "License :: MIT License",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    packages=setuptools.find_packages(),
    long_description=open('readme.md').read(),
    install_requires=[  
        #"Numpy >= 1.9.2",
        #"Scikit-learn  >= 0.14.1",
        #"Matplotlib >= 1.1.1",
        #"Biopython >= 1.66"
    ],
    #cmdclass={'install': Install},
    package_data={
        'ppanini' :['tests/data/*']},
    entry_points={
        'console_scripts': [
            'ppanini = ppanini.ppanini:_main',
            'preppanini = ppanini.utils.preppanini:main',
            'ppanini_eval_roc = ppanini.utils.plot_roc:main',
            'ppanini_test = ppanini.tests.ppanini_test:main'
        ]},
    test_suite= 'ppanini.tests.ppanini_test.get_unittests()',
    zip_safe = False
 )

