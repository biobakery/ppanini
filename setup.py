
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
VERSION = "0.7.1"
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
    version=VERSION,
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
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
        ],
    packages=setuptools.find_packages(),
    long_description="PPANINI (Prioritization and Prediction of functional Annotations" + \
    " for Novel and Important genes via automated data Network Integration) provides" + \
    " a computational pipeline to prioritize microbial genes based on their metagenomic" + \
    " properties (e.g. prevalence and abundance). The resulting prioritized list of gene" + \
    " candidates can then be analyzed further using our visualization tools." + \
    " PPANINI prioritizes important genes in a microbial community based on presenceabsence" + \
    " and abundance from metagenomic data. Sequencing a metagenome typically produces millions" + \
    " of short DNA/RNA reads. PPANINI takes a genes abundances table for all the samples in a study," + \
    " it ranks the important genes",
    install_requires=[  
        #Numpy >= 1.9.2",
        #"Scipy >= 0.13.0"
        #"Scikit-learn  >= 0.14.1",
        #"Matplotlib >= 2.0.2",
        #"Biopython >= 1.66",
        #"pandas >= 0.18.1"
    ],
    #cmdclass={'install': Install},
    include_package_data=True,
    package_data={
        #'ppanini' :['data/*'],
        'ppanini' :['tests/data/*','data/*'],
        },
    entry_points={
        'console_scripts': [
            'ppanini = ppanini.ppanini:_main',
            #'ppanini_preprocess = ppanini.tools.preppanini:main',
            'ppanini_gene_caller = ppanini.tools.ppanini_gene_caller:main',
            'ppanini_test = ppanini.tests.ppanini_test:main',
            'ppanini_scatterplot = ppanini.tools.ppanini_scatterplot:main',
            'ppanini_barplot = ppanini.tools.ppanini_barplot:main',
            'ppanini_rocplot = ppanini.tools.ppanini_rocplot:main',
            'ppanini_rev_uniref_mapper = ppanini.tools.attach_GO:rev_load_polymap',
            'ppanini_join_tables = ppanini.tools.ppanini_join_tables:main',
            'ppanini_rename_contigs = ppanini.tools.ppanini_rename_contigs:main',
            'ppanini_infer_gene = ppanini.tools.ppanini_infer_gene:main',
            'ppanini_cluster2genes = ppanini.tools.ppanini_cluster2genes:main',
            'ppanini_press = ppanini.tools.ppanini_press:main',
            'ppanini_fasta_select = ppanini.tools.fasta_select:main'            
        ]},
    test_suite= 'ppanini_test.get_unittests()',
    zip_safe = False
 )

