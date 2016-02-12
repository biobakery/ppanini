
try:
    from setuptools import setup, find_packages
except ImportError:
    sys.exit("Please install setuptools.")


setup(
    name="ppanini",
    version="0.6.1",
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
    long_description=open('readme.md').read(),
    install_requires=[  
        #"Numpy >= 1.9.2",
        #"Scikit-learn  >= 0.14.1",
        #"Matplotlib >= 1.1.1",
        #"Biopython >= 1.66"
    ],
    packages=find_packages(),
    #cmdclass={'install': Install},
    package_data={
        'ppanini' :[]},
    entry_points={
        'console_scripts': [
            'ppanini = ppanini.ppanini:_main',
            'preppanini = ppanini.utils.preppanini:main',
            'ppanini_eval_roc = ppanini.utils.plot_roc:main',
            'ppanini_visualizer = ppanini.utils.ppanini_visualizer:main',
            'ppanini_plot_metagenome_genome = ppanini.utils.plot_metagenome_genome:main',
            'ppanini_plot_genome_hits = ppanini.utils.plot_genome_hits:main',
            'ppanini_normalize_table = ppanini.utils.normalize_table:main',
            'ppanini_join_tables = ppanini.utils.join_tables:main',
            'ppanini_imp_centroids_prabXtract = ppanini.utils.imp_centroids_prabXtract:main',
            'ppanini_imp_centroids_extracter = ppanini.utils.imp_centroids_extracter:main',
            'ppanini_create_mapper = ppanini.utils.create_mapper:main',
            'ppanini_write_mapper = ppanini.utils.write_mapper:main'
        ]},
    test_suite= 'ppanini.tests.ppanini_test.main',
    zip_safe = False
 )

