
try:
    from setuptools import setup, find_packages
except ImportError:
    sys.exit("Please install setuptools.")


setup(
    name="ppanini",
    version="0.5.0",
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
        "Numpy >= 1.9.2",
        "Scikit-learn  >= 0.14.1",
        "Matplotlib >= 1.1.1",
        "Biopython >= 1.66"
    ],
    packages=find_packages(),
    #cmdclass={'install': Install},
    package_data={
        'ppanini' :[]},
    entry_points={
        'console_scripts': [
            'ppanini = ppanini.ppanini:_main',
            'prepannini = ppanini.utils.preppanini:main',
            'ppanini_eval_roc = ppanini.utils.plot_roc:evaluation_multi_roc'
        ]},
    test_suite= 'ppanini.tests.ppanini_test.main',
    zip_safe = False
 )

