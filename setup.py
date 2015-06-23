from setuptools import setup

setup(
    name='dots_for_microarrays',
    version='0.1.2',
    author='Sandy Macdonald',
    author_email='sandyjmacdonald@gmail.com',
    packages=['dots_backend', 'dots_tests'],
    scripts=['dots_scripts/dots_workflow.py'],
    url='https://github.com/sandyjmacdonald/dots_for_microarrays',
    download_url='https://github.com/sandyjmacdonald/dots_for_microarrays/tarball/v0.1.2',
    license='MIT',
    description='Simple analysis of Agilent one-color arrays.',
    keywords = ['bioinformatics', 'genomics', 'microarrays'],
    classifiers = [
    'License :: OSI Approved :: MIT License',
    'Operating System :: Unix',
    'Operating System :: POSIX :: Linux',
    'Operating System :: MacOS :: MacOS X',
    'Programming Language :: Python :: 2.7',
    'Topic :: Software Development :: Libraries :: Python Modules',
    'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'statsmodels',
        'scikit-learn',
        'bokeh',
        'pillow',
        'brewer2mpl',
        'nose'
    ],
    setup_requires=[
        'nose',
        'numpy',
    ]
)