from setuptools import setup

setup(
    name='dots_for_microarrays',
    version='0.1.1',
    author='Sandy Macdonald',
    author_email='sandyjmacdonald@gmail.com',
    packages=['dots_backend', 'dots_tests'],
    scripts=['dots_scripts/dots_workflow.py'],
    url='https://github.com/sandyjmacdonald/dots_for_microarrays',
    download_url='https://github.com/sandyjmacdonald/dots_for_microarrays/tarball/v0.1.1',
    license='LICENSE.txt',
    description='Simple analysis of Agilent one-color arrays.',
    keywords = ['bioinformatics', 'genomics', 'microarrays'],
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