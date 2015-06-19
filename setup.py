from setuptools import setup

setup(
    name='dots_for_microarrays',
    version='0.1.0',
    author='Sandy J. Macdonald',
    author_email='sandyjmacdonald@gmail.com',
    packages=['dots_backend', 'dots_tests'],
    scripts=['dots_scripts/dots_workflow.py'],
    url='https://github.com/sandyjmacdonald/dots_for_microarrays',
    license='LICENSE',
    description='Simple analysis of Agilent one-color arrays.',
    long_description=open('README').read(),
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