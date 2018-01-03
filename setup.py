from setuptools import setup, find_packages
import re, io


if __name__ == '__main__':
    setup(
        name='HiCEnterprise',
        version='0.1.0',
        author='Hania Kranas',
        packages=["HiCEnterprise"],
        # url='', # TODO add github address
        license='LICENSE.txt',
        description='Scripts for prediction of interaction between regions/domains based on Hi-C maps',
        long_description=open('README.md').read(),
        entry_points={'console_scripts': ['HiCEnterprise=HiCEnterprise.__main__:main']},
        keywords=['bioinformatics','hi-c','interactions','regions','domains','chromatin'],
        # classifiers=[], # TODO add classifiers etc
         install_requires=[
             "numpy",
             "scipy",
             "statsmodels",
             "matplotlib"
         ],
        setup_requires=['pytest-runner'],
        tests_require=['pytest'],

    )
