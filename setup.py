from setuptools import setup, find_packages

if __name__ == '__main__':
    setup(
        name='HiCEnterprise',
        version='0.1.2-1',
        author='Hania Kranas',
        packages=find_packages(exclude=['tests']),
        include_package_data=True,
        url='https://github.com/hansiu/HiCEnterprise',
        license='LICENSE.txt',
        description='Scripts for prediction of interaction between regions/domains based on Hi-C maps',
        long_description=open('README.md').read(),
        entry_points={'console_scripts': ['HiCEnterprise=HiCEnterprise.__main__:main']},
        keywords=['bioinformatics', 'hi-c', 'interactions', 'regions', 'domains', 'chromatin'],
        classifiers=[
            'Development Status :: 2 - Pre-Alpha'
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
        install_requires=[
            "numpy",
            "scipy",
            "statsmodels",
            "matplotlib"
        ],
        setup_requires=['pytest-runner'],
        tests_require=['pytest'],

    )
