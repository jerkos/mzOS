from setuptools import setup, find_packages


entry_points = {
    "console_scripts": [
        'mzos = mzos.scripts.mzos_script:main'
    ]
}

setup(
    name='mzOS',
    version='0.1.3',
    packages=find_packages(),
    package_data={'mzos.ressources': ['*'],
                  'mzos.third_party.emass': ['*']
                  },
    url='http://github.com/jerkos/mzOS',
    license='MIT',
    author='Marco Dubois',
    author_email='cram@hotmail.fr',
    description='Heuristic based feature annotations/identifications of LC-MS metabolomics dataset.',
    entry_points=entry_points,
    long_description=open('README.md').read(),
    requires=['bioservices', 'scipy', 'numpy', 'six', 'sklearn', 'pandas'],
    classifiers=['Development Status :: 3 - Alpha',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3.3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 'Topic :: Scientific/Engineering :: Bio-Informatics'
                 ]
)
