from setuptools import setup, find_packages

setup(
    name='mzOS',
    version='0.1.1',
    packages=find_packages(),
    include_package_data=True,
    url='http://github.com/jerkos/mzOS',
    license='MIT',
    author='Marco',
    author_email='cram@hotmail.fr',
    description='Heuristic based feature annotations/identifications of LC-MS metabolomics dataset.',
    long_description=open('README.md').read(),
    requires=['bioservices', 'scipy', 'numpy', 'six', 'sklearn'],
    classifiers=['Development Status :: 3 - Alpha',
                 'License :: OSI Approved :: MIT License',
                 'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3.3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 'Topic :: Scientific/Engineering :: Bio-Informatics'
                 ]
)
