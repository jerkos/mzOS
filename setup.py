from distutils.core import setup

setup(
    name='mzOS',
    version='0.1',
    packages=['mzos', 'mzos.tests', 'mzos.ressources'],
    packages_data=['mzos.tests.data'],
    url='http://github.com/jerkos/mzOS',
    license='',
    author='Marco',
    author_email='marc.dubois@omics-services.com',
    description='', requires=['bioservices', 'scipy', 'numpy', 'pandas']
)
