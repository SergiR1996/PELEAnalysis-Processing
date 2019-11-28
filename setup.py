from setuptools import find_packages, setup


setup(name='PELEAnalysis-Processing',
      packages=find_packages(),
      version='1.0.0',
      license='MIT',
      description='Tools for MD & PELE analysis, preparation of the system, mutation of PDB files using Schrodinger python API, tools to study esterases, and more',
      author='Sergi Rodà',
      author_email='sergi.rodallordes@bsc.es',
      url='https://github.com/SergiR1996/PELEAnalysis-Processing',
      install_requires=['argparse', 'math', 'matplotlib', 'numpy', 'sklearn', 'pandas', 'mdtraj', 'seaborn', 'scipy', 'pickle', 'mlxtend', 'itertools', 'multiprocessing', 'time'],
     ) 
