from setuptools import setup

setup(name='motif_scan',
      version='0.7.0',
      description='Scan RBP motifs and secondary structure from PFMs',
      url='http://github.com/kcha/motif_scan',
      author='Kevin Ha',
      author_email='k.ha@mail.utoronto.ca',
      license='MIT',
      packages=['motif_scan'],
      scripts=['bin/motif_scan', 'bin/combine_pfms'],
      install_requires=['setuptools',
                        'pandas >= 0.17',
                        'numpy >= 1.10.0',
                        'biopython >= 1.66'],
      zip_safe=False
      )
