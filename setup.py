import os.path
import sys
from setuptools import find_packages
from distutils.core import setup, Extension

if sys.version_info < (2, 7):
      sys.stderr.write("rnascan requires Python 2.7, or Python 3.5 or later. "
                       "Python %d.%d detected.\n" % sys.version_info[:2])
      sys.exit(1)
elif sys.version_info[0] == 3 and sys.version_info[:2] < (3, 5):
      sys.stderr.write("rnascan requires Python 2.7, or Python 3.5 or later. "
                       "Python %d.%d detected.\n" % sys.version_info[:2])
      sys.exit(1)

# Borrowing setup.py code from Biopython

def is_pypy():
    import platform
    try:
        if platform.python_implementation() == 'PyPy':
            return True
    except AttributeError:
        # New in Python 2.6, not in Jython yet either
        pass
    return False


def can_import(module_name):
    """can_import(module_name) -> module or None"""
    try:
        return __import__(module_name)
    except ImportError:
        return None


def is_Numpy_installed():
    if is_pypy():
        return False
    return bool(can_import("numpy"))

EXTENSIONS = []

if is_Numpy_installed():
    import numpy
    numpy_include_dir = numpy.get_include()
    EXTENSIONS.append(
        Extension('rnascan.BioAddons.motifs._pwm',
                  ["rnascan/BioAddons/motifs/_pwm.c"],
                  include_dirs=[numpy_include_dir],
))

here = os.path.abspath(os.path.dirname(__file__))
exec(open(os.path.join(here, 'rnascan/version.py')).read())


setup(name='rnascan',
      version=__version__,
      description='Scan RBP motifs and secondary structure from SSMs',
      url='http://github.com/morrislab/rnascan',
      author='Kevin Ha, Kate Cook, Kaitlin Laverty',
      author_email='k.ha@mail.utoronto.ca, kate.cook@gmail.com, kaitlin.laverty@mail.utoronto.ca',
      license='AGPLv3',
      packages=find_packages(),
      scripts=['scripts/run_folding'],
      install_requires=['setuptools',
                        'pandas >= 0.24',
                        'numpy >= 1.10.0',
                        'biopython >= 1.66'],
      entry_points={
            'console_scripts': [
                'rnascan = rnascan.rnascan:main'
            ]
      },
      ext_modules=EXTENSIONS,
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose']
      )
