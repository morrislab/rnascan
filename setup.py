from setuptools import find_packages
from distutils.core import setup, Extension

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
        Extension('motif_scan.BioAddons.motifs._pwm',
                  ["motif_scan/BioAddons/motifs/_pwm.c"],
                  include_dirs=[numpy_include_dir],
))

setup(name='motif_scan',
      version='0.7.1',
      description='Scan RBP motifs and secondary structure from PFMs',
      url='http://github.com/kcha/motif_scan',
      author='Kevin Ha',
      author_email='k.ha@mail.utoronto.ca',
      license='MIT',
      packages=find_packages(),
      scripts=['bin/motif_scan', 'bin/combine_pfms'],
      install_requires=['setuptools',
                        'pandas >= 0.17',
                        'numpy >= 1.10.0',
                        'biopython >= 1.66'],
      ext_modules=EXTENSIONS,
      zip_safe=False
      )
