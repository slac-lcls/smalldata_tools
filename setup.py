
"""
Setup script for pscache.
"""

from glob import glob


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


setup(name='smalldata_tools',
      version='0.1',
      author="Silke Nelson",
      author_email="snelson@stanford.edu",
      description='tools for creating and analyzing LCLS small data',
      packages=["smalldata_tools"],
      package_dir={"smalldata_tools": "smalldata_tools"},
      scripts=[s for s in glob('scripts/*') if not s.endswith('__.py')],
      test_suite="test")

