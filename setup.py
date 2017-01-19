from setuptools import setup
from distutils.core import setup,Extension
from Cython.Build import cythonize
import os
import sys
setup(
   name="VaLoR",
   version="0.0.1",
   author="Xuefang Zhao, University of Michigan",
   author_email="xuefzhao@umich.edu",
   description="Long read based genomic structural variants validator.",
   packages=["valor_vali"],
   scripts=["valor_vali/valor"],
   package_data={
       "valor_vali": [
           "templates/pred_config",
           "templates/truth_config",
           "contrib/SMCScoring.py",
       ],
   },
    ext_modules = cythonize("valor_vali/*.pyx"),
     install_requires=[
     'cython', 'numpy','scipy','matplotlib','sklearn'
      ],
    license="Propriety",
    keywords="Long read",
    classifiers=[
       "Development Status :: 1 - Planning",
       "Environment :: Console",
       "Intended Audience :: Information Technology",
       "Intended Audience :: Science/Research",
       "License :: Other/Proprietary License",
       "Natural Language :: English",
       "Programming Language :: Python :: 2",
       "Programming Language :: Python :: 2 :: Only",
       "Topic :: Scientific/Engineering :: Bio-Informatics",
   ],
)

