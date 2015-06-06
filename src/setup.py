#!/usr/bin/env python
#
# Setup script for the eqparse library
# $Id: setup.py
#
# Usage: python setup.py install
#

from os.path import abspath as abs_path
#from sys.path import append as sys_path
import sys
sys.path.append(abs_path('src'))

from distutils.core import setup

try:
    # add download_url syntax to distutils
    from distutils.dist import DistributionMetadata
    DistributionMetadata.classifiers = None
    DistributionMetadata.download_url = None
except:
    pass


DESCRIPTION="A mathematical system parser - from one generic master metadescription to user-defined language- or system- dependent formats"

LONG_DESCRIPTION="""\
The eqparse library is a library designed to manipulate a master 
metadescription file of a mathematical system  (which may consist 
of time-dependent equations, ordinary differential equations, 
parameters) into a particular file format (e.g. CPP, C, MATLAB, XPPAUT). 
This is useful for environments where multiple software applications/
languages are required for model development through prototyping to 
large scale model analysis (e.g. parameter and stability investigation, 
preparation of unit models to be organised into class objects for 
multi-unit model development etc). This is to (A) reduce the error 
in translating mathematical descriptions of code into different 
formats, (B) increase your time on on other research once parser has 
been implemented. Currently the input file for master metadescription 
is CSV (coma delimited) and example for Matlab, C++ and XPPAUT file 
output parsers are available as examples. EqParse can also 
automatically format any equations or dynamic variables (such as ODEs) 
into vectorised names (e.g. dydt[0], dydt[1], p[1], pstore[2]) 

Note: certain functions of the MATLAB parser prepares a file making 
use of another application created by the author called SYSSLI - 
a dynamic parameter slider for investigating mathematical systems in 
real time (MATLAB only)
"""
#url="personalpages.manchester.ac.uk/postgrad/haroon.arshad-2/",

setup(
    name="eqparse",
    version="1.0 beta",
    author="Haroon Arshad",
    author_email="haroon.arshad@student.manchester.ac.uk",
	url="https://github.com/h1arshad/EqParse",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    download_url="personalpages.manchester.ac.uk/postgrad/haroon.arshad-2/",
    license="Python (MIT style)",
    packages=["eqparse"],
    platforms="Python 1.5 and later.",
    classifiers=(
        "Development Status :: premature",
        "Operating System :: OS Independent",
        "Topic :: Text Processing :: Parser",
		"License :: OSI Approved :: MIT License"
        )
    )

# check install directory install serialisation file
# redfine this to include check where eqparse package is installed
#from eqparse import __eqp_dir, __eqp_memory_file, __eqp_memory_file_comb, __eqp_mac_address as mac_address
from eqparse import __eqp_dir, __eqp_memory_file, __eqp_memory_file_comb, __eqp_mac_address as mac_address
from pickle import dump as pdump
from pickle import load as pload

for file_name in [__eqp_memory_file, __eqp_memory_file_comb]:
	file = open(file_name,'wb')
	file.truncate()
	file.close()
	fprop = {};
	pdump(fprop,open(file_name,'wb'))

#file = open(__eqp_memory_file,'wb')
#file.truncate()
#file.close()
#fprop = {};
#pdump(fprop,open(__eqp_memory_file,'wb'))


#file = open(__eqp_memory_file_comb,'wb')
#file.truncate()
#fcomb = {}
#file.close()
#pdump(fcomb,open(__eqp_memory_file_comb,'wb'))

 



