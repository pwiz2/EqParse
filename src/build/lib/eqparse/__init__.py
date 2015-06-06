# $Id: __init__.py
# --------------------------------------------------------------------
# The smc toolkit is
#
# Copyright (c) 2013-2015 by Haroon Arshad
#
# By obtaining, using, and/or copying this software and/or its
# associated documentation, you agree that you have read, understood,
# and will comply with the following terms and conditions:
#
# import eqparse
# Permission to use, copy, modify, and distribute this software and
# its associated documentation for any purpose and without fee is
# hereby granted, provided that the above copyright notice appears in
# all copies, and that both that copyright notice and this permission
# notice appear in supporting documentation, and that the name of
# Secret Labs AB or the author not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# SECRET LABS AB AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD
# TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANT-
# ABILITY AND FITNESS.  IN NO EVENT SHALL SECRET LABS AB OR THE AUTHOR
# BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY
# DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
# WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
# ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
# OF THIS SOFTWARE.
# --------------------------------------------------------------------

"""
eqparse import module initialise
"""

from site import getsitepackages
__eqp_dir = getsitepackages()[0] + "/eqparse/"
__eqp_memory_file =  __eqp_dir + 'eqpparseMemory.dat'
__eqp_memory_file_comb = __eqp_dir + 'eqpparseMemoryCombination.dat'
print(__eqp_dir)

# get mac address for identity of computer
from uuid import getnode as get_mac
__eqp_mac_address = get_mac()


__version__ = '0.1.0'
__license__ = 'MIT'
__author__ = 'Haroon Arshad'


import sys
import os
sys.path.append(os.path.dirname(__file__))

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))

#from .eqparse import Timer, CreateLibrary, BaseParse, CppParser, XppautParser, MatlabParser, SBParser, LatexParser

__all__ = [
	'Timer',
	'CreateLibrary',
	'BaseParse',
	'CppParser',
	'XppautParser',
	'MatlabParser',
	'SBParser',
	'LatexParser',
	'error'
]

import smc_helper_functions

import timer

import createlibrary

import baseparse

import cppparser

import xppautparser

import matlabparser

import sbparser

import latexparser



