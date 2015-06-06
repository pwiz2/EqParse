#! /usr/bin/env python

# initialise all

#from os.path import realpath as smc_realpath
#from os.path import join as smc_join
#from os.path import join as smc_split
#from sys import path as smc_path
#import inspect
#if smc_realpath(os.path.abspath(smc_join(smc_split(inspect.getfile( inspect.currentframe() ))[0],"smc"))) not in smc_path:
#    smc_path.insert(0, cmd_subfolder)

#from sys.path import append as append_path
#from os.path import abspath as abs_path
#append_path(abs_path('smc'))


from sys import exit as do_exit
"""
	Error reporting
	@arg err_str - string to output once error is instantiated from program
"""
def error (err_str):
	if isinstance (err_str,str) == False: error ("Error message not set properly")
	do_exit("ERROR: " + err_str + "\n")


# python 2.x to python 3.x compatibility	
try:
	basestring
except NameError:
	basestring = str
