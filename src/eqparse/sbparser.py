#! /usr/bin/env python

from baseparse import *
import time
"""
SBParser PARSER
"""

class SBParser (BaseParse):

	def __init__ (self, Lib):
		super(SBParser, self).__init__()
		self.l_enclose = "()"
		self.r_enclose = ")"
		self.operations = ['\nif(', ")\n", "else{\n\t", "}\nelse if{\n", "}\n", "!=", "((", ")^(", "))", "sqrt(("]

		self.vec_counter_start = 0
		self.comment_prefix = "**********"		
		self.fileExt = "txt"
		
		self.initialise_library (Lib)
		self.temp_name = ""
		self.temp_filename = ""
		
		#def set_directory(self, dir_str):
		#	""" Overridden XPPAUT Set directory where file specific parser files are saved to, 
		#	the directory can be a relative path and not necessarily an absolute path
		#	
		#	:param dir_str: string of path to directory
		#	"""
		#	if isinstance(dir_str, str) is False: error("set_directory argument not a string")
		#	elif not dir_str: pass
		#	elif (os.path.isdir(dir_str)) is False: error("directory does not exist")
		#	else: 
		#		self.directory = dir_str	
	
	def set_temp_ic_name (self) :
		#print "dir:" + self.directory
		#print "tit:" + self.title 
		#.replace('/','_').replace('\\','_') + self.title.replace('/','_').replace('\\','_')
		self.temp_name = "tdat_245345_" + self.directory.replace('/','_').replace('\\','_') + self.title.replace('/','_').replace('\\','_') + str(time.time())[-4:1]
		self.temp_filename = self.temp_name + ".py"
		
		#print self.directory
		#print self.title
		#print self.temp_name

		
	def printarray (self, c, specific = None):
		if specific is None:
			for elem in self.info['range']:
				print(self.names[elem] + ": " + str(c[elem]))
			print("\n")
		else:
			print(specific + ": " + str(c[self.names.index(specific)]))			

	def truncate_ics (self):
		self.set_temp_ic_name()

		#print "hello"
		absic_file = open(self.temp_filename,"w+")
		absic_file.truncate()
		
		#self.PERM_I_LIST = ['IF[[', "]]:", "ELSE:", "ELIF:", "END:", "!=", "POW[[", "||", "]]"]	### put in self. vars, helper

		#o_ops =  ["POW[[", "||", "]]","exp(","log(","log10("]
		#sub_ops = ["(", ")**(", ")","math.exp(","math.log(","math.log10("]
		o_ops =  [")^(","exp(","log(","log10(","sqrt("]
		sub_ops = [")**(","math.exp(","math.log(","math.log10(","math.sqrt("]
		repl_set = dict(zip(o_ops,sub_ops))
		abs_out = [x['init-value'] for x in self.data]
		mIndex = self.order_index(self.info['range'])
		#print mIndex
		absic_file.write("import math\n\n")
		absic_file.write("from math import sqrt as sqrt\n\n")
		absic_file.write("def get_value():\n")
		
		for key, val in repl_set.iteritems():
			for it in self.info['range']:
				abs_out[it] = abs_out[it].replace(key, val)
		#print len(mIndex)
		#JOJO=[]
		#print abs_out	
		for elem in mIndex:
			#print self.names[elem]
			#JOJO.append(self.names[elem])
			absic_file.write('\t' + self.names[elem] + ' = ' + abs_out[elem] + '\n')
			#absic_file.write('\tprint "' + self.names[elem] + '" + ": " + str(' + self.names[elem] + ')\n')
		
		#print list(set(self.names)-set(JOJO))
		#print set.(self.names).symmetric_difference(set(JOJO))

		absic_file.write("\treturn [")
		#for elem in mIndex[:-2]:
		for elem in self.info['range']:
			absic_file.write(self.names[elem] + ",")
		absic_file.write(self.names[-1])
			#print elem
		
		absic_file.write("]")
		absic_file.close()


		iModule = __import__(self.temp_name)
		new_absic = map(str, iModule.get_value())
		#print len(self.names)		
		#print len(new_absic)
		#print list(set(self.names)-set(JOJO))
		#print self.temp_name
		
	def parse_one (self):
		
		# data vars
		mfile = self.title + "/" + "sbmodel"
		mcomment = self.write_comment
		mwrite = self.write

		# open file
		self.open_file (mfile)
		
		# A : gather the function related index
		## note: j for i, is converting list-of-lists to list
		#def get_list (self, mKey, mLib = None, mPrefix ="", mPostfix ="", mIndex = None):
		#indexrhs = self.new_index(self.index['rhs']['yes']) # get index of non ode rhs
		#grhs = self.get_list ("rhs", self.data, "", "", indexrhs)
		#func_names = [j for i in [x['rhs'] for x in self.info['function']] for j in i] # ORIGINAL
		func_names = [self.names[elem] for elem in self.index['rhs']['yes']]
		#print func_names
		#print len(func_names)
		
		func_names = [self.names.index(x) for x in func_names]
		func_names.sort(key=int)
		is_func = {'yes':func_names,'no':[r for r in range(len(self.names)) if r not in func_names]}
		#print is_func
		# B : initial conditions get absolute numbers (need to run and create external file)
		self.truncate_ics()
		#startFile = "tdat_245345_" + self.title + ".py"
		#import startFile as semiSolve
		iModule = __import__(self.temp_name)
		#try:
		#	os.remove(self.temp_filename)
		#except OSError:
		#	pass
		#try:
		#	os.remove(self.temp_filename+"c")		
		#except OSError:
		#	pass
		new_absic = map(str, iModule.get_value())
		
		
		
		########  HEADER
		#mcomment (mfile, "XPPaut ode file for " + mfile)
		#mcomment (mfile, "H. Arshad 2013\n")
		#nnn = ['0']*len(arg_rhs)
	
		mcomment (mfile, "MODEL NAME","\n","")
		mwrite (mfile, self.title)
	
		#print len(new_absic)
		#print new_absic
		#print self.names
		#for i in range(len(self.names)):
		#	print str(i) + ": " + str(new_absic[i]) + " <-" + str(self.names[i])
		#	#print str(new_absic[i]) + ": " + str(self.names[i])
		#print new_absic[-1]
		
		mcomment (mfile, "MODEL NOTES","\n","\n")
		mwrite (mfile, "A note of this model will be here: currently not implemented")	

		mcomment (mfile, "MODEL STATES","\n","\n")
		self.pattern_write (mfile, ["name", "rhs"],["d/dt(",") = "," %"], self.data, self.index['ode']['yes'])
		mwrite (mfile,"")
		self.pattern_write (mfile, ["name", new_absic],["","(0) = ",""], self.data, self.index['ode']['yes'])

		mcomment (mfile, "MODEL PARAMETERS","\n","\n")
		#mIndex = self.new_index(self.index['var']['yes'], self.index['print']['no'], self.index['ode']['no'])
		#for i in mIndex:
		#	print self.names[i]
		self.pattern_write (mfile, ["name", new_absic],[""," = "," %parameters"], self.data, self.index['rhs']['no'])
		
		mcomment (mfile, "MODEL VARIABLES","\n","\n")
		#self.pattern_write (mfile, ["name", new_absic],[""," = "," %parameters"], self.data, self.index['rhs']['no'])


		mcomment (mfile, "MODEL REACTIONS","\n","\n")
		mIndex = self.new_index(self.index['rhs']['yes'], self.index['ode']['no'])
		self.pattern_write (mfile, ["name", "rhs"],[""," = "," %reaction"],self.data,mIndex)
		
		mcomment (mfile, "MODEL FUNCTIONS","\n","\n")


		mcomment (mfile, "MODEL EVENTS","\n","\n")


		mcomment (mfile, "MODEL MATLAB FUNCTIONS","\n","\n")


		self.close_file (mfile)