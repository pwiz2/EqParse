#! /usr/bin/env python

from baseparse import *
import time
"""
Xppaut PARSER
"""

class XppautParser (BaseParse):

	def __init__ (self, Lib):
		super(XppautParser, self).__init__()
		self.l_enclose = "()"
		self.r_enclose = ")"
		self.operations = ['\nif(', ")\n", "else{\n\t", "}\nelse if{\n", "}\n", "!=", "((", ")^(", "))", "sqrt(("]

		self.vec_counter_start = 0
		self.comment_prefix = "# "		
		self.fileExt = "ode"
		
		self.initialise_library (Lib)
		self.temp_name = ""
		self.temp_filename = ""
		
		# default xpp settings (array of strings - each array is new line
		# note: this is temporary, will need to better redefine these
		self.xpp_settings = ["@ total=200000,dt=50,method=stiff,xlo=0,xhi=180000,ylo=0,yhi=0.001,bounds=1000,nout=1,maxstor=50000,bell=0,back=White",
"@autovar=c,ntst=150,nmax=400,npr=100,dsmin=0.000001,dsmax=0.5,ds=0.01,parmin=0,parmax=50,epsu=0.000001,epss=0.000001",
"@ autoxmin=0,autoxmax=10,autoymin=0,autoymax=1"]

	
	def set_temp_ic_name (self) :
		self.temp_name = "tdat_245345_" + self.directory.replace('/','_').replace('\\','_') + self.title.replace('/','_').replace('\\','_') + str(time.time())[-4:1]
		self.temp_filename = self.temp_name + ".py"
	    
	def xpp_search_and_define (self) :
		pass
		
	def set_xpp_settings (self, input) :
		self.xpp_settings = input;
		
		
	def printarray (self, c, specific = None):
		if specific is None:
			for elem in self.info['range']:
				print(self.names[elem] + ": " + str(c[elem]))
			print("\n")
		else:
			print(specific + ": " + str(c[self.names.index(specific)]))			

	def truncate_ics (self):
	
		self.set_temp_ic_name()

		absic_file = open(self.temp_filename,"w+")
		absic_file.truncate()

		o_ops =  [")^(","exp(","log(","log10(","sqrt("]
		sub_ops = [")**(","math.exp(","math.log(","math.log10(","math.sqrt("]
		repl_set = dict(zip(o_ops,sub_ops))
		abs_out = [x['init-value'] for x in self.data]
		mIndex = self.order_index(self.info['range'])
		absic_file.write("import math\n"+
		                  "from math import sqrt as sqrt\n"+
		                  "def get_value():\n")
		
		for key, val in repl_set.items():
			for it in self.info['range']:
				abs_out[it] = abs_out[it].replace(key, val)

		for elem in mIndex:
			absic_file.write('\t' + self.names[elem] + ' = ' + abs_out[elem] + '\n')

		absic_file.write("\treturn [")
		for elem in self.info['range']:
			absic_file.write(self.names[elem] + ",")
		absic_file.write(self.names[-1])
		
		absic_file.write("]")
		absic_file.close()


		iModule = __import__(self.temp_name)
		new_absic = map(str, iModule.get_value())
		
	def parse_one (self, concateq_unordered = [], parameterCompat_change = []):
		
		mfile = self.title
		mcomment = self.write_comment
		mwrite = self.write

		# open file
		self.open_file (mfile)
		
		# A : gather the function related index
		func_names = [self.names[elem] for elem in self.index['rhs']['yes']]
		
		func_names = [self.names.index(x) for x in func_names]
		func_names.sort(key=int)
		is_func = {'yes':func_names,'no':[r for r in range(len(self.names)) if r not in func_names]}

		# B : initial conditions get absolute numbers (need to run and create external file)
		self.truncate_ics()
		iModule = __import__(self.temp_name)
		try:
			os.remove(self.temp_filename)
		except OSError:
			pass
		try:
			os.remove(self.temp_filename+"c")		
		except OSError:
			pass
		new_absic = map(str, iModule.get_value())

		
		# C : gather new NAME and RHS for ode/funcs arguments based on aternative names
		# this is because XPPAUT don't allow certain names due to the implementation
		# in order to save editing one by one, a user is allowed to enter a dictionary
		# of parameter/function/ode name changes by passing parameterCompat_change []
		# note the following might be slow since [x['names'] for x in self.data] creates a new list. Is this necessary?
		if parameterCompat_change:
			pnew_orignames = list(parameterCompat_change.keys()) # need to change this? ppp.py
			pnew_newnames = list(parameterCompat_change.values())
		else:
			pnew_orignames = []
			pnew_newnames = []

		modNames = self.search_and_replace(self.get_list('name'), pnew_newnames, pnew_orignames)
		modRhs = self.search_and_replace(self.get_list('rhs'), pnew_newnames, pnew_orignames)
		#print modRhs
		
		# order concateq_unordered according to self.order.init_value
		concateq = []
		for elem in self.index['order']['init-value']:
			#-- note self.names in if loop since concateq_unordered is user defined with original values
			if self.names[elem] in concateq_unordered:
				concateq.append(modNames[elem])
		
		#-- get the index of the element in original names list order (since it is not the whole list)
		concat_index = [modNames.index(elem) for elem in concateq]

		a = ['' if elem in concat_index else modRhs[elem] for elem in self.info['range']]
			
		concateq_rhs_sub = ['(' + modRhs[elem] + ')' for elem in concat_index]



		#### alternative search and replace (slow but maintains changes)
		for elem in range(len(concateq_rhs_sub)):	
			line = re.split('(\W)',concateq_rhs_sub[elem])
			for word_pos in range(len(line)):
				 if line[word_pos] == concateq[elem]:
				 	pass
				 elif line[word_pos] in concateq:
				 	line[word_pos] = concateq_rhs_sub[concateq.index(line[word_pos])]
			concateq_rhs_sub[elem] = ''.join(line)

		
		"""get rhs of ODE: if argument was supplied, the ODE"""
		a = [''.join(line) for line in [[concateq_rhs_sub[concateq.index(word)] if word in concateq else word for word in mline] for mline in [re.split('(\W)',line) for line in a]]]
		a = [re.split('(\W)',line) for line in a]


		""" get each variable-dependent rhs and check for ode dependencies """
		ode_list = self.get_only_names (modNames, self.index['ode']['yes'])
		arg_list = [list(set(ode_list).intersection(elem)) for elem in a]
		#-- which is faster? set version or not? which one?
		#set_arg_list = [set(ode_list).intersection(elem) for elem in a]
		#print arg_list
		#print set_arg_list
		#print arg_list[self.names.index("I_Na_NSCCb")]
		
		""" get all function variables which are dependent on ODE variables (at first glance) """
		#-- note: the following may need to be changed for speed. Looks ugly too!
		# see SOF: "remove empty strings from a list of strings python"
		names_ode_dep = []
		for iter in range(len(arg_list)):
			if arg_list[iter]:
				names_ode_dep.append(modNames[iter])


		""" get all variables which are dependent on ode variables (indirect) """		
		arg_list_indirect = [list(set(names_ode_dep).intersection(elem)) for elem in a]

		arg_list_indirect = [[arg_list[modNames.index(elem)] for elem in row] if row is not None else [] for row in arg_list_indirect]
	
		arg_list_indirect = [list(set([item for sublist in row for item in sublist])) for row in arg_list_indirect]
		
		arg_list = [arg_list[elem] + arg_list_indirect[elem] for elem in self.info['range'] ]
		arg_list = [list(set(row)) for row in arg_list]	
		
		arg_list = [','.join(elem) for elem in arg_list]
		arg_list = ['('+elem+')' if elem else elem for elem in arg_list]

		arg_names = [modNames[elem]+arg_list[elem] if elem in self.index['ode']['no'] else modNames[elem] for elem in self.info['range']]
		
				
		#arg_rhs = self.search_and_replace ([x['rhs'] for x in self.data], arg_names)
		a = [''.join(line) for line in a]
		arg_rhs = self.search_and_replace (a, arg_names)


		""" ARG1 concatenate certain formulas into equation """
		#### alternative search and replace (slow but maintains changes)
		for elem in range(len(arg_rhs)):	
			line = re.split('(\W)',arg_rhs[elem])
			for word_pos in range(len(line)):
				 if line[word_pos] in concateq:
				 	line[word_pos] = concateq_rhs_sub[concateq.index(line[word_pos])]
			arg_rhs[elem] = ''.join(line)
			

		
		
		
		
		
		########  HEADER
		mcomment (mfile, "XPPaut ode file for " + mfile)
		mcomment (mfile, "H. Arshad 2013\n")
		nnn = ['0']*len(arg_rhs)
	
		# INITIAL CONDITIONS OF SYSTEM OF ODES	
		mcomment (mfile, ["#","# Initial conditions","#"], "\n", "\n\n")
		mIndex = self.order_index(self.index['ode']['yes'])

		self.pattern_write (mfile, [modNames, list(new_absic)], ["", "(0)=", ""], self.data, mIndex)

		mcomment (mfile, ["#","# Define parameter values","#"], "\n", "\n\n")
		mIndex =  self.order_index(self.index['ode']['no'],is_func['no'])
		self.pattern_write (mfile, [modNames, list(new_absic)], [ "param ", "=", ""], self.data, mIndex)
	
		# DEFINE THE FUNCTIONS IN MODEL	
		mcomment (mfile, ["#","# Functions used in the model","#"], "\n", "\n\n")
		mIndex = self.order_index(self.index['ode']['no'],self.new_index(is_func['yes'],list(set(self.info['range']).difference(concat_index))))
		self.pattern_write (mfile, [arg_names, arg_rhs], [ "", "=", ""], self.data, mIndex)
			
		# ODE DECLARATIONS (differential equations
		mcomment (mfile, ["#","# ODE declarations/differential equations","#"], "\n", "\n\n")	
		mIndex =  self.order_index(self.index['ode']['yes'])
		self.pattern_write (mfile, [arg_names, arg_rhs], [ "d", "/dt = ", ""], self.data, mIndex)

		
		# define properties of XPPaut system
		mcomment (mfile, self.xpp_settings,"\n","\n\n")

		# footer and close file
		mwrite (mfile, "\n\ndone")
		self.close_file (mfile)