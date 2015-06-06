#! /usr/bin/env python

from baseparse import *
import time
"""
Xppaut PARSER
"""

class LatexParser (BaseParse):

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
		
	def parse_one (self, concateq_unordered = [], parameterCompat_change = []):
		
		# data vars
		mfile = self.title
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
		try:
			os.remove(self.temp_filename)
		except OSError:
			pass
		try:
			os.remove(self.temp_filename+"c")		
		except OSError:
			pass
		new_absic = map(str, iModule.get_value())
		#print new_absic
		#for elem in mIndex:
		#	print self.names[elem] + " + " + str(new_absic[elem])
		##print self.directory
		##print self.title
		##print self.temp_name	
		
		# C : gather new NAME and RHS for ode/funcs arguments based on aternative names
		# this is because XPPAUT don't allow certain names due to the implementation
		# in order to save editing one by one, a user is allowed to enter a dictionary
		# of parameter/function/ode name changes by passing parameterCompat_change []
		# note the following might be slow since [x['names'] for x in self.data] creates a new list. Is this necessary?
		if parameterCompat_change:
			pnew_orignames = parameterCompat_change.keys()
			pnew_newnames = parameterCompat_change.values()
		else:
			pnew_orignames = []
			pnew_newnames = []
		#print pnew_orignames
		#print pnew_newnames

		modNames = self.search_and_replace(self.get_list('name'), pnew_newnames, pnew_orignames)
		modRhs = self.search_and_replace(self.get_list('rhs'), pnew_newnames, pnew_orignames)
		#print modRhs
		
		# order concateq_unordered according to self.order.init_value
		concateq = []
		for elem in self.index['order']['init-value']:
			#-- note self.names in if loop since concateq_unordered is user defined with original values
			if self.names[elem] in concateq_unordered:
				concateq.append(modNames[elem])
		#print concateq
		
		#-- get the index of the element in original names list order (since it is not the whole list)
		concat_index = [modNames.index(elem) for elem in concateq]
		#a = self.get_list('rhs') #[x['rhs'] for x in self.data]
		#-- may need to change below - repeating the same function above
		#a = self.search_and_replace(self.get_list('rhs'), pnew_newnames, pnew_orignames)
		#-- change this to array form [] list concatenate
		#for elem in concat_index:
		#	a[elem] = ''
		#print len(a)
		#for elem in self.info['range']:
		#	if elem in concat_index:
		#		a[elem] = ''
		#	else:
		#		a[elem] = modRhs[elem]
		a = ['' if elem in concat_index else modRhs[elem] for elem in self.info['range']]
		#print len(a)
			
		# get the rhs of the VARS you want to embed into the equation
		#concateq_rhs_sub = []
		#for elem in concateq:
		#	concateq_rhs_sub.append('(' + self.data[self.names.index(elem)]['rhs'] + ')')			
		concateq_rhs_sub = ['(' + modRhs[elem] + ')' for elem in concat_index]
		#for i in range(len(concateq_rhs_sub)):
		#	print str(concateq_rhs_sub[i]) + "\n" + str(concateq_rhs_new[i]) + "\n" + str(mnew[i]) + "\n\n" 
		
		#print concateq_rhs_sub
		#print concateq[concateq.index("I_NaKCl_Na")]
		#print concateq_rhs_sub[concateq.index("I_NaKCl_Na")]

		#concateq_rhs_sub = [''.join(line) for line in [ [concateq_rhs_sub[concateq.index(word)] if word in concateq else word for word in mline] for mline in [re.split('(\W)',line) for line in concateq_rhs_sub]]]
		#print concateq_rhs_sub[concateq.index("I_NaKCl_Na")]
		
		#import time
		#start = time.time()

		#### alternative search and replace (slow but maintains changes)
		for elem in range(len(concateq_rhs_sub)):	
			line = re.split('(\W)',concateq_rhs_sub[elem])
			#print 'line: ' + str(line)
			for word_pos in range(len(line)):
				 if line[word_pos] == concateq[elem]:
				 	pass
				 elif line[word_pos] in concateq:
				 	#print ':::::' + line[word_pos] + ' : ' + concateq_rhs_sub[concateq.index(line[word_pos])]
				 	line[word_pos] = concateq_rhs_sub[concateq.index(line[word_pos])]
				 	#print '>>>>' + line[word_pos]
			concateq_rhs_sub[elem] = ''.join(line)

		#concateq_rhs_sub = [''.join(line) for line in [ [self.data[self.names.index(word)]['rhs'] if word in concateq else word for word in mline] for mline in [re.split('(\W)',line) for line in concateq_rhs_sub]]]
		##print self.data[self.names.index("I_NaKCl_Cl")]['rhs']
		#self.printarray(a,"I_Natotm")
		#print concateq_rhs_sub[concateq.index("I_NaKCl_Na")]
		
		"""get rhs of ODE: if argument was supplied, the ODE"""
		a = [''.join(line) for line in [[concateq_rhs_sub[concateq.index(word)] if word in concateq else word for word in mline] for mline in [re.split('(\W)',line) for line in a]]]
		a = [re.split('(\W)',line) for line in a]
		#print a
		#repl_set = dict(zip(concateq,concateq_rhs_sub))
		#for key, val in repl_set.iteritems():	# replace concat equations embed into others
		#	for it in self.info['range']:
		#		a[it] = a[it].replace(key, val)
				
		#ode_list = self.new_names (self.names, self.index['ode']['no'],["",""],True)
		""" get each variable-dependent rhs and check for ode dependencies """
		#ode_list = self.get_only_names (self.names, self.index['ode']['yes'])
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
				##names_ode_dep.append(self.names[iter])
		#ntemp = [self.names[elem] for elem in self.index['ode']['yes']]
		#print names_ode_dep
		#print arg_list

		""" get all variables which are dependent on ode variables (indirect) """		
		arg_list_indirect = [list(set(names_ode_dep).intersection(elem)) for elem in a]
		#print arg_list_indirect	
		#self.printarray(arg_list_indirect,"I_Na_leak")
		#self.printarray(arg_list,"I_Na_leak")
		
		##arg_list_indirect = [[arg_list[self.names.index(elem)] for elem in row] if row is not None else [] for row in arg_list_indirect]
		arg_list_indirect = [[arg_list[modNames.index(elem)] for elem in row] if row is not None else [] for row in arg_list_indirect]

		#print arg_list_indirect
		
		#self.printarray(arg_list_indirect,"I_Na_leak")
		#print arg_list_indirect
		
		arg_list_indirect = [list(set([item for sublist in row for item in sublist])) for row in arg_list_indirect]
		#self.printarray(arg_list_indirect,"I_Na_leak")
		#self.printarray(arg_list,"I_Na_leak")

		
		arg_list = [arg_list[elem] + arg_list_indirect[elem] for elem in self.info['range'] ]
		arg_list = [list(set(row)) for row in arg_list]
		
		#print "\n"
		#print arg_list 
		#print arg_list_indirect[self.names.index("I_BKCa")]
		
		#print arg_list_indirect[self.names.index("I_BKCa")]

		# check if other ode dep index is a part of it
		#index_ode_dep = [self.info['range'][arg_list.index(x)] for x in arg_list if x]
		#names_ode_dep = [self.names[elem] for elem in index_ode_dep]
		#names_ode_dep = [self.info['range'][arg_list.index(x)] for x in arg_list if x]
		#names_ode_dep = []
		#for iter in range(len(arg_list)):
		#	if arg_list[iter]:
		#		names_ode_dep.append(self.names[iter])
		#rhs_arg_list = [list(set(names_ode_dep).intersection(elem)) for elem in a]
		#print rhs_arg_list

		

		#rhs_arg_list = [arg_list[self.names.index(row)] for row in [elem for elem in rhs_arg_list] if row]
		#rhs_arg_list = [[arg_list[self.names.index(elem)] for elem in row] if row is not None else [] for row in rhs_arg_list]
		#rhs_arg_list = [set([item for sublist in row for item in sublist]) for row in rhs_arg_list]

		#for elem in rhs_arg_list
		#print rhs_arg_list
				
		
		
		arg_list = [','.join(elem) for elem in arg_list]
		arg_list = ['('+elem+')' if elem else elem for elem in arg_list]
		#arg_names = []
		#for elem in self.info['range']:
		#	if elem in self.index['ode']['no']:
		#		arg_names.append(modNames[elem]+arg_list[elem])
		#		##arg_names.append(self.names[elem]+arg_list[elem])
		#	else:
		#		arg_names.append(modNames[elem])
		#		##arg_names.append(self.names[elem])
		arg_names = [modNames[elem]+arg_list[elem] if elem in self.index['ode']['no'] else modNames[elem] for elem in self.info['range']]
		#for i in range(len(arg_names_new)):
		#	print str(arg_names[i]) + ":::" + str(arg_names_new[i])
		
				
		#arg_rhs = self.search_and_replace ([x['rhs'] for x in self.data], arg_names)
		a = [''.join(line) for line in a]
		arg_rhs = self.search_and_replace (a, arg_names)

		#print arg_list
		#print "END OF arg_list"		
		#print arg_rhs
		#print "END OF arg_rhs"
		#print new_absic
		#print "END OF new_absic"
		#print arg_names
		#print "END OF arg_names"
		
		#for elem in range(len(self.names)):
		#	print str(self.names[elem]) + ' = ' + str(arg_names[elem])
		
		""" ARG1 concatenate certain formulas into equation """
		#### alternative search and replace (slow but maintains changes)
		for elem in range(len(arg_rhs)):	
			line = re.split('(\W)',arg_rhs[elem])
			#print 'line: ' + str(line)
			for word_pos in range(len(line)):
				 #if line[word_pos] == concateq[elem]:
				 #	pass
				 if line[word_pos] in concateq:
				 	#print ':::::' + line[word_pos] + ' : ' + concateq_rhs_sub[concateq.index(line[word_pos])]
				 	line[word_pos] = concateq_rhs_sub[concateq.index(line[word_pos])]
				 	#print '>>>>' + line[word_pos]
			arg_rhs[elem] = ''.join(line)
			
		
		#self.printarray(a,"I_Na_leak")
		#self.printarray(arg_list,"I_Na_leak")
		#self.printarray(arg_names,"I_Na_leak")
		#print self.temp_filename
		#print arg_list
		#print arg_names
		
		#print "part 1:  %.2f" % (time.time()-start), ' seconds.'
	
		
		#start = time.time()
		
		################### NEW CHANGE DEPENDENCIES

		#a = [re.split('(\W)',line) for line in a]
		#a get_list(mKey, mLib=None, mPrefix='', mPostfix='', mIndex=None)[source]
		
		#""" get index of rhs in function order """
		#narg_index = []
		#for i in range(len(self.info['function'])):
		#	narg_index = narg_index + self.info['function'][i]['index']
		#
		#""" get rhs (narg_rhs) and original names (narg_names) for the index of function order (narg_index)
		#in addition get the ode names list (narg_ode_names)"""	
		#narg_rhs = self.get_list('rhs',self.data,'','',narg_index)
		#narg_names = self.get_only_names(self.names, narg_index)
#
#		
		#narg_ode_names = self.get_only_names(self.names, self.index['ode']['yes'])
		#narg_ode_rhs = self.get_list('rhs',self.data,'','',self.index['ode']['yes']) # none function order	
#
		#narg_nonode_names = self.get_only_names(self.names, self.new_index(self.index['ode']['no'],self.index['rhs']['yes']))				
		#narg_nonode_rhs = self.get_list('rhs',self.data,'','',self.new_index(self.index['ode']['no'],self.index['rhs']['yes'])) # none function order
		#
#
		#
		#""" get index dependecy order for non-ode functions """		"""    """
		#index_nonode = self.new_dependancy_index(narg_nonode_rhs, self.get_only_names(self.names, self.new_index(self.index['ode']['no'],self.index['rhs']['yes'])))
		#index_nonode = self.new_dependancy_index(self.new_index(self.index['ode']['no'],self.index['rhs']['yes']),'rhs')
#
		#
		#""" set new order for non-ode as names and rhs """
		#narg_nonode_rhs = [narg_nonode_rhs[elem] for elem in index_nonode]
		#narg_nonode_names = [narg_nonode_names[elem] for elem in index_nonode]

		#""" split version of RHS to compare vaiables"""
		#rhs_split = [re.split('(\W)',line) for line in narg_nonode_rhs]
		#rhs_split_ode = [re.split('(\W)',line) for line in narg_ode_rhs]
#
		#		
		#""" get if each has any ode dependence (direct) """
		#DVAR = [set(narg_ode_names).intersection(elem) for elem in rhs_split] 
#
		#""" Check which functions are dependent on each other (by iter value) """
		#""" compare values and return index: find a better way, prob need to use itertools (?) """
		#indices = []
		#names_dict = dict((k,i) for i,k in enumerate(narg_nonode_names)) # dict for iter
		#for line in rhs_split:
		#	#print names_dict.values()
		#	dn = set(names_dict.keys()).intersection(line) # dn = dependent names
		#	indices.append([names_dict[elem] for elem in dn])
		#
		#"""combine each dependency with [can be written in list comprehension]"""
		#for it in range(0,len(indices)):
		#	for collec in indices[it]:
		#		DVAR[it] = set(DVAR[collec]).union(DVAR[it])
##
		#""" Do the same for each ODE, by finding out which of the ODEs ar dependent on each other"""
		#rhs_split_ode = [re.split('(\W)',line) for line in narg_ode_rhs]
		#for it in range(0,len(narg_ode_names)):
		#	func_dep_odes = set(narg_names).intersection(rhs_split_ode[it])
		#	#narg_names.in
		#		
		#"""organise new names depending on variable (replace function with argument)"""
		#narg_nonode_names_replaced = []
		#DVAR = ['('+','.join(s)+')' if s else '' for s in list(DVAR)]
		#for it in range(0,len(DVAR)):
		#	narg_nonode_names_replaced.append(narg_nonode_names[it] + DVAR[it])
		#
		#""" Now replace rhs with new names (need better performing one)""" 
		#for line_iter in range(0,len(rhs_split)):
		#	for word_iter in range(0,len(rhs_split[line_iter])):
		#		for wname_iter in range(0,len(narg_nonode_names)):
		#			if rhs_split[line_iter][word_iter] == narg_nonode_names[wname_iter]:
		#				rhs_split[line_iter][word_iter] = narg_nonode_names_replaced[wname_iter]
		#names_nonode_rhs = [''.join(line) for line in rhs_split]
		#
		#for line_iter in range(0,len(rhs_split_ode)):
		#	for word_iter in range(0,len(rhs_split_ode[line_iter])):
		#		for wname_iter in range(0,len(narg_nonode_names)):
		#			if rhs_split_ode[line_iter][word_iter] == narg_nonode_names[wname_iter]:
		#				rhs_split_ode[line_iter][word_iter] = narg_nonode_names_replaced[wname_iter]
		#names_ode_rhs = [''.join(line) for line in rhs_split_ode]
		#print names_ode_rhs


		#print "part 2:  %.2f" % (time.time()-start), ' seconds.'
		
		#narg_nonode_names = [narg_nonode_names[iter] + DVAR[iter] for iter in range(0,len(narg_nonode_names))]
		
		#for i in range(0,len(narg_index)):
		#	print str(narg_index[i]) + " : " + str(narg_names[i]) + " = " + str(narg_rhs[i])
		#print index_nonode
		#for i in range(0,len(narg_nonode_names)):
		#	print str(i) + " : " + str(narg_nonode_names[i]) + " = " + str(names_nonode_rhs[i])
		#for i in range(0,len(narg_nonode_names)):
		#	print str(narg_nonode_names[i]) + " = " + str(DVAR[i])
		#print DVAR
		################### END OF CHANGE DEPENDENCIES
		
		
		
		
		
		########  HEADER
		mcomment (mfile, "XPPaut ode file for " + mfile)
		mcomment (mfile, "H. Arshad 2013\n")
		nnn = ['0']*len(arg_rhs)
	
		# INITIAL CONDITIONS OF SYSTEM OF ODES	
		mcomment (mfile, ["#","# Initial conditions","#"], "\n", "\n\n")
		#mIndex =  self.order_index(self.index['ode']['yes'])
		mIndex = self.order_index(self.index['ode']['yes'])
		#print "HELLO::" + str(len(new_absic))
		#print "HELLO::" + str(len(mIndex))
		#print mIndex
		#print self.title
		#print "new_absic size :" + str(len(new_absic))
		#print "self.names size :" + str(len(self.names))
		#print self.temp_filename
		#print len(new_absic)
		#print len(self.names)
		##self.pattern_write (mfile, ["name", new_absic], ["", "(0)=", ""], self.data, mIndex)
		self.pattern_write (mfile, [modNames, new_absic], ["", "(0)=", ""], self.data, mIndex)

		# DEFINE PARAMS
		#print self.info['function']
		#print self.index['ode']['yes']
		mcomment (mfile, ["#","# Define parameter values","#"], "\n", "\n\n")
		mIndex =  self.order_index(self.index['ode']['no'],is_func['no'])
		##self.pattern_write (mfile, ["name", new_absic], [ "param ", "=", ""], self.data, mIndex)
		self.pattern_write (mfile, [modNames, new_absic], [ "param ", "=", ""], self.data, mIndex)
	
		# DEFINE THE FUNCTIONS IN MODEL	
		mcomment (mfile, ["#","# Functions used in the model","#"], "\n", "\n\n")
		#print set(self.info['range']).difference(concat_index)
		mIndex = self.order_index(self.index['ode']['no'],self.new_index(is_func['yes'],list(set(self.info['range']).difference(concat_index))))
		#concat_index = sorted(concat_index,reverse=True)
		#print concat_index
		#for elem in sorted(concat_index,reverse=True):
		#	mIndex.pop(elem)
		#print mIndex
		self.pattern_write (mfile, [arg_names, arg_rhs], [ "", "=", ""], self.data, mIndex)
		#self.pattern_write (mfile, [arg_names, nnn], [ "", "=", ""], self.data, mIndex)
			
		# ODE DECLARATIONS (differential equations
		mcomment (mfile, ["#","# ODE declarations/differential equations","#"], "\n", "\n\n")	
		mIndex =  self.order_index(self.index['ode']['yes'])
		self.pattern_write (mfile, [arg_names, arg_rhs], [ "d", "/dt = ", ""], self.data, mIndex)
		#self.pattern_write (mfile, [arg_names, nnn], [ "d", "/dt = ", ""], self.data, mIndex)

		
		# define properties of XPPaut system
		mcomment (mfile, self.xpp_settings,"\n","\n\n")
		#mcomment (mfile, ["@total=200000","@dt=50","@method=stiff","@xlo=0","@xhi=180000","@ylo=0","@yhi=0.001","@bounds=1000","@nout=1","@maxstor=50000","@bell=0","@back=White","@autovar=c","@ntst=150","@nmax=400","@npr=100","@dsmin=0.000001","@dsmax=0.5","@ds=0.01","@parmin=0","@parmax=50","@epsu=0.000001","@epss=0.000001","@autoxmin=0","@autoxmax=10","@autoymin=0","@autoymax=1"])

		# footer and close file
		mwrite (mfile, "\n\ndone")
		self.close_file (mfile)