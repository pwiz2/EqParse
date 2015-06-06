#! /usr/bin/env python

import csv
import re
import copy
import errno

from sys import platform as _platform
import os
from smc_helper_functions import error




class BaseParse (object):
	""" 
	Base generic parse class containing attributes and functions necessary to
	derive new parser types
	
	Description: The main objective of :class:`BaseParse` is to organise universal
	functions and member variables that are inherited in other class 
	parser modules. Current supplied parsers include the :class:`MatlabParser`,
	:class:`CppParser` and :class:`XPPautParser` as child modules.
	
	Attributes:
	
	"""
	
	
	def __init__ (self):
		""" Initialise universal variables usually important in most or all
		language syntax of the following. if a module inherits this class, 
		the following variables usually need to be redefined according to its
		language syntax. The following parameters are automatically defined as
		__init__ is called
		
		:param l_enclose: The left character to access a piece of memory is usually a left square bracket (C++) or a  normal left bracket (MATLAB
		:param r_enclose: Similar to :param l_enclose:
		:param comment_prefix: character used for single-line comment
		:param vec_counter_start: starting number to access first element of a container/array (default: 0)
		:param container_type: a list of container types used
		:param PERM_I_LIST: the input file contain preconditioned statements in csv input file such as if statements or power function.
		:param operations: the required result during parse of the :param PERM_I_LIST: list to the respective language syntax
		"""
		
		
		self.l_enclose = ""
		
		self.r_enclose = ""
		
		self.comment_prefix = ""
		
		self.vec_counter_start = 0
		
		self.var_prefix = ""			# DECLARE_VAR_PREFIX
		
		self.container_type = ""		# CPP only []
		
		self.operations = []			# COND_LIST
		
		self.PERM_I_LIST = ['IF[[', "]]:", "ELSE:", "ELIF:", "END:", "!=", "POW[[", "||", "]]", "SQRT[["]	### put in self. vars, helper
		
		self.file_ext = ""
		
		self.directory = ""
		
		self.open_files = {}
		#pass
		#self.lib = CreateLibrary("asmc_mod.csv")	# share actual data
		
		# toold for pickle library
		#self
		self.DOPARSE = 1;
			
			
			
	def initialise_library (self, lib): #, k = True):
		"""In order to initialise library
		
		:param lib: The created library from the :class:`createlibrary` class
		:returns: void
		"""

		self.__storageOfLibraryPriv = lib
		self.data_vec = self.__storageOfLibraryPriv.copy_lib_vec()
		self.data = self.__storageOfLibraryPriv.copy_lib()
		self.info = self.__storageOfLibraryPriv.info
		self.index = self.__storageOfLibraryPriv.info['index']
		self.names = self.__storageOfLibraryPriv.names
		self.size_vars = self.__storageOfLibraryPriv.size_vars
		self.title = self.info['title']
		self.null = self.__storageOfLibraryPriv.null
		self.directory = self.__storageOfLibraryPriv.directory

		self.DOPARSE = self.info['doWeParse'];
					
		self.replace_operations()
		# vectorize data
		vnames = self.vectorise_name(
			self.new_index(self.index['var']['yes'],self.index['print']['yes']),
			"p_store")
		
		vnames = self.vectorise_name(
			self.new_index(self.index['var']['yes'],self.index['print']['no']), 
			"p",
			vnames)
		
		vnames = self.vectorise_name (
			self.index['ode']['yes'], 
			"dydt", 
			vnames)
		self.search_and_replace (self.data_vec, vnames)


	def get_directory (self):
		return self.directory
		
		
	def open_file (self, file_n_str_key):
		if self.directory:
			self.open_files[file_n_str_key] = open (self.directory + "/" + file_n_str_key + "." + self.file_ext, "w")
		else:	
			self.open_files[file_n_str_key] = open (file_n_str_key + "." + self.file_ext, "w")
		if not self.open_files[file_n_str_key]: error("no files created in parser class " + self.__class__.__name__)

	
	def close_file (self, file_n_str):
		if self.open_files[file_n_str]:
			self.open_files[file_n_str].close()
		else: error("File does not exist")
	
	def end (self):
		# check if all files are closed
		for of in self.open_files:
			if of: error("file " + of.name + " is still open")
	
	def replace_operations(self, ops = None):
		if ops is None:
			ops = self.operations
		repl_set = dict(zip(self.PERM_I_LIST,ops))
		for elem in ['rhs','init-value']:
			for key, val in repl_set.items(): # python 2.x: .iteritems():
				for it in self.info['range']:
					self.data[it][elem] = self.data[it][elem].replace(key, val)
					self.data_vec[it][elem] = self.data_vec[it][elem].replace(key, val)


	
	# @[universal] fast version
	def search_and_replace (self, origContainer, repNames, origNames = None, return_it = True):
		"""Search and replace the variables of a container of equations (normally the 
		RHS or initial condition) to specified new names. All containers
		should have same length and coincide to variable  name properties
		
		:param origContainer: data of formular/equations containing variables one wants to change
		:param repNames: the names one wants to replace to parallel to :param `origNames`:
		:param origNames: the original names to replace (default - original variable names)
		"""		
		if origNames is None:
			origNames = self.names	
			
			
		if isinstance (origContainer[0], str):
			return  [''.join(line) for line in [[repNames[origNames.index(word)] if word in origNames else word for word in mline] for mline in [re.split('(\W)',line) for line in origContainer]]]

		elif isinstance (origContainer[0], dict):	# return_it not available here
			if origContainer is self.data_vec:
				# here we mutate actual code
				stored_keys = ['name', 'init-value', 'rhs']		
					
				for value in self.info['range']:
					for key in stored_keys:	
						origContainer[value][key] = ''.join([repNames[origNames.index(elem)] if elem in self.names else elem for elem in re.split('(\W)',origContainer[value][key])])
			
			elif origContainer is self.data:
				error("You can't search and replace the original container (forbidden)")
			
			else:
				error("error?? why cant we have other dictionaries then")



	# @[universal] get_list, with index if required
	# can be used to get simple prefixed/postfixes new names for original names ()
	# @[arg] mKey: key to find
	# @[arg] mKey: key to find
	def get_list (self, mKey, mLib = None, mPrefix ="", mPostfix ="", mIndex = None):
		"""get a list of a particular property of the variables. The list of properties depends on a key which associates a particular property of a variables library (e.g. comment, rhs). 
		
		:param mKey: the property type of a variable
		:param mLib: library to use e.g self.data (normal format) self.data_vec (vectorised format) (default - self.data)
		:param mPrefix: pre-fix the returned property with a string.
		:param mPostfix: post-fix the returned property with a string.
		"""		
		if mKey not in self.info['keys']: error("key does not exist, for function get_list")
		if mLib is None: mLib = self.data
		#if not mPrefix and if not mPostfix:
		#	return 
		if mIndex is None: #mIndex = self.info['range']
			#if isinstance (mPrefix, str) and isinstance (mPostfix, str):
			new_l=[mPrefix + item[mKey] + mPostfix for item in mLib]
			#else: 
			
			return new_l
		else:
			#new_l = []
			#for line in mLib:
			#for which in mIndex:
				#print line #[mKey]
			#	new_l.append(mPrefix + mLib[which][mKey] + mPostfix)
			new_l=[mPrefix + mLib[which][mKey] + mPostfix for which in mIndex]
			return [mPrefix + mLib[which][mKey] + mPostfix for which in mIndex]

	
	def get_list_vectorised (self, mKey, mPrefix ="", mPostfix ="", mIndex = None):
		if mKey not in self.info['keys']: error("key does not exist, for function get_list")
		if mIndex is None: #mIndex = self.info['range']
			new_l=[mPrefix + item[mKey] + mPostfix for item in self.data_vec]			
			return new_l
		else:
			new_l=[mPrefix + mLib[which][mKey] + mPostfix for which in mIndex]
			return [mPrefix + mLib[which][mKey] + mPostfix for which in mIndex]


	# @[universal] pattern print: print selection of variables with same pattern
	# @[arg] w_file : the file to write too (can be a file object or string of file id key)
	# @[arg] w_data_tuple : a tuple(?) of lists in order of print, can be string (a key of data)
	#		or a new data set
	# @[arg] w_str_tuple : same string to print on either side of data
	# @[arg] w_lib : ((optional)) library of variables to use (e.g. vectorised version or otherwise)
	#		note: have to define if w_data_tuple is list of strings
	# @[arg] indexes : ((optional)) index of values to print from_vars
	# note: by definition w_str_tuple has one more element than w_data_tuple to
	#       fill the data in 
	def pattern_write (self, w_file, w_data_tuple, w_str_tuple, w_lib = None, mIndex = None, end_wrr_ = "\n", beg_wrr_ = ""):
		"""Write to file output each of the required variables with a similar coded pattern
		
		:param str w_file: file to write to specified by string used in open_file function
		:param list w_data_tuple: combination of strings (if same for each variable) and lists (each with same size) that will be printed in consecutive order
		:param str w_str_tuple: 
		:param w_lib: post-fix the returned property with a string.
		:param list mIndex: specify index of variables to go through.		
		"""		
		# perform checks
		fname = " pattern_write "
		if isinstance(w_data_tuple, list) and isinstance(w_str_tuple, list):
			if len(w_data_tuple) != len(w_str_tuple) - 1:
				error("not right sizes in function " + fname)
			if len(w_data_tuple) == 1:	# if only element is one list
				w_data_tuple = [w_data_tuple]
			#if len(w_str_tuple) == 2 and isinstance(w_data_tuple,str):
			#	w_data_tuple = [w_data_tuple]
		elif isinstance(w_data_tuple, str) and isinstance(w_str_tuple, list) and len(w_str_tuple) == 2:
			w_data_tuple = [w_data_tuple]
		else: error("not correct types in function " + fname)

		if mIndex is None:
			if w_lib is None: mIndex = range(len(w_data_tuple[0]))
			else: mIndex = self.info['range']		
		# convert and get all required lists
		d_range = range(len(w_data_tuple))
		for elem in d_range:
			if isinstance(w_data_tuple[elem],str):
				w_data_tuple[elem] = self.get_list (w_data_tuple[elem], w_lib)
			elif isinstance(w_data_tuple[elem],list): pass	# allow user-def list
			else: error(" unknown types for w_data_tuple")
		
		# write out
		for i in mIndex:
			f_str = ""
			for elem in d_range:
				f_str += w_str_tuple[elem]
				f_str += str(w_data_tuple[elem][i])
			f_str += w_str_tuple[-1]
			self.write(w_file, f_str, end_wrr_, beg_wrr_)		

	# change list according to a pattern
	# could also use RHS, to define (ODE) * dt etc
	# non-return is not defined yet - safe not to,to avoid overwritting self.data
	#def pattern_change (self, oList, mIndex, specifier = ["",""], return_it = True):
	def pattern_change (self, gIndex , orig,  *argv):
		if gIndex == None:
			mIndex = range (len(orig))
		else:
			mIndex = gIndex
		a = []
		if isinstance(argv[0],list): m_size = len(argv[0])
		if isinstance(argv[1],list): m_size = len(argv[1])
		for it in range(m_size):
			if it in mIndex:
				f_str = ""
				for elem in argv:
					if isinstance(elem,str):
						f_str += elem
					else: f_str += elem[it]
				a.append(f_str)
			else:
				a.append(orig[it])
		return a
			
	def new_index (self, *argv):
		return  list(set(argv[0]).intersection(*argv))
	
	def get_index (self, index_to_get):
		return [self.names.index(elem) for elem in index_to_get]
		
	
	# %[unviversal] index_order, according to another ordered index range
	def order_index (self, index_to_order, order_by = None):
		if order_by is None:
			order_by = self.index['order']['init-value']
		return list(filter(lambda elem: elem in index_to_order, order_by))
		

	def new_dependancy_index(self, index_to_order, key_to_order, is_return = True):
		ko = re.split
		
		# restructure index
		#if is_return is True:
		#	r_index = copy.deepcopy(index_to_order)
		#else:
		#	r_index = index_to_order
		#print index_to_order
		if is_return is False:
			r_index = index_to_order
		
		r_index = list(filter(None, [iter if self.data[iter][key_to_order] != self.null else None for iter in index_to_order]))

		d_name = self.get_list ("name", self.data,"","", r_index)

		dependencies = 	self.get_list (key_to_order, self.data,"","", r_index)

		##### dependencies = [[x for x in y if not x[0].isdigit()] for y in \
		#		[filter(None, elem) for elem in \
		#		[ko(r'(?:\.|\^| |\+|\*|-|/|\)|\(|\]\]|\[\[|\|\||,|log10|exp|log|| )',elem) for elem in dependencies]]]
		dependencies = [[x for x in y if not x[0].isdigit()] for y in \
				[[i for i in elem if i] for elem in \
				[ko(r'(?:\.|\^| |\+|\*|-|/|\)|\(|\]\]|\[\[|\|\||,|log10|exp|log|| )',elem) for elem in dependencies]]]


		counter_d = 0
		for i in range(len(r_index)-1,-1,-1):
			if dependencies[i]:
				r_index.append (r_index.pop(i))
				d_name.append (d_name.pop(i))
				dependencies.append (dependencies.pop(i))
				counter_d += 1

		counter_d = len(dependencies) - counter_d	
		while counter_d < len(r_index):
			if not set(dependencies[counter_d]).intersection (d_name[counter_d:]):
				counter_d += 1
			else:
				r_index.append (r_index.pop(counter_d))
				d_name.append (d_name.pop(counter_d))
				dependencies.append (dependencies.pop(counter_d))

		return r_index

		
	# @[universal] define a new preserved order	
	def new_order (self, *argv):
		check = []
		for vName in argv:
			if vName not in self.names:
				check.append(vName)
		if check:
			error("in function new_order, these names are not defined in list " + str(check))	
		return [self.names.index(vName) for vName in argv]

	# only works for names in data set???
	def new_names (self, oNames, mIndex, specifier = ["",""], return_it = True):
		names_list = []
		if len(specifier) == 1: specifier[1] = ""
		for i in range(len(oNames)):
			if i in mIndex:	
				names_list.append(specifier[0] + oNames[i] + specifier[1])
			else:			
				names_list.append(oNames[i])
		return a

		
	def new_names_restricted (self, oNames, mIndex, specifier = ["",""], return_it = True):
		a = []
		if len(specifier) == 1: specifier[1] = ""
		for i in range(len(oNames)):
			if i in mIndex:	
				a.append(specifier[0] + oNames[i] + specifier[1])
		return a
				
	def get_only_names (self, oNames, mIndex):
		a = []
		for elem in mIndex:
			a.append(oNames[elem])
		return a
		
	def vectorise_name (self, mIndex, specifier, names_list = None, return_it = True):
		
		counter = range(len(mIndex),self.vec_counter_start-1,-1)
		
		# change name
		if names_list is None:	# if there is no names list

			vnames = [specifier + self.l_enclose + str((list(counter)).pop()) + self.r_enclose if x in mIndex else self.names[x] for x in self.info['range']]

		else:

			vnames = [specifier + self.l_enclose + str((list(counter)).pop()) + self.r_enclose if x in mIndex else names_list[x] for x in self.info['range']]			
			
		if return_it is True:
			return vnames;
	
	def write (self, fileAttr, strAttr, end_wrr_ = "\n", beg_wrr_ = ""):
		#if isinstance (fileAttr, file): f = fileAttr
		f = fileAttr
		if isinstance (fileAttr, str): f = self.open_files [fileAttr]
		else: error("fillAttr")
		if isinstance (strAttr, str):
			comb_str = beg_wrr_ + strAttr + end_wrr_
			f.write(comb_str) # could use f.write(bytes(comb_str,'UTF-8'))
		elif isinstance (strAttr, list) and isinstance (strAttr[0], str):
			f.write( beg_wrr_)
			for strElem in strAttr:
				if strElem is not "": 	f.write( strElem + "\n")
				else: 					f.write( "\n")
			f.write( end_wrr_ )
								
	def write_comment (self, fileAttr, strAttr, end_wrr_ = "\n", beg_wrr_ = ""):
		if isinstance (strAttr, str): 	strAttr = self.comment_prefix + " " + strAttr
		elif isinstance (strAttr, str):	strAttr = [self.comment_prefix + " " + elem for elem in strAttr]
		self.write (fileAttr, strAttr, end_wrr_, beg_wrr_)
		
	
	def get_names (self):
		return self.names

	def rasterise_dict_modifiers (self, refi):
		return None
		
	#@classmethod
	#def info(cls):
	#""" Print description of parser to screen"""
	#	print  
	
	#@clasmethod
	#def make_parser(cls):