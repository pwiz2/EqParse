#! /usr/bin/env python

"""
	create Library of variables : Read saved model file (current formats: csv, )
"""

import csv
import re
import copy

from smc_helper_functions import *
from operator import itemgetter

from collections import defaultdict
from collections import Counter
from collections import OrderedDict

import os.path
import time

genabspath = os.path.abspath
getmtime = os.path.getmtime
ctime = time.ctime


# change a key (vKey) of an element name (vName) from variable library (Lib) to something else (vToo)
#def change_a_type (Lib, vName, vKey, vToo):
#	index = Lib.names.index(vName)
#	Lib.data[index][vKey] = vToo;
#	Lib.data_vec[index][vKey] = vToo;


from __init__ import __eqp_mac_address as macaddress
from __init__ import __eqp_memory_file as memfilename
from __init__ import __eqp_memory_file_comb as memfilename_comb

import pickle

 
class CreateLibrary (object):
	"""
	This is the reader for the mathematical system parser. 
	Before calling a parser module, one must call this class with the csv
	files to be parsed, organised and managed into a dictionary for easy look-up
	and readability for the different parsers. Currently this module only accepts 
	CSV files.
	The future long-term goal is the creation of other possible formats such as XML and
	possibly a GUI interface.
	
	usage: 
	"""
	#__metaclass__ = ReadOnlyClass
	
	
	table = []
	def __init__(self, model_name=None, model_file_list=None):
		
		if model_file_list is None:
			error("Must include name of file list")

		if model_file_list is None:
			error("Must include file(s) of master model file")
		
		self.full_path = genabspath(model_name)
		self.model_title = model_name.rsplit('/')[-1]
		self.directory = genabspath(model_name.rsplit('/',1)[0])
		self.nofl = model_file_list	
		self.null = 'NULL'
		self.table = []		

		#self.ischanged = range(len(model_file_list)) # change to single line code - all 1 value
		
		#for i in list(range(len(self.ischanged))):
		#	self.ischanged[i] = 0
		self.ischanged = [0] * len(model_file_list)
		
		# determines if we regenerate the parsed files if they have been edited
		# total value must equal len of input files if we choose NOT to update
		REGENERATE = 0 # default: set to not generate files

		fprop = pickle.load(open(memfilename, 'rb')) #-- get dictionary of same name values
		fcomb = pickle.load(open(memfilename_comb, 'rb')) #-- get dictionary of combinations done
		if fprop is None:
			fprop = dict()
			pdump(fprop,memfilename) 

		if fcomb is None:
			fcomb = dict()
			pdump(fcomb,memfilename_comb) 
			
		for i in range(len(model_file_list)):
			current_fprop_path = genabspath(genabspath(model_file_list[i])).replace('/','_').replace('/','_')[1:]
			#-- get time stamp values
			current_fprop_tstamp =  getmtime(genabspath(model_file_list[i]))
			if current_fprop_path in fprop: # if this item exist in dict, compare to current value
				if fprop[current_fprop_path] == current_fprop_tstamp:
					REGENERATE = REGENERATE + 1
					self.ischanged[i] = 1
					#print "same time stamp..."
				else: # !vital! remove all combinations of model_file_list[i] as this file has changed
					for k in fcomb.keys():
						if model_file_list[i] in list(fcomb[k]):
							del fcomb[k]
			else:
				fprop[current_fprop_path] = current_fprop_tstamp;
			
		#-- gather variables from the files combined - remove duplicates (how does it do this?)
		try:				
			for i in range(len(model_file_list)):
				fn_extension = model_file_list[i].rsplit('.',1)[1]
				if fn_extension == "csv":
					self.table = self.table + [row for row in csv.reader(open(model_file_list[i],'rU'), delimiter=',') if row != []]
				elif fn_extension == "another_extension":
					pass
				else:
					error ("incorrect extension")
		except:
			error ("Error: incorrect filename and/or extension given.\n:: " + str(model_file_list) + " does not exist")
		
		# 0. 
		# SOF: How do you remove duplicates from a list in python whilst preserving order?
		#self.table
		allthenames=[row[1] for row in self.table if row[0][0] != '%']		
		
		# 1. store title, 
		self.info = {}; self.names = []
		if REGENERATE == len(model_file_list) and (self.full_path in fcomb):
			self.info['doWeParse'] = 0;
			print(":: NOT Regenerated")
		else:
			self.info['doWeParse'] = 1;
			print(":: Regenerated")
			
		pickle.dump(fprop,open(memfilename, 'wb'))
		pickle.dump(fcomb,open(memfilename_comb, 'wb'))

		self.info['title'] = self.model_title #self.table[0][1]
		self.info['keys'] = ['type', 'name', 'comment', 'dimension', 'init-value', 'print', 'userdef', 'rhs']
		self.info['const_keys'] = ["type", "name"]
		self.info['defualt'] = {'comment': "{unknown}",  # default value for each key
							'dimension': "n/a", 
							'init-value': '0', 
							'print': 'NO', 
							'userdef': 'NO', 
							'rhs': 'NULL'} 
		self.info['function'] = []	# function order
		
		# 2. generate data as a list of dictionaries
		self.data_store = [self.create_data()]
		
		self.data = self.data_store[0]
		
		self.d_counter = False
		self.size_vars = 0		#2aa size_name is size of total variables
		for elem in self.data:
			self.names.append(elem['name']) # 2ab
			self.size_vars += 1	# 2aa
		self.info['range'] = range(len(self.names))
		
		self.vec = [] # vectorized data
		self.info['total'] = len(self.data)
		self.info['index'] = {'order': {}}
		self.add_index (['ode','yes'],{'type': 'ODE'})
		self.add_index (['ode','no'],{'type': 'VAR'})
		self.add_index (['var','yes'],{'type': 'VAR'})
		self.add_index (['var','no'],{'type': 'ODE'})
		self.add_index (['userdef','yes'],{'userdef': 'YES'})
		self.add_index (['userdef','no'],{'userdef': 'NO'})
		self.add_index (['print','yes'],{'print': 'YES'})
		self.add_index (['print','no'],{'print': 'NO'})
		self.add_index (['rhs','yes'],{},{'rhs': 'NULL'})
		self.add_index (['rhs','no'],{'rhs': 'NULL'})		
		self.add_index (['init-value','yes'],{'init-value': 'NULL'})
		self.add_index (['init-value','no'],{},{'init-value': 'NULL'})
		
		self.info['index']['order']['init-value'] = self.get_index_dependency ('init-value', 'name')

		pass
		
	def change_variable (self, mname, mval):
		m = self.names.index(mname)
		self.data[m]['init-value'] = str(mval)
		for elem in range(len(self.data_store)):
			self.data_store[elem][m]['init-value'] = str(mval)	
		
	# add current time, one can add only before defining parsers	
	def add_variable (self, mDict):
		orig_keys = self.info['keys']
		
		# check if there are key properties which dont exit in our library
		if set(mDict.keys()).issubset(orig_keys) is False: 
			error("Keys defined in var dict are not part of the original values\nMust be part of the following attributes: " + str(orig_keys))
			
		# check if values are correct,
		for key in mDict.keys():
			if key == 'type':
				if key not in ['ODE','VAR'] : 
					error('var creation: \'type\' (is \'ODE\' or \'VAR\')')
			elif key == 'name':
				if mDict['name'][0].isdigit(): 
					error('var creation: \'name\'')

	
	def function_order (self, *argv):
		if not argv:
			None
		else: 
			names_rhs = []
			for it in self.info['index']['rhs']['yes']:
				names_rhs.append (self.names[it])
			arg_names = [item for sublist in argv for item in sublist['rhs']]
			new_occ =  list(set(arg_names) - set(names_rhs))	# these items are not declared in RHS
			old_occ = list(set(names_rhs) - set(arg_names))		# these items are apparent new items
			if new_occ:
				error('(function_order) - The following items are not declared as a rhs function in library, please check the spelling, or if this variable has a rhs defined: ' + str(new_occ))
			elif old_occ:
				error('(function_order) - The following items you defined does not exist: ' + str(old_occ))
			
			for i in range(len(argv)):
				argv[i]['index'] = [self.names.index(elem) for elem in argv[i]['rhs']]
			self.info['function'] = argv

	# there may be faster methods to find duplicates, as in this SOF page
	# http://stackoverflow.com/questions/5419204/index-of-duplicates-items-in-a-python-list
	# http://stackoverflow.com/questions/716477/join-list-of-lists-in-python
	# this one is modified
	def dup2 (self, n): #n="123123123"
		counter = Counter(n) #{'1': 3, '3': 3, '2': 3}
		dups = [i for i in counter if counter[i]!=1] #['1','3','2']
		return sorted([inner for outer in [i[1:] for i in 
				[[i for i,j in enumerate(n) if j==item] for item in dups]] for inner in outer],
				key=int,reverse=True)
                        
	def create_data (self):
		a = [dict(zip(self.info['keys'],values)) for values in self.table if not values[0][0] == "%"]
	
		t_names = [x['name'] for x in a]
		b2 = self.dup2(t_names)
		for elem in b2:
			a.pop(elem)
		return a
		
	def copy_lib (self):		
		self.data_store.append(self.create_data())
		return self.data_store[-1]
		
	def copy_lib_vec (self):		
		mkeys = self.info['keys']
		self.vec.append (self.create_data())
		for elem in self.info['index']['userdef']['yes']:
			self.vec[-1][elem]['init-value'] = self.data[elem]['init-value']
		return self.vec[-1]
		
	def get_index_dependency (self, 
			m_key_to_order = 'init-value', 
			resp_order_val = 'name', 
			cond_remove = {}, 
			cond_include = {}):
		# SOF:: How to find all occurenceies of an element in list?   <-- figure how to do this in list coprehension (?)
		
		if isinstance(m_key_to_order, basestring):
			key_to_order = m_key_to_order
		ko = re.split
		dependencies = filter(None, [entity[key_to_order] if entity[key_to_order] != self.null else '' for entity in self.data])

		dependencies = [[x for x in y if not x[0].isdigit()] for y in \
				[filter(None, elem) for elem in \
				[ko(r'(?:\.|\^| |\+|\*|-|/|\)|\(|\]\]|\[\[|\|\||,|log10|exp|log|| )',elem) for elem in dependencies]]]
		d_name = [entity[resp_order_val] if entity[key_to_order] != self.null else '' for entity in self.data]
		r_index = []
		for elem in d_name:
			if elem != '':
				r_index.append(d_name.index(elem))
		if cond_remove:
			counter_d = 0
			r_index_len = len(r_index);
			while counter_d < r_index_len:
				for mkey in cond_remove.keys():
					if str(self.data[r_index[counter_d]][mkey]) == str(cond_remove[mkey]):
						r_index.pop (counter_d)
						dependencies.pop (counter_d)
						d_name.pop (counter_d)
						r_index_len -= 1
					else:
						counter_d += 1

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

		
	def add_index (self, define_dict_keys , cond_on = {}, cond_not_on = {}):
		
		indexed_values = self.info['index'];
		for elem in define_dict_keys[:-1]:
			if elem not in indexed_values.keys():
				indexed_values[elem] = {}
			indexed_values = indexed_values[elem]			
		indexed_values[define_dict_keys[-1]] = []
		
		indexes_on = []; indexes_not_on = []
		
		seen = set(); seen_add = seen.add		
		for key in cond_on.keys():
			if key in self.data[0].keys():
				indexes_on += [i for i, x in enumerate(self.data) if x[key] == cond_on[key]]
			else: error("stop")
		indexes_on = [x for x in indexes_on if x not in seen and not seen_add(x)]
		
		seen = set(); seen_add = seen.add
		for key in cond_not_on.keys():
			if key in self.data[0].keys():
				indexes_not_on += [i for i, x in enumerate(self.data) if x[key] == cond_not_on[key]]
			else: error("stop")
		indexes_not_on = [x for x in indexes_not_on if x not in seen and not seen_add(x)]			
		
		if not cond_on:
			indexes_on = range(0,self.info['total'])
		
		#for elem in indexes_not_on:
		#	if elem in indexes_on:
		#		indexes_on.remove(elem)
		indexes_on = [elem for elem in indexes_on if elem not in indexes_not_on]
		indexed_values[define_dict_keys[-1]] = indexes_on


	# note: need to run this after ever library is finished
	"""
	overide set_directory : 
	"""	
	def set_directory(self):
		"""Set directory where file specific parser files are saved to, 
		the directory can be a relative path and not necessarily an absolute path
		
		:param dir_str: string of path to directory
		"""
		if isinstance(self.directory, str) is False: error("set_directory argument not a string")
		elif not self.directory: pass
		else: 
			if not os.path.exists(self.directory + "/" + self.model_title):
				os.makedirs(self.directory + "/" + self.model_title)	

				
	def complete(self):
		fcomb = pickle.load(open(memfilename_comb, 'rb')) #-- get dictionary of combinations done
		fcomb[self.full_path] = self.nofl
		pickle.dump(fcomb,open(memfilename_comb, 'wb'))
		if self.info['doWeParse'] == 1:
			fprop = pickle.load(open(memfilename, 'rb'))
			for i in range(len(self.nofl)):
				#if self.ischanged[i]==1
				current_fprop_path = genabspath(genabspath(self.nofl[i])).replace('/','_').replace('/','_')[1:]
				current_fprop_tstamp =  getmtime(genabspath(self.nofl[i]))
				fprop[current_fprop_path] = current_fprop_tstamp;
			pickle.dump(fprop, open(memfilename, 'wb'))
