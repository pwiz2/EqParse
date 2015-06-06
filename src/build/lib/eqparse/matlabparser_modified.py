#! /usr/bin/env python

from baseparse import *

"""
MATLAB PARSER
"""

class MatlabParserModified (BaseParse):
	"""
	To create matlab parsed file
	"""
	
	#__metaclass__ = ReadOnlyClass


	def __init__ (self, Lib):
		super(MatlabParserModified, self).__init__()
		self.l_enclose = "("
		self.r_enclose = ")"
		self.operations = ['\nif(', ")\n", "else\n\t", "elif\n", "end\n", "~=", "(", ").^(", ")","sqrt("]
		self.vec_counter_start = 1
		self.comment_prefix = "%"		
		self.fileExt = "m"		
		
		self.initialise_library (Lib)
		# set title as folder name in matlab directory,...
		# ... and create if not exist
		#s_direct = self.directory + "/" + self.title;
		#if not os.path.exists (s_direct):
		#	print "::::::::::::wtf" + s_direct
		#	os.makedirs (s_direct)
		#else:
		#	print "WARNING: over writing folder " + self.title
		#pass
	
	
		
	"""
	CREATE SINGLE-FILE MATLAB FILE FUNCTION
	"""
	def solve_single(self):
		"""generate single MATLAB (single_generic.m) file to solve equations. Runs the Euler finite differece scheme to model with  dt=5 (ms)	
		"""	
		
		# simplified variables
		mfile = "single_generic"
		mfio = self.title + "/" + mfile
		mcomment = self.write_comment
		mwrite = self.write
		set_p = "set_param"
		total_time = "tot_secs"
		s_mat = "e_a"					# storage matrix
		mget_result = "S"
		
		# open single file
		self.open_file (mfio)
		
		# write header
		mwrite (mfio, ["function " + mget_result + " = " + mfile + " (" + set_p + "," + total_time + ")", \
							"close all"], "\n\n")
		
		# get_list only used for 
		#self.get_list (self, mKey, mLib = None, mPrefix ="", mPostfix ="", mIndex = None):
		mcomment (mfio, "gather arguments from list and define settings")
		mwrite (mfio, ["dt = 5;", "t =0:dt:" + total_time + "*1000;" ])
		
		# X = set_param.X and init-values direct
		mcomment (mfio, "gather user defined variables")
		self.pattern_write (mfio, ["name", "name", "comment"], ["", " = " + set_p +".", "; % ", ""],
									self.data, self.index['userdef']['yes']) # ode can't be user de	
		mIndex = self.order_index (self.index['userdef']['no'])	
		self.pattern_write (mfio, ["name", "init-value", "comment"], ["", " = ", "; % ", ""],
									self.data, mIndex)	
		#rhs_iter = self.get_list ("rhs", self.data, "", "", None)
		#rhs_init = self.get_list ("init-value", self.data, "", "", None)
		#init_list = [rhs_iter[it] if it in self.index['ode']['no'] else rhs_init[it] for it in self.info['range']]

		# allocate storage							
		mcomment (mfio, "allocate space for printed variables (including ODEs)", " \n\n")
		mwrite (mfio, s_mat + " = zeros (1, size(t,2)-1);")
		
		mcomment (mfio, "First all the functions...", " \n", "\n")
		mIndex = self.new_index (self.index['ode']['no'], self.index['print']['yes'])
		self.pattern_write (mfio, ["name","name"], ["", " = [ ", ", " + s_mat + "];"], 
									self.data, mIndex)
		
		mcomment (mfio, "... then ODE storage", " \n", "\n")
		mIndex = self.new_index (self.index['ode']['yes'], self.index['print']['yes'])
		self.pattern_write (mfio, ["name","name"], ["", " = [ ", ", " + s_mat + ", 0];"], 
									self.data, mIndex)
		
		# perform iterated simple step
		mwrite (mfio, ["max_iter = size(t,2);",
						"for iter = 1:max_iter"], "\n", "\n\n")
							
		index_print_nonode = self.new_index (self.index['ode']['no'], self.index['print']['yes'])
		lhs_names_list = self.new_names(self.new_names(self.names, index_print_nonode, ["","(iter)"]), self.index['ode']['yes'], ["","(iter+1)"])		
		for i in range(len(self.info['function'])):
			f_index = self.info['function'][i]['index']
			f_var_name = self.info['function'][i]['rhs']
			lhs_name = self.get_only_names (lhs_names_list, f_index)
			rhs_names = self.new_names(self.names, self.index['print']['yes'], ["","(iter)"])
			rhs_list = self.get_list ("rhs", self.data, "", "", f_index)
			rhs_list = self.search_and_replace (rhs_list, rhs_names + ["t(iter)"], self.names + ["t"], False)	
			for e in range(len(f_var_name)):
				if f_index[e] in self.index['ode']['yes']:
					rhs_list[e] = f_var_name[e] + "(iter) + ( " + rhs_list[e] + ")*dt"
			#print rhs_list		
			mcomment (mfio, self.info['function'][i]['comment'],"\n","\n")
			self.pattern_write (mfio, [lhs_name,rhs_list], ["", " = ", ";"])
			#def pattern_change (self, mIndex , orig,  *argv):
			
		## print rhs_list
		#index_print_nonode = self.new_index (self.index['ode']['no'], self.index['print']['yes'])
		#lhs_names = self.new_names(self.new_names(self.names, index_print_nonode, ["","(iter)"]), self.index['ode']['yes'], ["","(iter+1)"])
		## print rhs_list
		#rhs_names = self.new_names(self.names, self.index['print']['yes'], ["","(iter)"])
		#rhs_list = self.get_list ("rhs", self.data, "", "", None)

		#rhs_list = self.search_and_replace (rhs_list, rhs_names + ["t(iter)"], self.names + ["t"], False)
		#rhs_list = self.pattern_change (self.index['ode']['yes'], rhs_list, self.names , "(iter) + ( ", rhs_list ," )*dt") ### this only works for all elements

		#self.pattern_write (mfio, [lhs_names,rhs_list], ["", " = ", ";"],	self.data, self.index['rhs']['yes'])
									
		# reshape ode and time storage - well defined		
		mwrite (mfio, "end", "\n\n", "\n")
		self.pattern_write (mfio, ["name","name"], ["", " = ", "(1:end-1);"], 
									self.data, mIndex)
										
		mwrite (mfio, mget_result + ".('t') = t./1000; " + self.comment_prefix + " change to seconds")
		mIndex = self.order_index (self.index['print']['yes'])
		self.pattern_write (mfio, ["name", "name"], [ mget_result + ".('", "') = ", ";"], self.data, mIndex)
									
		# close single file
		self.close_file (mfio)



	"""
	CREATE RUN-FILE MATLAB FILE FUNCTION
	contains the actual set-up and organisation of cell system
	in a structured order, is also main file
	"""
	def multi_runfile (self, newSubDirectory = "", 
						file_ext = "", 
						no_guess = None, 
						yes_guess = None, 
						wr_modifiers = ""):
		"""Create run file (run.m) that uses the ode file created from running alongside with :py:mod:`multi_odefile` module. This format is currently set for MATLAB defined implicit ode function - ode15s. 
		
		:param newSubDirectory:
		:param file_ext: 
		:param no_guess: 
		:param yes_guess: 
		:param wr_modifiers:	
		"""	
								
		#mfile = "run" + file_ext
		mfile =  "run" + file_ext
		mfio = self.title + "/" + newSubDirectory + "run" + file_ext
		mcomment = self.write_comment
		mwrite = self.write
		m_print = "p_store"	# VC_NAME_PRINT
		m_noprint = "p"	# VEC_NAME
		m_ode = "dydt"	# ode var ODE_NAME
		m_out = "out"
		m_otime= 'out_time'
		t_end = "t_end"
		set_p = "set_params"
		ode_method = "ode15s"
		m_ode_fname = "ode"
					
		# open file
		self.open_file (mfio)
		
		# define header
		if yes_guess is None: mp = m_out
		else: mp = "b_out"
		mwrite (mfio, "function " + mp + " = " + mfile + " (" \
							+ set_p + ", " + t_end + ", b_in)", "\n\n")	
		
		
		# preallocate storage
		mI_user_no_ode = self.index['print']['no']
		mcomment (mfio, "preallocate storage for data")
		mwrite (mfio, [m_noprint + " = zeros(1," + str(len(mI_user_no_ode)) + ");",
						m_print + " = zeros(1," + str(len(self.index['print']['yes'])) + ");",
						m_ode + " = zeros(1," + str(len(self.index['ode']['yes'])) +");"])

		# if parameterizing it, create initial data from external function
		# and define parameter initial values also
		if no_guess is not None:
			mcomment (mfio, "Define sample data to fit model here ()","\n","\n")
			mwrite (mfio, "[tdata,ydata] = create_ic (t_end);", "\n\n", "")
			mwrite (mfio, "\nif nargin<3")
			mName = [elem for elem in self.names]; cntr = 1;
			for elem in yes_guess:
				mName[elem] = "b(" + str(cntr) + ")"
				cntr = cntr + 1
			self.pattern_write (mfio, [mName, "init-value", "comment"], [ " ", " = ", "; \t% ", ""], self.data_vec, yes_guess)
			mwrite (mfio, "else\n\tb=b_in\nend")
			
		# ..user defined args
		
		##struct_vec_names = self.get_list ("name", self.data_vec)
		##t_end_vec = struct_vec_names[self.names.index(t_end)]
		#struct_vec_names[index_tend] = t_end
		struct_vec_names = "name"
		
		struct_names = self.new_names(self.names, self.index['userdef']['yes'], [set_p + ".",""])
		self.pattern_write (mfio, [struct_vec_names, struct_names, "comment"], [ " ", " = ", ";% ", ""], self.data_vec, self.index['userdef']['yes'])		
		# .. function dependent
		mIndex = self.order_index(self.index['userdef']['no'])
		self.pattern_write (mfio, ["name", "init-value", "comment"], [ " ", " = ", ";% ", ""], self.data_vec, mIndex)

		# predefine time, if no input given
		mwrite (mfio, ["if nargin<2", "\t" + t_end + " = " + set_p + "." + t_end + ";", "end\n"])

		# start ode
		if no_guess is None:
			mwrite (mfio, ["options = odeset('InitialStep',0.02,'MaxStep',25, 'RelTol', 1e-9, 'AbsTol', 1e-9);",
							"[t," + m_ode + "] = " + ode_method + "(@(tt,y)" + m_ode_fname + "(tt,y," 
									+ m_noprint + "," + m_print + ")," +
								"[0," + t_end + "*1000]," + m_ode + ",options);"],"\n","\n\n")
			# output resultsinto a struct
			mwrite (mfio, ["[~," + m_print + " ] = " + m_ode_fname + "([],[],[],[]);",
					"[~, index_max_t_] = max(" + m_print + "(:,end));",
					#m_print + " = " + m_print + "(1:index_max_t_+1,:);",
					m_print + " = " + m_print + "(1:index_max_t_,:);",
					m_out + "= struct;",
					m_out + ".ode = t./1000;",	#struct;",
					m_out + ".var = " + m_print + "(:,end)./1000;"])
		else:	# is parameterised
			mwrite (mfio, "[bmin, Smin] = fminsearch(@(b)" + self.title + "_guess" + "(b, tdata, ydata, "+ m_ode + "," + \
			 						m_noprint+","+m_print + "),b);\ndisp(Smin)","\n","\n\n")										
			mwrite (mfio, mp + " = struct;")
			for elem in yes_guess:
				mName[elem] = "bmin" + mName[elem][mName[elem].find('('):]
			self.pattern_write (mfio, ["name",mName],[mp + ".('","') = ","; "], self.data, yes_guess)
			
		
		lhs = self.get_list ("name", self.data_vec, "", "", None)	# get a column form dict list
		num_iter = [filter(str.isdigit, elem) for elem in lhs]		# get number from derived list
		sto_iter = [elem[:elem.index('(')] for elem in lhs]		# get the name of storage for each
		if no_guess is None:
			# output resultsinto a struct
			#mwrite (mfio, ["[~," + m_ode + "] = " + self.title + "_ode" + "([],[],[],[],[],[]);",
			#		m_out + "= struct;",
			#		m_out + ".ode = t./1000;"])		
			index_ode = self.new_index(self.index['ode']['yes'], self.index['print']['yes'])		
			self.pattern_write (mfio, ["name", sto_iter, num_iter], [ m_out + ".('", "') = ", "(:,",");"], self.data, index_ode)
		#mwrite (mfio, m_out + ".ode.('t') = t./1000;", "\n\n")		
			index_var = self.new_index(self.index['ode']['no'], self.index['print']['yes'])
			self.pattern_write (mfio, ["name", sto_iter, num_iter], [ m_out + ".('", "') = ", "(:,",");"], self.data, index_var)		
		#mwrite (mfio, m_out + ".var.('t') = " + m_print + "(:,end)./1000;", "\n\n")
		
		self.close_file (mfio)





	"""
	CREATE ODE-FILE MATLAB FILE FUNCTION
	the actual ODE system set up as a separate m file
	for passing into an ODE system
	"""	
	def multi_odefile (self, newSubDirectory = "",
						file_ext = "", 
						no_guess = None, 
						yes_guess = None, 
						wr_modifiers = ""):
		"""Create matlab ODE file (ode.m) to be used with :py:mod:`multi_runfile` and :py:mod:`multi_runfile_slider` modules"""

						
		mfile = "ode" + file_ext
		mfio = self.title + "/" + newSubDirectory +  "ode" + file_ext
		mcomment = self.write_comment
		mwrite = self.write
		m_print = "p_store"	# VC_NAME_PRINT
		m_noprint = "p"	# VEC_NAME
		m_ode = "dydt"	# ode var ODE_NAME
		m_ode_rhs = "y" # rhs dydt as y
		mpersistent_out = "out_" + m_print
		mpersistent_temp = "pers_" + m_print
		len_of_odes = len(self.new_index(self.index['ode']['yes']))
		
		# open file
		self.open_file (mfio)

		# define header (function name etc)
		if yes_guess is None: new_arg = ""; temp = ""; jo = "," + mpersistent_out
		else: new_arg = "b"; temp = "," + new_arg; jo = ""
		mwrite (mfio, "function [" + m_ode +  jo + "] = " + mfile + 
										"(t," + m_ode_rhs + "," + m_noprint + "," + m_print + temp + ")", "\n\n")
										

		
		
		if yes_guess is None:								
			mwrite (mfio, ["persistent " + mpersistent_temp + " counter",
							"if isempty(" + mpersistent_temp + ")", 
							"\t" + mpersistent_temp + " = zeros(200,size(" + m_print + ",2));\n\tcounter = 1;",
							"end"], "\n" )
			mwrite (mfio, ["if mod(counter-1,200) == 0",
							"\t" + mpersistent_temp + " = [" + mpersistent_temp + "; zeros(200,size(" + m_print + ",2))];",
							"end"], "\n")		
		
			mcomment (mfio, "solve the ode step...")
			
			mwrite (mfio, [m_ode + " = zeros(" + str(len_of_odes) + ",1);"])
			mwrite (mfio, "if (nargout == 1)")

		if yes_guess is None:
		# if guess
			oderhs_names = self.vectorise_name (self.new_index(self.index['var']['yes'],self.index['print']['yes']), "p_store")
			oderhs_names = self.vectorise_name (self.new_index(self.index['var']['yes'],self.index['print']['no']), "p", oderhs_names)
			oderhs_names = self.vectorise_name (self.index['ode']['yes'], "y", oderhs_names)
			m_rhs = self.search_and_replace (self.get_list ("rhs", self.data, "", "", None), oderhs_names, self.names)
		else:
			yes_guess_names = self.get_only_names(self.names, yes_guess)
			no_guess_names = self.get_only_names(self.names, yes_guess)
			# vectorize data
			vnames = self.vectorise_name (self.new_index(self.index['var']['yes'],self.index['print']['yes']), "p_store")
			vnames = self.vectorise_name (self.new_index(self.index['var']['yes'],self.index['print']['no']), "p", vnames)
			vnames = self.vectorise_name (self.index['ode']['yes'], "dydt", vnames)
			cntr = self.vec_counter_start
			for elem in yes_guess:
				vnames[elem] = "b(" + str(cntr) + ")"
				cntr = cntr + 1			
			m_rhs = self.search_and_replace (self.get_list ("rhs", self.data, "", "", None), vnames, self.names)
			
			
		mwrite (mfio, m_print + "(end) = t;")
		setODEFormatAs = 2
		# 1: old format, based on function
		# 2: new format, ordering dependent
		if setODEFormatAs == 1:
		# ORIGINAL ODE/PARAM UPDATE
			for i in range(len(self.info['function'])):
				mcomment (mfio, self.info['function'][i]['name'],"\n","\n")
				self.pattern_write (mfio, ["name", m_rhs, "comment"], [ " ", " = ", "; %", ""], self.data_vec, self.info['function'][i]['index'])
		elif setODEFormatAs == 2:
			index_nonode = self.new_dependancy_index(self.new_index(self.index['ode']['no'],self.index['rhs']['yes']),'rhs')
			#for indexer in index_nonode:
			#	print self.data[indexer]['rhs']
			# print all function related dependencies in front of order
			self.pattern_write (mfio, ["name", m_rhs, "comment"], [ " ", " = ", "; %", ""], self.data_vec, index_nonode)
			# print all ODE related dependencies out to file in last order
			self.pattern_write (mfio, ["name", m_rhs, "comment"], [ " ", " = ", "; %", ""], self.data_vec, self.index['ode']['yes'])			
			#new_dependancy_index	
			#mIndex_nonode = self.order_index (self.new_index(self.index['ode']['no'],self.index['rhs']['yes']))
			#mIndex_ode = self.order_index (self.index['ode']['yes'])
			#self.pattern_write (mfio, ["name", m_rhs, self.names, "comment"], [ " ", " = ", "; % {", "} ", " "], self.data_vec, mIndex_nonode)
			#self.pattern_write (mfio, ["name", m_rhs, self.names, "comment"], [ " ", " = ", "; % {", "} ", " "], self.data_vec, mIndex_ode)
		else:
			None		
		
		
		if yes_guess is None:
			mcomment (mfio, "...otherwise sort the persistent storage out", "\n", "\n\n")
			mwrite (mfio, [mpersistent_temp + "(counter,:) = " + m_print + ";",
							"counter = counter + 1;"])
			mwrite (mfio, ["else",
						"\t " + mpersistent_out + " = " + mpersistent_temp + ";",
                                                "\t" + mpersistent_temp + " = [];",
						"end"])
						
			mwrite (mfio, "end % end ode function")		
				
		self.close_file (mfio)




		
		






	"""
	CREATE PARAM-FILE MATLAB FILE FUNCTION
	is a function which returns all the parameters defined to pass into
	other functions where parameters need defining
	Currently used in single-file and multi-run-file m files
	"""		
	def multi_paramfile (self, newSubDirectory = "", file_ext = "", no_guess = None, yes_guess = None):
		
		#if no_guess is None: no_guess = self.['range']
		#if no_guess is not None:
		#	print len(no_guess)
		# helper variables
		#mfile = "param" + file_ext
		mfile = "param" + file_ext
		mfio = self.title + "/" + newSubDirectory + "param" + file_ext
		mcomment = self.write_comment
		mwrite = self.write
		scalor_p = "set_param"
		info_p = "print_info"
		
		# open file
		self.open_file (mfio)
		
		mcomment (mfio, "initialize the original parameter output", "\n",)
		mwrite (mfio, ["function [" + scalor_p + "," + info_p +  "] = " + mfile, 
						scalor_p + " = struct;",
						info_p + " = struct;"])
		# create initialized parameter struct
		#mIndex = self.index['userdef']['yes']
		#self.pattern_write (mfio, ["name","init-value"], [scalor_p + ".('", "') = ", ";"], 
		#							self.data, mIndex)
		
		if no_guess	is None: 	mIndex = self.index['order']['init-value'];
		else:					mIndex = self.order_index (no_guess)
		p_list = self.get_list ('name', self.data, scalor_p + ".('", "')") 
		p_rhs = self.search_and_replace(self.get_list ('init-value'), p_list , self.names )
		self.pattern_write (mfio, [p_list,p_rhs], ["", " = ", ";"], 
									self.data, mIndex)
		#rhs_list = self.search_and_replace (rhs_list, rhs_names + ["t(iter)"], self.names + ["t"], False)
		
		# store variables that are printed as ode type or var type (for time sake)
		index_ode = self.new_index(self.index['ode']['yes'], self.index['print']['yes'])
		self.pattern_write (mfio, ["name", "comment"], [ info_p + ".('", "') = 'ode'; % ", ";"], self.data, index_ode)

		index_var = self.new_index(self.index['ode']['no'], self.index['print']['yes'], mIndex)
		self.pattern_write (mfio, ["name", "comment"], [ info_p + ".('", "') = 'var'; % ", ""], self.data, index_var)
		
		# close file
		self.close_file (mfio)
		
		
		
	"""
	CREATE PARAMETER ESTIMATE MATLAB FILE
	For parameter estimates, put into a separate folder
	with required files. To run, use the actual run.m file
	"""		
	def p_est_main (self, param_to_est):
		
		s_direct = self.directory + "/" + self.title + "/" + "guess";
		print(s_direct)
		if not os.path.exists (s_direct):
			os.makedirs (s_direct)
		else:
			print("WARNING: over writing folder " + self.title)
		
		# check if param_to_est is contained in value names
		temp = set(param_to_est).difference(set(self.names))
		if temp:
			error("These values do not exist for parameter estimate " + str(temp))
		n_guess_names = list(set(self.names).difference(set(param_to_est)))
		y_guess_names = list(param_to_est)
		n_guess = self.get_index(n_guess_names)
		y_guess = self.get_index(y_guess_names)
		#mIndex = self.new_index (self.index['ode']['no'], self.index['print']['yes'])
			
		
		#### DEFINE MAIN EST FILE
		#mfileo = "guess"
		#mfile = mfileo
		mfile = "guess"
		mfio = self.title + "/" + "guess" + "/" + "guess"
		mcomment = self.write_comment
		mwrite = self.write
		mp_est = "b"
		mp_err = "S"
		m_ode = "dydt"
		
		# open file
		self.open_file (mfio)
		
		# start function and create header
		mwrite (mfio, "function S = " + mfile + " (b, tdata, ydata, dydt, p, p_store)","\n\n","")
						#"if (nargout) == 2",
						#"\tout = " + m_ode,
						#"else"],"\n\n")
		mcomment (mfio, "numerical integration setup")
		mwrite (mfio, "tspan = 0:0.02:max(tdata);")
		mwrite (mfio, "options = odeset('InitialStep',0.02,'MaxStep',25, 'RelTol', 1e-9, 'AbsTol', 1e-9);")		
		mwrite (mfio, "[t,ysol] = ode15s(@(tt,y)ode_guess(tt,y,p,p_store,b),tspan,dydt,options);")
		
		mcomment (mfio, "find predicted values","\n","\n\n")
		mwrite (mfio, "ypred = interp1 (t,ysol,tdata);")
		#plot_names_index = self.get_index(["V_m", "Ca_i", "K_i", "Na_i", "Cl_i", "R_01"])
		mwrite (mfio, "figure(1)\n")
		cntr = 1
		for i in ["V_m", "Ca_i", "Na_i", "K_i", "Cl_i", "R_01"]:
			m = self.names.index(i);
			k = self.data_vec[m]['name'][self.data_vec[m]['name'].find('(')+1:self.data_vec[m]['name'].find(')')]
			mwrite (mfio, ["subplot(3,2," + str(cntr) + ")\nplot(tdata,ydata(:," + str(k) + "),'x','MarkerSize',10)",
							"hold on",
							"plot(t,ysol(:," + str(k) + "))",
							"hold off"])
			cntr = cntr + 1
		mwrite (mfio, "drawnow")
						#"plot(tdata,ydata(:,end),'x','MarkerSize',10)",
						#"hold on",
						#"plot(t,ysol(:,end))",
						#"hold off",
						#"drawnow"])
		
		mcomment (mfio, "compute total error","\n","\n\n")
		mwrite (mfio, ["S = 0;",
						"for i = 1:length(tdata)",
						"\tS = S +sum((ypred(i,:)-ydata(i,:)).^2);",
						"end"])
			
					
		mwrite (mfio, "end","","\n\n")	
		#mcomment (mfio, "Define sample data to fit model here","\n","\n\n")	
		#mcomment (mfio, "Get initial conditions here","\n","\n\n")	
		#mcomment (mfio, "Initial guess of parameters here","\n","\n\n")		
		#mcomment (mfio, "Minimization step","\n","\n\n")
		#mcomment (mfio, "output estimated parameters","\n","\n\n")
						
		self.close_file (mfio)
		
		### DEFINE PARAM FILE
		self.multi_paramfile ("guess/", "_g",n_guess, y_guess)
		
		### DEFINE ODE FILE
		self.multi_odefile ("guess/", "_g", n_guess, y_guess)
		
		### DEFINE RUN FILE
		self.multi_runfile ("guess/", "_g", n_guess, y_guess)
		
		### IC file, based on ODE and initial condition (given steady state is in )	
		create_file_name = self.title + "/guess/" + "create_ic"	
		self.open_file(create_file_name)
		mwrite (create_file_name, "function [tdata, ydata] = create_ic (t_end)","\n\n")
		mcomment (create_file_name, "create default time and data (assumed to be at steady state)")
		mwrite (create_file_name, "tdata = 0:0.5:t_end;")
		mIndex = self.order_index(self.info['range'])
		self.pattern_write (create_file_name, ["name", "init-value"], [ " ", " = ", ";"], self.data, mIndex,"")
		mwrite (create_file_name, "\n")
		p_list = self.get_list ("name")
		cntr = 1
		for it in self.index['ode']['yes']:
			p_list[it] = "ydata(:," + str(cntr) + ")"
			cntr = cntr + 1
		self.pattern_write (create_file_name, [p_list, 'init-value'],[" ", " = ", ".*ones(size(tdata)); "], self.data, self.index['ode']['yes'])

		mwrite (create_file_name, "\n")
		self.close_file(create_file_name)
		





			

		
	""""""	
	""""""	
	""""""	
	""""""
	""""""	
	""""""	
	""""""	
	""""""	
	""""""	
	""""""	
	""""""	
	""""""			
	""""""	
	"""			
	CREATE SLIDER_RUN_FILE MATLAB FILE FUNCTION
	contains the actual set-up and organisation of cell system
	in a structured order, is also main file
	"""
	def multi_runfile_slider (self, newSubDirectory = "", file_ext = "", no_guess = None, yes_guess = None):

		#mfile = "run" + file_ext
		mfile =  "run_slider" + file_ext
		mfio = self.title + "/" + newSubDirectory + "run_slider" + file_ext
		mcomment = self.write_comment
		mwrite = self.write
		m_print = "p_store"	# VC_NAME_PRINT
		m_noprint = "p"	# VEC_NAME
		m_ode = "dydt"	# ode var ODE_NAME
		m_out = "out"
		m_otime= 'out_time'
		t_end = "t_end"
		set_p = "initial"
		ode_method = "ode15s"
		m_ode_fname = "ode"
		m_add_events = "passEvents"
		m_struct = "options"
		# open file
		self.open_file (mfio)
		
		# define header
		if yes_guess is None: mp = m_out
		else: mp = "b_out"
		mwrite (mfio, "function [" + mp + "," + set_p + "] = " + mfile + " (" \
							+ set_p + "," + m_add_events + "," + m_struct \
							+ ")\n\n")
		#					+ set_p + ", " + t_end + ", b_in)", "\n\n")	
		
		#
		# persistent variables
		mwrite (mfio, "persistent name_only type_only iter_only","\n\n")
		
		
		# set options as nothing here
		
		#
		# NOTE: this is the extra commands for slider version of run
		#
		mwrite (mfio,"if nargout==2","\n")
				
		mwrite (mfio,"\tout = {};","\n")		
		# initial params : get initial params only (for slider)
		mcomment (mfio, "\tinitialize the original parameter output", "\n",)
		mwrite (mfio, '\t' + set_p + ' = get_ic;')

		#if no_guess	is None: 	mIndex = self.order_index(self.index['userdef']['yes']) #self.index['order']['init-value'];
		#else:					mIndex = self.order_index (no_guess)
		
		mIndex = self.order_index(self.index['userdef']['yes'])
		mIndex_odeyes = self.order_index(mIndex,self.index['ode']['yes'])
		mIndex_odeno = self.order_index(mIndex,self.index['ode']['no'])
		#for elem in yes_guess:
		#	mName[elem] = "b(" + str(cntr) + ")"
		#	cntr = cntr + 1
		
		##p_list = self.get_list ('name', self.data, set_p + ".","")
		p_list = self.get_list ('name')
		for elem in mIndex_odeyes:
			p_list[elem] = set_p + p_list[elem]
		for elem in mIndex_odeno:
			p_list[elem] = set_p + p_list[elem]
		
		p_rhs = self.search_and_replace(self.get_list ('init-value'), p_list , self.names )
		#self.pattern_write (mfio, [p_list,p_rhs], ["", " = ", ";"], 
		#							self.data, mIndex)
		
		# print info not used in this file. is a waste for in house slider [?]
		# store variables that are printed as ode type or var type (for time sake)
		#index_ode = self.new_index(self.index['ode']['yes'], self.index['print']['yes'])
		#self.pattern_write (mfio, ["name", "comment"], [ info_p + ".('", "') = 'ode'; % ", ";"], self.data, index_ode)
		#index_var = self.new_index(self.index['ode']['no'], self.index['print']['yes'], mIndex)
		#self.pattern_write (mfio, ["name", "comment"], [ info_p + ".('", "') = 'var'; % ", ""], self.data, index_var)
		
		mwrite (mfio, "else", "\n\n")
		#
		#
		#
		#print "wooooooo"
		#print str(len(self.index['ode']['yes']))
		#print str(len(self.new_index(self.index['rhs']['no'])))
		#print str(len(self.new_index(self.index['ode']['no'],self.index['rhs']['yes'])))
		
		## preallocate storage
		#mI_user_no_ode = self.new_index(self.index['print']['no'],self.index['ode']['no'])
		#mI_uder_no_ode_but_print = self.new_index(self.index['print']['yes'],self.index['ode']['no'])
		#mcomment (mfio, "preallocate storage for data")
		#mwrite (mfio, [m_noprint + " = zeros(1," + str(len(mI_user_no_ode)) + ");",
		#				m_print + " = zeros(1," + str(len(mI_uder_no_ode_but_print)+1) + ");",
		#				m_ode + " = zeros(1," + str(len(self.index['ode']['yes'])) +");"])

		# if parameterizing it, create initial data from external function
		# and define parameter initial values also
		if no_guess is not None:
			mcomment (mfio, "Define sample data to fit model here ()","\n","\n")
			mwrite (mfio, "[tdata,ydata] = create_ic (t_end);", "\n\n", "")
			mwrite (mfio, "\nif nargin<3")
			mName = [elem for elem in self.names]; cntr = 1;
			for elem in yes_guess:
				mName[elem] = "b(" + str(cntr) + ")"
				cntr = cntr + 1
			self.pattern_write (mfio, [mName, "init-value", "comment"], [ " ", " = ", "; \t% ", ""], self.data_vec, yes_guess)
			mwrite (mfio, "else\n\tb=b_in\nend")

									
		# ..user defined args
		
		##struct_vec_names = self.get_list ("name", self.data_vec)
		##t_end_vec = struct_vec_names[self.names.index(t_end)]
		#struct_vec_names[index_tend] = t_end
		#struct_vec_names = "name"
		
		#struct_names = self.new_names(self.names, self.index['userdef']['yes'], [set_p + ".",""])
		#self.pattern_write (mfio, ["name", p_list, "comment"], [ " ", " = ", "; % ", ""], self.data_vec, self.index['userdef']['yes'])	
		## .. function dependent
		#mIndex = self.order_index(self.index['userdef']['no'])
		#self.pattern_write (mfio, ["name", "init-value", "comment"], [ " ", " = ", "; % ", ""], self.data_vec, mIndex)
		# this is for the time output
		#mwrite (mfio, [m_print + "(" + str(len(mI_uder_no_ode_but_print) + 1) + ") = 0;" ])

		# initialise data for system run
		mcomment (mfio, "initialise the data prior to ode solver run")
		mwrite (mfio,"[" + m_print + "," + m_noprint + "," + m_ode + "] = initialise_data(" + set_p + ");\n")		

		#A = self.get_list ("name");
		#Avec = self.get_list ("name", self.data_vec);
		##index_debug = self.new_index(self.index['ode']['no'],self.index['print']['no']) + \
		#index_debug = self.new_index(self.index['print']['yes']) #+ \
		#	#self.new_index(self.index['ode']['yes']);
		
		#for elem in index_debug:
		#	print A[elem] + ":" + Avec[elem]
			
				
		#vec_names = []
		#for elem in self.data_vec:
		#	vec_names.append(elem['name']) # get new vec names list
		#vec_names = [elem['name'] for elem in self.data_vec]
		
		#index_only_p = self.new_index(self.index['ode']['no'], self.index['userdef']['yes'])
		#iter_only_p = self.new_names_restricted(vec_names,index_only_p)
		#names_only_p = self.new_names_restricted(self.names,index_only_p)
		#print names_only_p
		#for i in range(len(names_only_p)):
		#	iter_only_p[i] = int(re.search('\((.*)\)', iter_only_p[i]).group(1))

		#index_only_dydt = self.new_index (self.index['ode']['yes'])
		#iter_only_dydt = self.new_names_restricted(vec_names, self.index['ode']['yes'])
		#names_only_dydt = self.new_names_restricted(self.names, self.index['ode']['yes'])
		#for i in range(len(names_only_dydt)):
		#	iter_only_dydt[i] = int(re.search('\((.*)\)', iter_only_dydt[i]).group(1))		
		#
		#mwrite (mfio, "\nnum_p = [" , "", "")
		#for i in range(len(names_only_dydt)-1):
		#	mwrite (mfio, "," + names_only_
		
		iter_only = [re.search('\((.*)\)', elem['name']).group(1) for elem in self.data_vec]
		type_only = [elem[:elem.find('(')] for elem in [elem['name'] for elem in self.data_vec]]		
		# print the names
		mwrite (mfio, "\nif isempty(name_only)\n\tname_only = {" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, "'" + self.names[i] + "'" + ",","")
		mwrite (mfio, "'" +  self.names[-1] + "'};\nend","","")
		# print the types
		mwrite (mfio, "\nif isempty(type_only)\n\ttype_only = {" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, "'" + type_only[i] + "'"  + ",","")
		mwrite (mfio, "'" + type_only[-1] + "'};\nend","","")
		# print the iters
		mwrite (mfio, "\nif isempty(iter_only)\n\titer_only = [" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, iter_only[i] + ",","")
		mwrite (mfio, iter_only[-1] + "];\nend")
				# myString[myString.find("!")+1:myString.find("@")]
					
# s[len(start):-len(end)]

		# predefine time, if no input given
		# mwrite (mfio, ["if nargin<2", "\t" + t_end + " = " + set_p + "." + t_end + ";", "end\n"])		
		mwrite (mfio, [t_end + " = " + set_p + "." + t_end + ";\n"])		
		# add events (m_add_events is passEvents)
		m_size_add_events = "size_" + m_add_events
		m_add_events_vec = "time_vec"
		mcomment (mfio, "set up static events module if not defined","\n")
		mwrite (mfio, "if nargin==1\n" + \
				"\t" + m_add_events + ".time = []; % time point step all < t_end [vector]\n" \
				"\t" + m_add_events + ".change = {[]}; % changes to parameter values [vector]\n" \
				"elseif nargin==2\n" \
				"\t" + m_add_events + ".change = {[]," + m_add_events + ".change{:}};\n" \
				"end\n\n")
		mwrite (mfio, m_add_events_vec + " = [0," + m_add_events + ".time,t_end];")
		
		# start ode
		m_print_out = m_print + "_all"
		if no_guess is None:
			mcomment (mfio, "solve ode set through event set","\n","\n")
			mwrite (mfio, "t=[];\n" + m_print_out + "=" + m_print + ";")
			mwrite (mfio, "options_ode = odeset('InitialStep',0.02,'MaxStep',250, 'RelTol', 1e-9, 'AbsTol', 1e-9);","\n")
			mwrite (mfio, "for i=1:numel(" + m_add_events_vec + ")-1\n","","")			
			mcomment (mfio, "change parameter","\n","\n\t")
			mwrite (mfio, "\tCHANGICALS = " + m_add_events + ".change{i};\n" \
					"    for j=1:numel(CHANGICALS)\n" \
					"        niter = strmatch(CHANGICALS{j}{1},name_only,'exact');\n" \
					"        switch type_only{niter}\n" \
					"            case '" + m_noprint + "'\n" \
					"                " + m_noprint + "(iter_only(niter)) = CHANGICALS{j}{2};\n" 
					"            case '" + m_ode + "'\n" \
					"                " + m_ode + "(end,iter_only(niter)) = CHANGICALS{j}{2};\n" 
					"            case '" + m_print + "'\n" \
					"                " + m_print + "(end,iter_only(niter)) = CHANGICALS{j}{2};\n" 
					"        end\n" \
					"    end\n" \
					,"\n\n","")
					#"\t" + 	"for j=1:numel(CHANGICALS)\n" \
					#"\t\t" +	"if strcmp(CHANGICALS{j}{2},'dydt')\n" \
					#"\t\t\t" +		"dydt{dydt_iter{CHANGICALS{j}{1}}} = CHANGICALS{j}{2});\n" \
					#"\t\t" + 	"elseif strcmp(CHANGICALS{j}{2},'p')\n" \
					#"\t\t\t" +		"p{p_iter{CHANGICALS{j}{1}}} = CHANGICALS{j}{2});\n" \
					#"\t\t" + 	"end\n" \
					#"\t" +	"end" \
					#	,"\n\n","")
			mwrite (mfio, "\t[t_temp," + m_ode + "_temp] = " + ode_method + "(@(tt,y)" + m_ode_fname + "(tt,y," 
						+ m_noprint + "," + m_print_out + "(end,:)),...\n\t\t"
						+ m_add_events_vec + "(i:i+1)'.*1000," + m_ode + "(end,:),options_ode);")
			mwrite (mfio, "\t[~," + m_print + "] = " + m_ode_fname + "([],[],[],[]);")
			mwrite (mfio, "\t[~, index_max_t_] = max(" + m_print + "(:,end));")
			mwrite (mfio, "\t" + m_print_out + " = [" + m_print_out + ";" + m_print + "(1:index_max_t_,:)]; %#ok<AGROW>")
			mwrite (mfio, "\t" + m_ode + " = [" + m_ode + ";" + m_ode + "_temp(1:end-1,:)]; %#ok<AGROW>")
			mwrite (mfio, "\tt = [t(1:end-1) t_temp'];")
			mwrite (mfio, "end","\n\n")
		#	mwrite (mfio, "t=[t',t_end*1000];","\n\n")
		#	mwrite (mfio, ["options = odeset('InitialStep',0.2,'MaxStep',500, 'RelTol', 1e-9, 'AbsTol', 1e-9);",
		#					"[t," + m_ode + "] = " + ode_method + "(@(tt,y)" + m_ode_fname + "(tt,y," 
		#							+ m_noprint + "," + m_print + ")," +
		#						"[0," + t_end + "*1000]," + m_ode + ",options);"],"\n","\n\n")
			# output resultsinto a struct
			
			# return organised data for output to slider		
			mwrite (mfio, m_out + "= return_organised(" + m_ode + "," + m_print + "_all);")
			mwrite (mfio, m_out + ".var = " + m_print + "_all(:,end)./1000;")
			mwrite (mfio, m_out + ".ode = t./1000; % converted to seconds") #,	#struct;",	
		
		
		
		
		mwrite(mfio,'end')
		self.close_file (mfio)
		
		self.supp_input_ic()
		self.vectorised_to_readable()	# return file to simplify for plotting
		self.supp_initialise_data()




	""""""	
	""""""	
	""" SUPPLEMENTARY FUNCTION """	
	""""""	
	""""""	
	def supp_initialise_data (self):
		mfile =  "initialise_data"
		mfio = self.title + "/" + mfile
		mcomment = self.write_comment
		mwrite = self.write
		#m_initial_ode = "initial.ode"
		m_initial_param = "init_param"
		m_print = "p_store"
		m_noprint = "p"
		m_ode = "dydt"
		
		# open file
		self.open_file (mfio)

		# initialise function
		mwrite (mfio, "function [" + m_print + "," + m_noprint + "," + m_ode + "] = " + mfile + "(" + m_initial_param + ")\n\n")
											
		# user defined indexes
		mIndex = self.order_index(self.index['userdef']['yes'])				# user defined ALL
		mIndex_odeyes = self.order_index(mIndex,self.index['ode']['yes'])	# user defined ODE
		mIndex_odeno = self.order_index(mIndex,self.index['ode']['no'])		# user defined VAR
		
		# user defined name list
		p_list = self.get_list ('name')
		for elem in mIndex_odeyes: p_list[elem] = m_initial_param + "." + p_list[elem]
		for elem in mIndex_odeno: p_list[elem] = m_initial_param + "." + p_list[elem]
		
		
		# preallocate storage
		mI_user_no_ode = self.new_index(self.index['print']['no'],self.index['ode']['no'])
		mI_uder_no_ode_but_print = self.new_index(self.index['print']['yes'],self.index['ode']['no'])
		
		mwrite (mfio, "if nargin>0")
		mcomment (mfio, "preallocate storage for data")
		mwrite (mfio, [m_noprint + " = zeros(1," + str(len(mI_user_no_ode)) + ");",
						m_print + " = zeros(1," + str(len(mI_uder_no_ode_but_print)+1) + ");",
						m_ode + " = zeros(1," + str(len(self.index['ode']['yes'])) +");"])
		# search and replace initial values with new p_list
		self.pattern_write (mfio, ["name", p_list, "comment"], [ " ", " = ", "; % ", ""], self.data_vec, self.index['userdef']['yes'])	
		# .. function dependent
		mIndex = self.order_index(self.index['userdef']['no'])
		self.pattern_write (mfio, ["name", "init-value", "comment"], [ " ", " = ", "; % ", ""], self.data_vec, mIndex)
		mwrite (mfio, [m_print + "(" + str(len(mI_uder_no_ode_but_print) + 1) + ") = 0;" ])
		
		mwrite (mfio, "else\n% get properties of parameters instead")
		
		# if we want to return the initial condition, nargout = 1
		mwrite (mfio, "if nargout == 1")
		
		#mwrite (mfio, "p = struct;")		
		# user defined indexes
		mIndex_yes = self.order_index(self.index['userdef']['yes'])				# user defined ALL
		mIndex_odeyes = self.order_index(mIndex_yes,self.index['ode']['yes'])	# user defined ODE
		mIndex_odeno = self.order_index(mIndex_yes,self.index['ode']['no'])		# user defined VAR
		
		# user defined name list
		p_list = self.get_list ('name')
		for elem in mIndex_odeyes: p_list[elem] = m_initial_param + "." + p_list[elem]
		for elem in mIndex_odeno: p_list[elem] = m_initial_param + "." + p_list[elem]
		
		# search and replace initial values with new p_list
		p_rhs = self.search_and_replace(self.get_list ('init-value'), p_list , self.names )
		self.pattern_write (mfio, [p_list,p_rhs], ["", " = ", ";"],self.data, mIndex_yes)
		mwrite (mfio, m_print + " = " + m_initial_param + ";")		

		
		mwrite (mfio, "else")
		# add parameter properties to data (as p_store)
		iter_only = [re.search('\((.*)\)', elem['name']).group(1) for elem in self.data_vec]
		type_only = [elem[:elem.find('(')] for elem in [elem['name'] for elem in self.data_vec]]		
		# print the name_only
		mwrite (mfio, m_print + " = {" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, "'" + self.names[i] + "'" + ",","")
		mwrite (mfio, "'" +  self.names[-1] + "'};\n","","")
		# print the type_only (as p)
		mwrite (mfio, m_noprint + " = {" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, "'" + type_only[i] + "'"  + ",","")
		mwrite (mfio, "'" + type_only[-1] + "'};\n","","")
		# print the iter_only (as dydt)
		mwrite (mfio, m_ode + " = [" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, iter_only[i] + ",","")
		mwrite (mfio, iter_only[-1] + "];\n")
		
		mwrite (mfio, "end")		
				
				
		mwrite (mfio, "end")										
		
		# close file
		self.close_file (mfio)


	""""""	
	""""""	
	""" SUPPLEMENTARY FUNCTION """	
	""""""	
	""""""	
	def supp_input_ic (self):
		mfile =  "get_ic"
		mfio = self.title + "/" + mfile
		mcomment = self.write_comment
		mwrite = self.write
		m_initial_ode = "ode"
		m_initial_param = "initparam"

		# open file
		self.open_file (mfio)

		# initialise function
		mwrite (mfio, "function " + m_initial_param + " = " + mfile + "()\n\n")
											
		# user defined indexes
		mIndex = self.order_index(self.index['userdef']['yes'])				# user defined ALL
		mIndex_odeyes = self.order_index(mIndex,self.index['ode']['yes'])	# user defined ODE
		mIndex_odeno = self.order_index(mIndex,self.index['ode']['no'])		# user defined VAR
		
		# user defined name list
		p_list = self.get_list ('name')
		for elem in mIndex_odeyes: p_list[elem] = m_initial_param + "." + p_list[elem]
		for elem in mIndex_odeno: p_list[elem] = m_initial_param + "." + p_list[elem]
		
		# search and replace initial values with new p_list
		p_rhs = self.search_and_replace(self.get_list ('init-value'), p_list , self.names )
		self.pattern_write (mfio, [p_list,p_rhs], ["", " = ", ";"], 
									self.data, mIndex)
										
		# close file
		self.close_file (mfio)
		
		
		
				


	"""
	"vectorised_to_readable"
	
	"""
	def vectorised_to_readable (self):

		#mfile = "run" + file_ext
		mfile =  "return_organised"
		mfio = self.title + "/" + mfile
		mcomment = self.write_comment
		mwrite = self.write
		m_print = "p_store"	# VC_NAME_PRINT
		m_noprint = "p"	# VEC_NAME
		m_ode = "dydt"	# ode var ODE_NAME
		m_out = "out"
		m_otime= 'out_time'
		t_end = "t_end"
		set_p = "set_params"
		ode_method = "ode15s"
		m_ode_fname = "ode"
		m_add_events = "passEvents"
		m_struct = "options"


		# open file
		self.open_file (mfio)
		
		# set function
		mwrite (mfio, "function " + m_out + " = " + mfile + " (" \
							+ m_ode + "," + m_print + ")\n\n")

		
		mwrite (mfio, [m_out + "= struct;"])
		
		lhs = self.get_list ("name", self.data_vec, "", "", None)	# get a column form dict list
		num_iter = [filter(str.isdigit, elem) for elem in lhs]		# get number from derived list
		sto_iter = [elem[:elem.index('(')] for elem in lhs]		# get the name of storage for each
		
		
		
		index_ode = self.new_index(self.index['ode']['yes'], self.index['print']['yes'])		
		self.pattern_write (mfio, ["name", sto_iter, num_iter, "comment"], [ m_out + ".", " = ", "(:,","); % ",""], self.data, index_ode)
		#mwrite (mfio, m_out + ".ode.('t') = t./1000;", "\n\n")
		
				
		index_var = self.new_index(self.index['ode']['no'], self.index['print']['yes'])
		self.pattern_write (mfio, ["name", sto_iter, num_iter, "comment"], [ m_out + ".", " = ", "(:,","); % ",""], self.data, index_var)		
		#mwrite (mfio, m_out + ".var.('t') = " + m_print + "(:,end)./1000;", "\n\n")
		
		
		"""
		# get ODE output
		index_ode = self.new_index(self.index['ode']['yes'], self.index['print']['yes'])		
		#self.pattern_write (mfio, ["name", sto_iter, num_iter], [ m_out + ".", " = ", "(:,",");"], self.data, index_var)	
		"""
		
		mwrite(mfio,'end')
		
		# close file
		self.close_file (mfio)
		
		
		
		
		
		
		
		
	"""			
	CREATE FSA_cvodes MATLAB FILE
	contains the actual set-up and organisation of cell system
	in a structured order, is also main file
	"""
	def FSA_cvodes (self):

		mfile =  "FSA_cvodes"
		mfio = self.title + "/" + mfile
		mcomment = self.write_comment
		mwrite = self.write
		m_time = "t"
		m_time_temp = "t_temp"
		m_dtout = "20" # 5 ms default
		m_ode = "dydt"
		m_ode_temp = "dydt_temp"
		m_outsens = "yS"
		m_outsens_temp = "yS_temp"
		m_sens_names = "pS"
		m_data = "pdata"
		m_print = "p_store"	# VC_NAME_PRINT
		m_noprint = "p"	# VEC_NAME
		#m_ode = "dydt"	# ode var ODE_NAME
		#m_out = "out"
		#m_otime= 'out_time'
		t_end = "t_end"
		set_p = "initial"
		#ode_method = "ode15s"
		#m_ode_fname = "ode"
		m_add_events = "passEvents"
		m_struct = "options"
		f_jacobian = "djacfn"
		f_ode = "FSA_ode"
					
		# open file
		self.open_file (mfio)
		
		# define header
		mwrite (mfio, "function [" +  m_time + "," + m_ode + "_out," + m_outsens + "_out," + set_p + "] = " \
					+ mfile +  " (" + set_p + "," + m_sens_names + "," + m_add_events + "," + m_struct + ")")
		mcomment (mfio, "note: This function does not store the variables p_store")
		mcomment (mfio, "which is used in run_slider example for extracting temporal")
		mcomment (mfio, "activity of stored functions","\n\n")
							
		# user data
		mcomment (mfio, "define the initial parameters for both parameters and ODE")
		
		# persistent variables
		mwrite (mfio, "persistent name_only type_only iter_only","\n\n")
		#
		# NOTE: this is the extra commands for slider version of run
		#
		mwrite (mfio,"if nargout==4","\n")
		mwrite (mfio,m_time + " = {}; " + m_ode + "_out = {}; " + m_outsens + "_out = {}; ","\n")	
		mwrite (mfio,set_p + ".param = [];")		
		mwrite (mfio,set_p + ".ode = [];")		
	
		mcomment (mfio, "initialize the original parameter output", "\n",)
	
		mIndex = self.order_index(self.index['userdef']['yes']) 	
		mIndex_odeyes = self.order_index(self.index['userdef']['yes'],self.index['ode']['yes'])
		mIndex_odeno = self.order_index(self.index['userdef']['yes'],self.index['ode']['no'])
		
		p_list = self.get_list ('name', self.data, set_p + ".","")
		p_list = self.get_list ('name')
		for elem in mIndex_odeyes:
			p_list[elem] = set_p + "." + p_list[elem]
		for elem in mIndex_odeno:
			p_list[elem] = set_p + "." + p_list[elem]
		p_rhs = self.search_and_replace(self.get_list ('init-value'), p_list , self.names )
		#self.pattern_write (mfio, [p_list,p_rhs], ["", " = ", ";"], 
		#							self.data, mIndex)
		mwrite (mfio, set_p + " = initialise_data;")
									
		mwrite (mfio, "else", "\n\n")

		# preallocate storage
		mI_user_no_ode = self.new_index(self.index['print']['no'],self.index['ode']['no'])
		mI_uder_no_ode_but_print = self.new_index(self.index['print']['yes'],self.index['ode']['no'])
		mcomment (mfio, "preallocate storage for data")
		mwrite (mfio, [m_noprint + " = zeros(" + str(len(mI_user_no_ode)) + ",1);",
						m_print + " = zeros(" + str(len(mI_uder_no_ode_but_print)) + ",1);",
						m_ode + " = zeros(" + str(len(self.index['ode']['yes'])) + ",1);"])
		
		####### redefine RHS
		# order concateq_unordered according to self.order.init_value
		#concateq = []
		#for elem in self.index['order']['init-value']:
		#	if self.names[elem] in concateq_unordered:
		#		concateq.append(self.names[elem])
		#print concateq
		concateq = self.new_names_restricted(self.names, mI_uder_no_ode_but_print)
		concat_index = [self.names.index(elem) for elem in concateq]
		#a = self.get_list('rhs') #[x['rhs'] for x in self.data]
		a = self.get_list ("rhs", self.data,"","", concat_index)
		#for elem in concat_index:
		#	a[elem] = ''
		#print a
			
 		#concateq_rhs_sub = []
 		#for elem in concateq:
 		#	concateq_rhs_sub.append('(' + self.data[self.names.index(elem)]['rhs'] + ')')
		concateq_rhs_sub = ['(' + elem + ')' for elem in a]
		#print concateq_rhs_sub
		######## this is not finished



		self.pattern_write (mfio, ["name", p_list, "comment"], [ " ", " = ", "; % ", ""], self.data_vec, self.index['userdef']['yes'])	
		# .. function dependent
		mIndex = self.order_index(self.index['userdef']['no'])
		self.pattern_write (mfio, ["name", "init-value", "comment"], [ " ", " = ", "; % ", ""], self.data_vec, mIndex)
		# this is the time output
		#mwrite (mfio, [m_print + "(" + str(len(mI_uder_no_ode_but_print) + 1) + ") = 0;" ])

		iter_only = [re.search('\((.*)\)', elem['name']).group(1) for elem in self.data_vec]
		type_only = [elem[:elem.find('(')] for elem in [elem['name'] for elem in self.data_vec]]		
		# print the names
		mwrite (mfio, "\nif isempty(name_only)\n\tname_only = {" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, "'" + self.names[i] + "'" + ",","")
		mwrite (mfio, "'" +  self.names[-1] + "'};\nend","","")
		# print the types
		mwrite (mfio, "\nif isempty(type_only)\n\ttype_only = {" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, "'" + type_only[i] + "'"  + ",","")
		mwrite (mfio, "'" + type_only[-1] + "'};\nend","","")
		# print the iters
		mwrite (mfio, "\nif isempty(iter_only)\n\titer_only = [" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, iter_only[i] + ",","")
		mwrite (mfio, iter_only[-1] + "];\nend")
		
		mwrite (mfio, [t_end + " = " + set_p + "." + t_end + ";\n"])
		
		
		# Check if chosen variables (::m_sens_names)
		mcomment (mfio, "check if chosen sensitive variables exist")
		mwrite (mfio, ["for iname = 1:length(" + m_sens_names + ")",
						"\tif isempty(strmatch(" + m_sens_names + "(iname), name_only,'exact'))",
						"\t\terror(['The sensitivity variable ' " + m_sens_names + "{iname}, ' does not exist'])",
						"\tend",
						"end"]) 
			
						
		#m_ode = "dydt"
		#m_ode_new = "dydt_coll"
		#m_outsens = "yS"
		#m_outsens_new = "yS_coll"
		#m_outsens_init = "yS0"
		#m_sens_names = "pS"
		#m_data = "pdata"
		
					 			
		# add events (m_add_events is passEvents)
		m_size_add_events = "size_" + m_add_events
		m_add_events_vec = "time_vec"
		mcomment (mfio, "set up static events module if not defined","\n","\n\n")
		mwrite (mfio, "if nargin==2\n" + \
				"\t" + m_add_events + ".time = []; % time point step all < t_end [vector]\n" \
				"\t" + m_add_events + ".change = {[]}; % changes to parameter values [vector]\n" \
				"elseif nargin==3\n" \
				"\t" + m_add_events + ".change = {[]," + m_add_events + ".change{:}};\n" \
				"end\n\n")

		mcomment (mfio, "Time vector when changes are made (x1000 convert to ms form s)")				
		mwrite (mfio, m_add_events_vec + " = [0," + m_add_events + ".time," + t_end + "].*1000;")
		mwrite (mfio, "Ns = length(" + m_sens_names + ");")
		mwrite (mfio, "def_ParamList = zeros(size(" + m_sens_names + "));")
		mwrite (mfio, "def_ParamScales = def_ParamList;")
		mwrite (mfio, "for elem = 1:length(" + m_sens_names + ")\n" +
						"\t def_ParamList(elem) = " + "iter_only(strmatch(" + m_sens_names + "(elem), name_only, 'exact'));\n" +
						"\t def_ParamScales(elem) = " + m_noprint + "(def_ParamList(elem));\n" +
						"end")
		
		#pS_index = 
		#mwrite (mfio, m_ode + " = " + m_ode + ";")		
		#mwrite (mfio, m_time + " = " + m_time + ";")	
			
		mcomment (mfio, "inital data")
		mwrite (mfio, m_outsens + " = zeros(" + str(len(self.index['ode']['yes'])) + ",Ns);")
		mwrite (mfio, m_time + " = 0.0;")  				
		
		# start ode
		#m_print_out = m_print + "_all"
		mcomment (mfio, "solve ode set through event set","\n","\n")
		mwrite (mfio, "for i=1:numel(" + m_add_events_vec + ")-1\n","","")

		######
		# FOR LOOP STARTS
		######
		
		mcomment (mfio, "change parameter set","\n","\n\t")
		mwrite (mfio, "\tCHANGICALS = " + m_add_events + ".change{i};\n" \
					"    for j=1:numel(CHANGICALS)\n" \
					"        niter = strmatch(CHANGICALS{j}{1},name_only,'exact');\n" \
					"        switch type_only{niter}\n" \
					"            case '" + m_noprint + "'\n" \
					"                " + m_noprint + "(iter_only(niter)) = CHANGICALS{j}{2};\n" 
					"            case '" + m_ode + "'\n" \
					"                " + m_ode + "(end,iter_only(niter)) = CHANGICALS{j}{2};\n" 
					"            case '" + m_print + "'\n" \
					"                " + m_print + "(end,iter_only(niter)) = CHANGICALS{j}{2};\n" 
					"        end\n" \
					"    end\n" \
					,"\n\n","")
		
		mcomment (mfio, "User data structure")
		mwrite (mfio, [m_data + "." + m_print + " = " + m_print + ";",
						m_data + "." + m_noprint + " = " + m_noprint + ";"])
			
		mcomment (mfio, "CVODE initialisation!")
		mwrite (mfio, "options_cvode = CVodeSetOptions('UserData'," + m_data + ",...\n" +
						"\t\t'RelTol'" + "," + "1.0e-4" + ",...\n" +
						"\t\t'AbsTol'" + "," + "ones(" + str(len(self.index['ode']['yes'])),",1).*1e-12" + ",...\n" +
						"\t\t'JacobianFn'," + "@" + f_jacobian +
						");")
		
		mcomment (mfio, "Moniter output for order etc (use with option?)","\n","\n\n")	
			
			
		mcomment (mfio, "initialise cvode now...","\n","\n\n")	
		mwrite (mfio, "CVodeInit(@" + f_ode + "," + "'BDF', 'Newton'," + m_time + "(end), "
						 + m_ode + "(:,end), options_cvode);")
						

		mcomment (mfio, "FSA initialisation options...","\n","\n\n")	
		mwrite (mfio, "FSAoptions = CVodeSensSetOptions('method','Simultaneous',...\n" +
						"\t\t'ErrControl'" + "," + "true" + ",...\n" +
						"\t\t'ParamField'" + ",'" + m_noprint + "',...\n" +
						"\t\t'ParamList'" + "," + "def_ParamList" + ",...\n" +
						"\t\t'ParamScales'" + "," + "def_ParamScales" + "...\n" +
						"\t);")
		mwrite (mfio, "CVodeSensInit(Ns, [], " + m_outsens + "(:,:,end), FSAoptions);\n")
				


	
		# do the noval FSA analysis
		mcomment (mfio, "FSA return and ODE cvodes solve will be done here","\n\n","\n\t")
		mwrite (mfio, m_time_temp + " = " + "[" + m_add_events_vec + "(i)" + "+" + m_dtout + ":" + m_dtout + ":" +  m_add_events_vec + "(i+1)];")
		mwrite (mfio, "[status, ~," + m_ode_temp + ", " + m_outsens_temp + "]" +
					  " = CVode (" + m_time_temp + ", " + "'Normal');")
		mwrite (mfio, ["if status==1 % error occured during cvodes simulaiton progression",
						"\terror('error in parsed FSA file')",
						"end"])
		# m_time : total time vector as arg out
		# m
		#m_time_temp
		# m_add_events_vec
		
		mcomment (mfio, "update variables")
		mwrite (mfio, m_time + " = [" + m_time + "," + m_time_temp + "];")
		mwrite (mfio, m_ode + " = cat(2," + m_ode + "," + m_ode_temp + ");")
		mwrite (mfio, m_outsens + " = cat(3," + m_outsens + "," + m_outsens_temp + ");")

		mcomment (mfio, "recover param and function values")
		#mwrite (mfio, "[~, ~, ~, " + m_print + ", " + m_noprint + "] = " \
		#			+ f_ode + "(" + m_time + "," + m_ode + "(1,:), " + m_data + ");")
		
		mwrite (mfio, "end % end of events for loop", "\n\n")
			
		mwrite (mfio, "si = CVodeGetStats;")
		mwrite (mfio, "CVodeFree;")		
		
		
				
		# change variables for output
		mcomment (mfio, "rewrite variables for easy output")
		mwrite (mfio, m_outsens + "_out = struct;")
		mwrite (mfio, m_ode + "_out = struct;")
		ode_list = self.new_names_restricted(self.names,self.index['ode']['yes'])
		[mwrite (mfio, m_outsens + "_out." + elem + " = struct; ","") for elem in ode_list]
		mwrite(mfio,"")
		[mwrite (mfio, m_ode + "_out." + elem + " = " + m_ode + "(" + str(ode_list.index(elem)+1) + ",:); ","") for elem in ode_list]		
		#mwrite (mfio, "\t" + m_ode + "_out.(" + m_sens_names + "(elem)) = " + m_ode)
		
		mwrite (mfio, "\nfor elem = 1:length(" + m_sens_names + ")")
		#mwrite (mfio, "\tfor iter = 1:" + str(len(self.index['ode']['yes'])))
		[mwrite (mfio, "\ttemp = " + m_outsens + "(" + str(ode_list.index(elem)+1) + ",elem,:); " + m_outsens + "_out." + elem + ".(" + m_sens_names + "{elem}) = temp(:);") for elem in ode_list]
		#mwrite (mfio, "\tend")
		mwrite (mfio, "end")
		#type_only = [elem[:elem.find('(')] for elem in [elem['name'] for elem in self.data_vec]]		

		mwrite (mfio, "end % end of all session", "\n\n","\n")





		
		mwrite (mfio, "\n\nend % end of function", "\n\n")
		
		
		
		# Jacobian function
		mcomment (mfio, "numerical jacobian using the file numjac.m","\n","\n\n")
		mwrite (mfio, "function [J, flag, new_data] = " + f_jacobian + \
						"(" + m_time + "," + m_ode + ",fy," + m_data + ")")
		mwrite (mfio, "J = numjac(@(tt,yy)" + f_ode + "(tt,yy," + m_data + ")," \
						+ m_time + "," + m_ode + "," + f_ode + "(" + m_time + "," \
						+ m_ode + "," + m_data + "),...\n" \
						+ "\t\t\tones(size(" + m_ode + ")) .* 1e-8,[],0);")	
		mwrite (mfio, "flag = 0;")
		mwrite (mfio, "new_data = [];")
		mwrite (mfio, "end\n\n")
		
		
		
		#mcomment (mfio, "new ode function which shows no p_store independence","\n","\n\n")
		#mwrite (mfio, "function [] = " + f_ode + \
		#				"(" + m_time + "," + m_ode + "," + m_data + ")")		
		#mwrite (mfio, "end\n\n")

		
		# close file
		self.close_file (mfio)
		
		
		
		
		
		
	"""
	CREATE FSA sensitivity ODE-FILE MATLAB FILE FUNCTION
	the actual ODE system set up as a separate m file
	for passing into an ODE system
	"""	
	def FSA_ode (self):
		"""Create FSA matlab ODE file (ode.m) to be used with :py:mod:`multi_runfile` and :py:mod:`multi_runfile_slider` modules"""

						
		mfile = "FSA_ode"
		mfio = self.title + "/" + mfile
		mcomment = self.write_comment
		mwrite = self.write
		m_data = "pdata"
		m_print = "p_store"	# VC_NAME_PRINT
		m_noprint = "p"	# VEC_NAME
		m_ode = "dydt"	# ode var ODE_NAME
		m_ode_rhs = "y" # rhs dydt as y
		mpersistent_out = "out_" + m_print
		#mpersistent_temp = "pers_" + m_print
		len_of_odes = len(self.new_index(self.index['ode']['yes']))
		
		# open file
		self.open_file (mfio)

		# define header (function name etc)
		#if yes_guess is None: new_arg = ""; temp = ""; jo = "," + mpersistent_out
		#else: new_arg = "b"; temp = "," + new_arg; jo = ""
		mwrite (mfio, "function [" + m_ode + ", flag, new_data] = " + mfile + "(t," + m_ode_rhs + ", " + m_data + ")\n\n")
										

		#mwrite (mfio, "persistent " + m_print + " " + m_noprint, "\n\n", "\n")		
		
		#mwrite (mfio, m_print + " = " + m_data + "." + m_print + ";")
		#mwrite (mfio, m_noprint + " = " + m_data + "." + m_noprint + ";")
		
		#if yes_guess is None:								
		#mwrite (mfio, ["persistent " + mpersistent_temp + " counter",
							#"if isempty(" + mpersistent_temp + ")", 
							#"\t" + mpersistent_temp + " = zeros(200,size(" + m_print + ",2));\n\tcounter = 1;",
							#"end"], "\n" )
		#mwrite (mfio, ["if mod(counter-1,200) == 0",
							#"\t" + mpersistent_temp + " = [" + mpersistent_temp + "; zeros(200,size(" + m_print + ",2))];",
							#"end"], "\n")		
		
		#mcomment (mfio, "solve the ode step...")
		
		#mwrite (mfio, [m_ode + " = zeros(" + str(len_of_odes) + ",1);"])
		#mwrite (mfio, "if (nargout == 1)")

		#if yes_guess is None:
		# if guess
		oderhs_names = self.vectorise_name (self.new_index(self.index['var']['yes'],self.index['print']['yes']), "p_store")
		oderhs_names = self.vectorise_name (self.new_index(self.index['var']['yes'],self.index['print']['no']), "p", oderhs_names)
		oderhs_names = self.vectorise_name (self.index['ode']['yes'], "y", oderhs_names)
		m_rhs = self.search_and_replace (self.get_list ("rhs", self.data, "", "", None), oderhs_names, self.names)

			
		#mwrite (mfio, m_print + "(end) = t;")
		mwrite (mfio, m_ode + " = zeros(" + str(len(self.index['ode']['yes'])) + ",1);")
		mwrite (mfio, "flag = 0;")
		mwrite (mfio, "new_data = [];")		
		
		#mwrite (mfio, "\nif nargout < 4")
		mwrite (mfio, "\t" + m_print + " = " + m_data + "." + m_print + ";")
		mwrite (mfio, "\t" + m_noprint + " = " + m_data + "." + m_noprint + ";")
		mwrite (mfio, "")			
		setODEFormatAs = 2
		# 1: old format, based on function
		# 2: new format, ordering dependent
		if setODEFormatAs == 1:
		# ORIGINAL ODE/PARAM UPDATE
			for i in range(len(self.info['function'])):
				mcomment (mfio, self.info['function'][i]['name'],"\n","\n")
				self.pattern_write (mfio, ["name", m_rhs, "comment"], [ "\t", " = ", "; %", ""], self.data_vec, self.info['function'][i]['index'])
		elif setODEFormatAs == 2:
			index_nonode = self.new_dependancy_index(self.new_index(self.index['ode']['no'],self.index['rhs']['yes']),'rhs')
			#for indexer in index_nonode:
			#	print self.data[indexer]['rhs']
			# print all function related dependencies in front of order
			self.pattern_write (mfio, ["name", m_rhs, "comment"], [ "\t", " = ", "; %", ""], self.data_vec, index_nonode)
			# print all ODE related dependencies out to file in last order
			self.pattern_write (mfio, ["name", m_rhs, "comment"], [ "\t", " = ", "; %", ""], self.data_vec, self.index['ode']['yes'])			
			#new_dependancy_index	
			#mIndex_nonode = self.order_index (self.new_index(self.index['ode']['no'],self.index['rhs']['yes']))
			#mIndex_ode = self.order_index (self.index['ode']['yes'])
			#self.pattern_write (mfio, ["name", m_rhs, self.names, "comment"], [ " ", " = ", "; % {", "} ", " "], self.data_vec, mIndex_nonode)
			#self.pattern_write (mfio, ["name", m_rhs, self.names, "comment"], [ " ", " = ", "; % {", "} ", " "], self.data_vec, mIndex_ode)
		else:
			None		
		
		#mwrite (mfio, "else")
		#mwrite (mfio, "\ts_" + m_print + " = " + m_print + ";")
		#mwrite (mfio, "\ts_" + m_noprint + " = " + m_noprint + ";")
		#mwrite (mfio, "end")
			

			
		mwrite (mfio, "\nend % end ode function")		
				
		self.close_file (mfio)






	""""""	
	""""""	
	""""""	
	""""""
	""""""	
	""""""	
	""""""	
	""""""	
	""""""	
	""""""	
	""""""	
	""""""			
	""""""	
	"""			
	CREATE SLIDER_RUN_FILE MATLAB FILE FUNCTION *** WITH INSET ODE ***
	contains the actual set-up and organisation of cell system
	in a structured order, is also main file
	"""
	def inset_runfile_slider (self, newSubDirectory = "", file_ext = "", no_guess = None, yes_guess = None):

		#mfile = "run" + file_ext
		mfile =  "run_slider_nested" + file_ext
		mfio = self.title + "/" + newSubDirectory + mfile
		mcomment = self.write_comment
		mwrite = self.write
		m_print = "p_store"	# VC_NAME_PRINT
		m_noprint = "p"	# VEC_NAME
		m_ode = "dydt"	# ode var ODE_NAME
		m_out = "out"
		m_otime= 'out_time'
		t_end = "t_end"
		set_p = "initial"
		ode_method = "ode15s"
		m_ode_fname = "ode"
		m_add_events = "passEvents"
		m_struct = "options"
		m_print_out = m_print + "_out"

		# open file
		self.open_file (mfio)
		
		# define header
		if yes_guess is None: mp = m_out
		else: mp = "b_out"
		mwrite (mfio, "function [" + mp + "," + set_p + "] = " + mfile + " (" \
							+ set_p + "," + m_add_events + "," + m_struct \
							+ ")\n\n")
		#					+ set_p + ", " + t_end + ", b_in)", "\n\n")	
		
		#
		# persistent variables
		mwrite (mfio, ["persistent name_only type_only iter_only\n",
						"counter = 1;","\n\n"])
		mwrite (mfio, m_print_out + " = 0;")
		
		
		
		#---------------------# DEFINE NESTED ODE
						
		m_ode_nested = "dydt"	# ode var ODE_NAME
		m_ode_rhs = "y" # rhs dydt as y
		#mpersistent_out = "out_" + m_print
		mpersistent_temp = m_print + "_out"
		len_of_odes = len(self.new_index(self.index['ode']['yes']))
		
		
		mcomment (mfio, "function to solve the ode step...")
		mwrite (mfio, "function " + m_ode + " = odesystem (t,y)")	
		mwrite (mfio, ["if mod(counter-1,200) == 0",
			"\t" + mpersistent_temp + " = [" + mpersistent_temp + "; zeros(200,size(" + m_print + ",2))];","end"], "\n")
		mwrite (mfio, [m_ode_nested + " = zeros(" + str(len_of_odes) + ",1);"])
		#mwrite (mfio, "if (nargout == 1)")

		oderhs_names = self.vectorise_name (self.new_index(self.index['var']['yes'],self.index['print']['yes']), "p_store")
		oderhs_names = self.vectorise_name (self.new_index(self.index['var']['yes'],self.index['print']['no']), "p", oderhs_names)
		oderhs_names = self.vectorise_name (self.index['ode']['yes'], "y", oderhs_names)
		m_rhs = self.search_and_replace (self.get_list ("rhs", self.data, "", "", None), oderhs_names, self.names)

			
		mwrite (mfio, m_print + "(end) = t;")
		setODEFormatAs = 2
		# 1: old format, based on function
		# 2: new format, ordering dependent
		if setODEFormatAs == 1:
		# ORIGINAL ODE/PARAM UPDATE
			for i in range(len(self.info['function'])):
				mcomment (mfio, self.info['function'][i]['name'],"\n","\n")
				self.pattern_write (mfio, ["name", m_rhs, "comment"], [ " ", " = ", "; %", ""], self.data_vec, self.info['function'][i]['index'])
		elif setODEFormatAs == 2:
			index_nonode = self.new_dependancy_index(self.new_index(self.index['ode']['no'],self.index['rhs']['yes']),'rhs')
			#for indexer in index_nonode:
			#	print self.data[indexer]['rhs']
			# print all function related dependencies in front of order
			self.pattern_write (mfio, ["name", m_rhs, "comment"], [ " ", " = ", "; %", ""], self.data_vec, index_nonode)
			# print all ODE related dependencies out to file in last order
			self.pattern_write (mfio, ["name", m_rhs, "comment"], [ " ", " = ", "; %", ""], self.data_vec, self.index['ode']['yes'])			

		else:
			None		
		mwrite (mfio, m_print_out + "(counter,:)=" + m_print + ";")
		mwrite (mfio, "counter = counter + 1;")
		mwrite (mfio, "end \n%--- end ode function")		
				
		#---------------------#  DEFINE NESTED ODE


		
		#
		# NOTE: this is the extra commands for slider version of run
		#
		mwrite (mfio,"if nargout==2","\n")
				
		mwrite (mfio,"\tout = {};","\n")		
		# initial params : get initial params only (for slider)
		mcomment (mfio, "\tinitialize the original parameter output", "\n",)
		mwrite (mfio, '\t' + set_p + ' = initialise_data;')

		#if no_guess	is None: 	mIndex = self.order_index(self.index['userdef']['yes']) #self.index['order']['init-value'];
		#else:					mIndex = self.order_index (no_guess)
		
		mIndex = self.order_index(self.index['userdef']['yes'])
		mIndex_odeyes = self.order_index(mIndex,self.index['ode']['yes'])
		mIndex_odeno = self.order_index(mIndex,self.index['ode']['no'])
		#for elem in yes_guess:
		#	mName[elem] = "b(" + str(cntr) + ")"
		#	cntr = cntr + 1
		
		##p_list = self.get_list ('name', self.data, set_p + ".","")
		p_list = self.get_list ('name')
		for elem in mIndex_odeyes:
			p_list[elem] = set_p + p_list[elem]
		for elem in mIndex_odeno:
			p_list[elem] = set_p + p_list[elem]
		
		p_rhs = self.search_and_replace(self.get_list ('init-value'), p_list , self.names )

		mwrite (mfio, "else", "\n\n")

		# if parameterizing it, create initial data from external function
		# and define parameter initial values also
		if no_guess is not None:
			mcomment (mfio, "Define sample data to fit model here ()","\n","\n")
			mwrite (mfio, "[tdata,ydata] = create_ic (t_end);", "\n\n", "")
			mwrite (mfio, "\nif nargin<3")
			mName = [elem for elem in self.names]; cntr = 1;
			for elem in yes_guess:
				mName[elem] = "b(" + str(cntr) + ")"
				cntr = cntr + 1
			self.pattern_write (mfio, [mName, "init-value", "comment"], [ " ", " = ", "; \t% ", ""], self.data_vec, yes_guess)
			mwrite (mfio, "else\n\tb=b_in\nend")

									
		# ..user defined args
		

		# initialise data for system run
		mcomment (mfio, "initialise the data prior to ode solver run")
		mwrite (mfio,"[" + m_print + "," + m_noprint + "," + m_ode + "] = initialise_data(" + set_p + ");\n")		
		mwrite (mfio, m_print_out + "=" + m_print + ";")

		
		
		iter_only = [re.search('\((.*)\)', elem['name']).group(1) for elem in self.data_vec]
		type_only = [elem[:elem.find('(')] for elem in [elem['name'] for elem in self.data_vec]]		
		# print the names
		mwrite (mfio, "\nif isempty(name_only)\n\tname_only = {" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, "'" + self.names[i] + "'" + ",","")
		mwrite (mfio, "'" +  self.names[-1] + "'};\nend","","")
		# print the types
		mwrite (mfio, "\nif isempty(type_only)\n\ttype_only = {" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, "'" + type_only[i] + "'"  + ",","")
		mwrite (mfio, "'" + type_only[-1] + "'};\nend","","")
		# print the iters
		mwrite (mfio, "\nif isempty(iter_only)\n\titer_only = [" , "", "")
		for i in range(len(iter_only)-1):
			mwrite (mfio, iter_only[i] + ",","")
		mwrite (mfio, iter_only[-1] + "];\nend")
				# myString[myString.find("!")+1:myString.find("@")]
					
# s[len(start):-len(end)]

		# predefine time, if no input given
		# mwrite (mfio, ["if nargin<2", "\t" + t_end + " = " + set_p + "." + t_end + ";", "end\n"])		
		mwrite (mfio, [t_end + " = " + set_p + "." + t_end + ";\n"])		
		# add events (m_add_events is passEvents)
		m_size_add_events = "size_" + m_add_events
		m_add_events_vec = "time_vec"
		mcomment (mfio, "set up static events module if not defined","\n")
		mwrite (mfio, "if nargin==1\n" + \
				"\t" + m_add_events + ".time = []; % time point step all < t_end [vector]\n" \
				"\t" + m_add_events + ".change = {[]}; % changes to parameter values [vector]\n" \
				"elseif nargin==2\n" \
				"\t" + m_add_events + ".change = {[]," + m_add_events + ".change{:}};\n" \
				"end\n\n")
		mwrite (mfio, m_add_events_vec + " = [0," + m_add_events + ".time,t_end];")
		
		# start ode
		if no_guess is None:
			mcomment (mfio, "solve ode set through event set","\n","\n")
			mwrite (mfio, "t=[];")
			mwrite (mfio, "% " + m_print_out + "=" + m_print + ";")
			mwrite (mfio, "options_ode = odeset('InitialStep',0.02,'MaxStep',250, 'RelTol', 1e-9, 'AbsTol', 1e-9);","\n")
			mwrite (mfio, "for i=1:numel(" + m_add_events_vec + ")-1\n","","")			
			mcomment (mfio, "change parameter","\n","\n\t")
			mwrite (mfio, "\tCHANGICALS = " + m_add_events + ".change{i};\n" \
					"    for j=1:numel(CHANGICALS)\n" \
					"        niter = strmatch(CHANGICALS{j}{1},name_only,'exact');\n" \
					"        switch type_only{niter}\n" \
					"            case '" + m_noprint + "'\n" \
					"                " + m_noprint + "(iter_only(niter)) = CHANGICALS{j}{2};\n" 
					"            case '" + m_ode + "'\n" \
					"                " + m_ode + "(end,iter_only(niter)) = CHANGICALS{j}{2};\n" 
					"            case '" + m_print + "'\n" \
					"                " + m_print + "(end,iter_only(niter)) = CHANGICALS{j}{2};\n" 
					"        end\n" \
					"    end\n" \
					,"\n\n","")

			mwrite (mfio, "\t[t_temp," + m_ode + "_temp] = " + ode_method + "(@" + m_ode_fname + "system,...\n\t\t"
						+ m_add_events_vec + "(i:i+1)'.*1000," + m_ode + "(end,:),options_ode);")
			#mwrite (mfio, "%\t[~," + m_print + "] = " + m_ode_fname + "([],[],[],[]);")
			#mwrite (mfio, "%\t[~, index_max_t_] = max(" + m_print + "(:,end));")
			#mwrite (mfio, "%\t" + m_print_out + " = [" + m_print_out + ";" + m_print + "(1:index_max_t_,:)]; %#ok<AGROW>")
			mwrite (mfio, "\t" + m_ode + " = [" + m_ode + ";" + m_ode + "_temp(1:end-1,:)]; %#ok<AGROW>")
			mwrite (mfio, "\tt = [t(1:end-1) t_temp'];")
			mwrite (mfio, "end","\n\n")

			mwrite (mfio, "[~,index_max_t] = max(" + m_print_out + "(:,end));")			
			mwrite (mfio, m_out + "= return_organised(" + m_ode + "," + m_print_out + "(1:index_max_t,:));")
			mwrite (mfio, m_out + ".var = " + m_print_out + "(1:index_max_t,end)./1000;")
			mwrite (mfio, m_out + ".ode = t./1000; % converted to seconds") #,	#struct;",

		
		
		
		
		mwrite(mfio,'end')
		
		
		mcomment (mfio, "--- end of all function run_slider_nested")
		mwrite(mfio,'end')
		self.close_file (mfio)
		
		#self.supp_input_ic()
		#self.vectorised_to_readable()	# return file to simplify for plotting
		#self.supp_initialise_data()


##### MATLAB - 'ALL TOGETHER' NESTED FUNCTION RUN_SLIDER
# this causes simulations to run faster as functions are nested and available locally

	#def inset_runfile_slider (self


