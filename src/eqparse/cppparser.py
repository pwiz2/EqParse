#! /usr/bin/env python

from baseparse import *

"""
C++ PARSER
"""

class CppParser (BaseParse):
	def __init__ (self, Lib):
		super(CppParser, self).__init__()
		self.l_enclose = "["
		self.r_enclose = "]"
		self.operations = ['\nif(', ")\n", "else{\n\t", "}\nelse if{\n", "}\n", "!=", "pow(", ",", ")", "sqrt"]

		self.vec_counter_start = 0
		self.comment_prefix = "//"		
		self.fileExt = "hpp"		
		
		self.initialise_library (Lib)	

	def cpp_original (self):
		# data vars
		mfile = self.title + "_o"
		mcomment = self.write_comment
		mwrite = self.write
				
		# open file
		self.open_file (mfile)
		
		# define header
		mwrite (mfile, ["#ifndef " + (mfile.upper() + "." + self.fileExt.upper()).replace('.','_'),
							"#define " + (mfile.upper() + "." + self.fileExt.upper()).replace('.','_')])
		mwrite (mfile, "#include \"pasmc_proxy.h\"","\n\n")
		
		# start the structure of class {cell}
		mwrite (mfile, ["struct " + mfile.replace(' ','_') + " : public cell_u_vars {",
						mfile.replace(' ','_') + "(){};"])
		mwrite (mfile, ["std::string DIRNAME = " + '"' + self.title + '";'],"\n\n")
		# define variables
		mIndex = self.order_index (self.info['range'])
		self.pattern_write (mfile, ["name","init-value", "dimension", "comment"],
									["double ", " = ", ";\t// {", "} ", ""], self.data, mIndex)
									# order DOES matter for some reason
		# define functions
		fNames = [elem['name'] for elem in self.info['function']]
		fComment = [elem['comment'] for elem in self.info['function']]

		self.pattern_write (mfile, [fNames,fComment],
									["void ", "();\t//",""],self.data ,range(len(fNames)))
		mwrite (mfile, "void operator() ();")
		#mwrite (mfile, ["void fParameterize();",
		#				"void operator() ();",
		#				"typedef void (" + mfile.replace(' ','_') + "::*func_advance) (void);",
		#				"func_advance step_it {&" + mfile.replace(' ','_') + "::fODEsolveAndUpdate};"])
		#mwrite (mfile, ["void parameterise_it (bool choice) {step_it = choice ? &" + mfile +"::fParameterize :",
		#				"\t\t&" + mfile +"::fODEsolveAndUpdate;};"])
		mwrite (mfile, ["void write_names_to_file (FILE* output_c);",
						"void write_data_to_file_text (FILE* output_c);"])		
		mwrite (mfile, "};","\n", "\n")

		
		# define an operator, to apply all funcions 9operator for ode itrations run
		mwrite (mfile, ["void " + mfile.replace(' ','_') + "::operator() ()", "{"], "", "\n\n")
		#self.pattern_write (mfile, [fNames[:-1], fComment[:-1]], [ "\t", "();\t//", ""], self.data, range(len(fNames[:-1]))) # was used to miss ode step
		self.pattern_write (mfile, [fNames, fComment], [ "\t", "();\t//", ""], self.data, range(len(fNames)))
		#mwrite (mfile, "\t(*this.*step_it)();")
		mwrite (mfile, "}")
		
		# define function to print out cpp initial conditions
		#mwrite (mfile,["void " + mfile.replace(' ','_') + "print_ic()"]) 
		#mwrite (mfile, "{")
		#mwrite (mfile, )
		#mwrite (mfile, "}")

		
		# define all function data
		#rhs_eular = self.get_list ('rhs', self.data, "(", ") * dt")
		#for i in self.index['ode']['yes']:
		#	rhs_eular[i] = self.names[i] + " + " + rhs_eular[i]
		##for i in range(len(self.info['function'])):
		#	mcomment (mfile, self.info['function'][i]['comment'],"\n","\n")
		#	mwrite (mfile, "void " + mfile.replace(' ','_') + "::" + self.info['function'][i]['name'] + "()\n{")
		#	self.pattern_write (mfile, ["name", rhs_eular], [ " ", " = ", ";"], self.data, self.info['function'][i]['index'])
		#	mwrite (mfile, "}", "\n")
		index_print_nonode = self.new_index (self.index['ode']['no'], self.index['print']['yes'])
		lhs_names_list = self.new_names(self.new_names(self.names, index_print_nonode), self.index['ode']['yes'])
		for i in range(len(self.info['function'])):
			mcomment (mfile, self.info['function'][i]['comment'],"\n","\n")
			mwrite (mfile, "void " + mfile.replace(' ','_') + "::" + self.info['function'][i]['name'] + "()\n{")
			f_index = self.info['function'][i]['index']
			f_var_name = self.info['function'][i]['rhs']
			lhs_name = self.get_only_names (lhs_names_list, f_index)
			rhs_list = self.get_list ("rhs", self.data, "", "", f_index)
			for e in range(len(f_var_name)):
				if f_index[e] in self.index['ode']['yes']:
					rhs_list[e] = lhs_name[e] + " + (" + rhs_list[e] + ") * dt"
			self.pattern_write (mfile, [lhs_name,rhs_list], ["", " = ", ";"])
			mwrite (mfile, "}", "\n")

		
		# write string names to file
		mwrite (mfile, "void " + mfile + "::write_names_to_file  (FILE* output_c)\n{", "\n", "\n")
		mwrite (mfile, "fprintf(output_c,\"tim \"); ")
		self.pattern_write (mfile, "name", [ "fprintf(output_c,\"%s \", \"", "\"); "], self.data, self.index['print']['yes'],"")
		mwrite (mfile, "fprintf(output_c,\"\\n\");\n}")
		mwrite (mfile, "void " + mfile + "::write_data_to_file_text  (FILE* output_c)\n{", "\n", "\n")
		mwrite (mfile, "fprintf(output_c, \"%f \",tim); ","","")
		self.pattern_write (mfile, "name", [ "fprintf(output_c,\"%.18f \",", "); "], self.data,  self.index['print']['yes'],"")
		mwrite (mfile, "fprintf(output_c,\"\\n\");\n}")
		
		# close file
		mwrite (mfile, "\n#endif")
		self.close_file (mfile)
