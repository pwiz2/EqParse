#! /usr/bin/env python


import time

"""
 	Timer class for timing between multiple functions
	note: if timing only one line function, se built in python command: 
"""	

class Timer:

	def __init__(self):
		pass
		
	def start (self):
		self.t_start = time.clock()
		
	def finish (self, m_str = None):
		self.t_finish = time.clock()
		if self.t_start:
			if m_str: 
				print("Time taken for " + m_str + " " + str(self.t_finish - self.t_start) + " seconds")
			else: print("Time taken " + str(self.t_finish - self.t_start) + " seconds")
			self.t_start = None
		else: error("Timer finish() used before start()")

