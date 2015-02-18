#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A constant-space parser for the GeneOntology OBO v1.2 format

Version 1.0
"""
from __future__ import with_statement
from collections import defaultdict


class OBOParser(object):

	def __init__(self,path_to_obo_file="../docs/go-basic.obo"): #Dangerous hard-coding
		self.path_to_obo_file = path_to_obo_file
		self.ontology_dictionary = {}
		self.load_obo()

	def __iter__(self):
		for key,value in self.ontology_dictionary.iteritems():
			yield key,value

	def load_obo(self):
		print 'bob'
		with open(self.path_to_obo_file, "r") as infile:
			print 'here'
			currentGOTerm = None
			for line in infile:
				print line
				line = line.strip()
				if not line: continue #Skip empty
				if line == "[Term]":
					if currentGOTerm: self.ontology_dictionary[currentGOTerm] =  self.processGOTerm(currentGOTerm)
					currentGOTerm = defaultdict(list)
				elif line == "[Typedef]":
					#Skip [Typedef sections]
					currentGOTerm = None
				else: #Not [Term]
					#Only process if we're inside a [Term] environment
					if currentGOTerm is None: continue
					key, sep, val = line.partition(":")
					currentGOTerm[key].append(val.strip())
			#Add last term
			if currentGOTerm is not None:
				self.ontology_dictionary[currentGOTerm] = self.processGOTerm(currentGOTerm)

	def processGOTerm(self, goTerm):
		"""
		In an object representing a GO term, replace single-element lists with
		their only member.
		Returns the modified object as a dictionary.
		"""
		ret = dict(goTerm) #Input is a defaultdict, might express unexpected behaviour
		for key, value in ret.iteritems():
			if len(value) == 1:
				ret[key] = value[0]
		return ret
