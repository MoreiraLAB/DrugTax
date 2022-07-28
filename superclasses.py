#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Call DrugTax on unique SMILEs
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "DrugTax"

import pandas as pd
import os
import sys
from drugtax import DrugTax

def download_unique_smiles(input_drugs_list, verbose = False):

	"""
	Download unique smiles from a list of input ligand names, making use of PubChempy
	"""
	import os.path
	from os import path
	from pubchempy import get_compounds
	smiles_dictionary = {}
	for current_drug in input_drugs_list:
		try:
			current_smile = get_compounds(current_drug, 'name')[0].isomeric_smiles
			smiles_dictionary[current_drug] = current_smile
			if verbose == True:
				print("Successfully downloaded:", current_drug)
		except:
			if verbose == True:
				print("Failed to download:", current_drug)
			else:
				continue
	return smiles_dictionary

def retrieve_taxonomic_class(input_smiles, input_mode = "file", target_column = "", \
								output_name = "", write_values = False, \
								input_sep = ""):
	"""
	Retrieve the taxonomic classes for an input list of SMILES. The input_mode parameter should be adapted if needed:
	- input_smiles: this argument is mandatory, should be used in conjunction with the other arguments depending on the input type
	can correspond an input file name (input_mode = "file"), a list of SMILEs (input_mode = "smiles_list") 
	or a list of drug names (input_mode = "drugs_list").
	- input_mode = "file": will open the file as if it was a table, 
	should be coupled with the target_column and input_sep arguments
		- target_column = "": should correspond to the name of the column listing the input smiles
		- input_sep = "": should be changed to the column separator used, when considering an input file
	- input_mode = "smiles_list": will compute a list of input SMILEs
	- input_mode = "drugs_list": will attempt to retrieve the SMILEs for input drug names
	- write_values = True: if set to True this outputs two tables, one with the DrugTax classification of each drug, 
	another with the superclass occurrence for each drug (necessary to generate UpSetPlot)
	- output_name = "": should be specified to yield proper outputs with write_values = True
	"""

	if input_mode == "file":
		input_table = pd.read_csv(input_smiles, sep = input_sep, header = 0)
		smile_column = list(input_table[target_column])

	if input_mode == "drugs_list":
		molecules_dictionary = download_unique_smiles(input_smiles)
		input_table = pd.DataFrame(molecules_dictionary.items(), columns = ["NAME","SMILE"])
		smile_column = list(input_table["SMILE"])

	if input_mode == "smiles_list":
		smile_column = input_smiles
		input_table = pd.DataFrame(smile_column, columns = ["SMILE"])

	tax_classes = []
	for index, current_smile in enumerate(smile_column):
		if ((index + 1) % 100) == 0:
			print("Currently:", index + 1, "/", len(smile_column))
		current_tax_class = str(DrugTax(current_smile))
		tax_classes.append(current_tax_class)
	
	input_table["Taxonomy"] = tax_classes
	assessment_table = input_table["Taxonomy"].value_counts()
	if write_values == True:
		input_table.to_csv(output_name + ".csv", sep = ",", index = None)
		assessment_table.to_csv(output_name.split(".")[0] + "_assess.csv", sep = ",")
	return input_table, assessment_table
