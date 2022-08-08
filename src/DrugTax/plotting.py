#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Given that different ligands can have multiple superclasses, stack these and plot them in a useful manner using upsetplots
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "DrugTax"

import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from upsetplot import plot as upset_plot
from upsetplot import from_memberships

def merge_categories(input_table, target_threshold = 10, \
						merging_mode = "min", report_changes = False, \
						report_name = ""):
	
	"""
	In order to avoid too large a list of unique categories, use this function to merge them if adequate
	- target_threshold = 10: change this to limit the number of superclass combinations 
	"""
	usable_superclasses, below_thresh_superclasses = {}, {}

	for index, current_entry in input_table.iterrows():
		current_name = current_entry["Superclasses"].split("-")
		current_count = int(current_entry["Count"])
		if current_count <= target_threshold:
			below_thresh_superclasses[current_entry["Superclasses"]] = current_count
		elif current_count > target_threshold:
			usable_superclasses[current_entry["Superclasses"]] = current_count

	import itertools
	convertion_table = []
	for current_key in below_thresh_superclasses.keys():
		split_key = current_key.split("-")
		max_length = len(split_key)
		possible_combinations = [list(x) for x in itertools.combinations(split_key, max_length - 1)]
		first = True

		for current_combination in possible_combinations:
			remerged_string = "-".join(current_combination)
			if remerged_string not in usable_superclasses.keys():
				continue

			current_fill = usable_superclasses[remerged_string] 
			if first == True:
				max_existing, min_existing = remerged_string, remerged_string 
				min_fill, max_fill = current_fill, current_fill
				first = False
			
			if current_fill < min_fill:
				min_fill = current_fill
				min_existing = remerged_string

			elif current_fill > max_fill:
				max_fill = current_fill
				max_existing = remerged_string
		
		if first == True:
			#If it was not possible to find any combination above the target one, keep the original
			usable_superclasses[current_key] = below_thresh_superclasses[current_key]
			continue

		if merging_mode == "max":
			usable_superclasses[max_existing] += below_thresh_superclasses[current_key]
			convertion_table.append([current_key, max_existing])

		if merging_mode == "min":
			usable_superclasses[min_existing] += below_thresh_superclasses[current_key]
			convertion_table.append([current_key, min_existing])

	if report_changes == True:
		convertion_dataframe = pd.DataFrame(convertion_table, columns = ["original", "merged"])
		convertion_dataframe.to_csv(report_name + "_report_changes.csv", index = False)

	cat_list, cat_count = [], []
	for proper_key in usable_superclasses.keys():
		split_proper_key = proper_key.split("-")
		cat_list.append(split_proper_key)
		cat_count.append(usable_superclasses[proper_key])
	return cat_list, cat_count

def plot_categories(input_file, threshold = 1, output_name = ""):

	"""
	Open the input file and analyse category overlap. Produces an upset plot with the category
	information as well as superclass occurrence analysis
	- threshold = 1: change the threshold to merge under populated categories.
	Categories will be aggregated, when possible, to a less specific category (less superclasses)
	- output_name = "": change the outputname of your png file
	"""
	opened_table = pd.read_csv(input_file, sep = ",", header = 0)
	opened_table.columns = ["Superclasses","Count"]
	organic = ["organoheterocyclic","organosulfur","lipids","allenes","benzenoids","phenylpropanoids_and_polyketides",\
				"organic_acids", "alkaloids", "organic_salts", "organic_nitrogen", "organic_oxygen", "organophosphorus", \
				"organohalogens", "organometallics", "nucleosides_nucleotides_analogues", "lignans_neolignans_and_related", \
				"organic_polymers", "hydrocarbon_derivatives", "hydrocarbon", "organic_zwitterions", "organic_cations", \
				"organic_anions", "carbenes", "organic_1_3_dipolar", "organopnictogens", "acetylides"]
	inorganic = ["homogenous_metal", "homogenous_non_metal", "mixed_metal_non_metal", "inorganic_salts", "miscellaneous_inorganic"]
	total = opened_table["Count"].sum()
	categories_list, categories_counts = merge_categories(opened_table, target_threshold = threshold, \
							report_changes = True, report_name = output_name, merging_mode = "min")
	data_object = from_memberships(categories_list, data = categories_counts)
	plt.figure(figsize = (50, 50))
	plt.legend(prop={'size': 12})
	upset_plot(data_object)
	plt.savefig(output_name + ".png")