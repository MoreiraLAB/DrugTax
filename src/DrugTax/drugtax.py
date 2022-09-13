#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Retrieve the kingdom and superclasses for an input smile.
Generate simple features
"""

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "DrugTax"

from .variables import *

def only_atoms_string(input_string):

	"""
	Retrieve only the atoms from an input smile
	"""
	import re
	characters_to_remove = ")(+-123456790-=#*.}{@"
	pattern = "[" + characters_to_remove + "]"
	return re.sub(pattern, "", input_string).replace("]","").replace("[","") 

def arbitrary_rings(input_string):

	"""
	Yield a molecule with all integers represented as &, that will allow its identification based on rings, regardless of the number of rings
	"""
	import re
	return re.sub("[1234567890]","&", input_string)

def character_count_features(input_smile, input_list =  CHARACTERS_LIST):

	"""
	Open a previously constructed text file with the available characters for SMILE format.
	Outputs dictionary with counts of listed characters 
	"""
	output_dictionary, usable_smile = {}, input_smile
	for entry in input_list:
		entry = entry.replace("\n","")
		output_dictionary["char_" + entry] = 0
		if (entry in input_smile):
			if len(entry) == 1:
				usable_smile = usable_smile.upper()
			output_dictionary["char_" + entry] = int(usable_smile.count(entry))
			usable_smile = usable_smile.replace(entry,"")
	return output_dictionary

def superclass_features_vector():

	"""
	Generate a feature vector that can be filled with superclass data
	"""
	return {"organic": 0, "inorganic": 0, "organoheterocyclic": 0, "benzenoid": 0, "organosulfur": 0, \
				"lipid": 0, "allene": 0, "phenylpropanoids_and_polyketides": 0, "carboxyl": 0, \
				"organic_acid": 0, "alkaloid": 0, "organic_salt": 0, "organic_nitrogen": 0, \
				"organic_oxygen": 0, "organophosphorus": 0, "aromatic_rings": 0, "organohalogens": 0, \
				"organometallics": 0, "nucleobases": 0, "sugars": 0, "homogenous_metal": 0, \
				"homogenous_non_metal": 0, "mixed_metal_non_metal": 0, "lignans_neolignans": 0, \
				"polymer_length": 0, "hydrocarbon_derivative": 0, "hydrocarbon": 0, "positive": 0, \
				"negative": 0, "carbene": 0, "organic_1_3_dipolar": 0, "inorganic_salt": 0, \
				"organopnictogen": 0, "acetylide": 0, "miscellaneous_inorganic": 0}

class DrugTax:

	def __init__(self, input_smile, input_type = "string"):

		"""
		Useful starting DrugTax attributes:
		- smile: the raw input SMILE
		- features: vector of features to be filled along superclass calculation
		- only_atoms_smile: removed non-atom characters from the input SMILE
		- carboxyl: whether or not there is a carboxyl group
		- superclasses: the vector onto which will be fed the superclass data
		- no_groups_smile: the SMILE without branched chemical groups
		- kingdom: organic or inorganic
		"""
		if input_type == "string":
			self.smile = input_smile
		elif input_type == "txt":
			self.smile = open(input_smile, "r").readlines()[0]

		self.features = {**superclass_features_vector(), **character_count_features(self.smile)} 
		self.smile = self.smile.replace("H","")
		self.only_atoms_smile = only_atoms_string(self.smile)
		self.carboxyl = False
		if ("C(=O)O" in self.smile):
			self.carboxyl = True
		self.identify_rings()
		self.identify_highlighted_molecules()
		self.superclasses = []
		import re
		self.no_groups_smile = "".join(re.split("\(|\)|\[|\]", self.smile)[::2])
		self.arbitrary_rings_smile = arbitrary_rings(self.smile)
		self.no_groups_arbitrary_rings_smile = arbitrary_rings(self.no_groups_smile)
		self.get_kingdom()
		if self.kingdom == "organic":
			self.get_organic_superclasses()
		elif self.kingdom == "inorganic":
			self.get_inorganic_superclasses()

	def get_kingdom(self):

		"""
		Organic compounds: 	compounds that contain at least one carbon atom, excluding fake organics
		Inorganic compounds: compounds that do not contain at least one carbon atom or fake organics
		Fake organics:
			- isocyanide/cyanide: [C-]#N non-hydrocarbyl derivatives
			- thiophosgene: C(=S)(Cl)Cl
			- carbon diselenide: C(=[Se])=[Se]
			- carbon monosulfide: [C-]#[S+]
			- carbon disulfide: C(=S)=S
			- carbon subsulfide: C(=C=S)=C=S
			- carbon monoxide: [C-]#[O+]  
			- carbon dioxide: C(=O)=O  
			- carbon suboxide: C(=C=O)=C=O
			- dicarbon monoxide: [C+]#C[O-]
		"""
		fake_organics = ["[C-]#N", "C(=S)(Cl)Cl", "C(=[Se])=[Se]", "[C-]#[S+]", \
							"C(=S)=S", "C(=C=S)=C=S", "[C-]#[O+]", "C(=O)=O", \
							"C(=C=O)=C=O", "C(=C=O)=C=O"]
		if (self.smile not in fake_organics) and ("C" in self.smile):
			self.kingdom = "organic"
			self.features["organic"] = 1
		else:
			self.kingdom = "inorganic"
			self.features["inorganic"] = 1

	def organo_hetero_cyclic(self):

		"""
		1 - Organoheterocyclic compounds: compounds containing a ring with least one carbon atom and one non-carbon atom.
		5 - Benzenoids: Aromatic compounds containing one or more benzene rings - only carbon atoms in ring. #For aromatic rings, the carbons are written in lower letter
		"""
		heterocyclic, benzenoid = False, False
		if len(self.rings) > 0:
			for iter_rings in self.rings.keys():
				only_atoms_ring_string = only_atoms_string(self.rings[iter_rings])
				if ((len(only_atoms_ring_string.replace("C","")) > 0) or len(only_atoms_ring_string.replace("c","")) > 0) and \
					((only_atoms_ring_string.count("C") >= 1) or (only_atoms_ring_string.count("c") >= 1)) and (heterocyclic == False):
					self.superclasses.append("organoheterocyclic")
					self.features["organoheterocyclic"] = 1
					heterocyclic = True
				elif (len(only_atoms_ring_string.replace("c","")) == 0) and \
					(benzenoid == False):
					self.superclasses.append("benzenoids")
					self.features["benzenoid"] = 1
					benzenoid = True
				

	def sulfur_bond(self):

		"""
		2 - Organosulfur compounds: compounds containing a carbon-sulfur bond.
		""" 
		registered_superclass = False
		for current_sulf_bond in ["SC", "CS"]:
			if current_sulf_bond in self.only_atoms_smile:
				if registered_superclass == False: 
					self.superclasses.append("organosulfur")
					registered_superclass = True
				self.features["organosulfur"] += self.only_atoms_smile.count(current_sulf_bond)

	def lipids(self):

		"""
		3 - Lipids and lipid-like molecules: Fatty acids and their derivatives, and substances related biosynthetically or functionally to these compounds.
		https://www.lipotype.com/lipidomics-services/fatty-acyl-analysis/
		Excluded short chain fatty-acids
		"""

		saturation_blind_smile = self.smile.replace("=","").replace("#","")
		#Fatty Acids
		if (self.carboxyl == True) and ("C"*4 in saturation_blind_smile):
			self.superclasses.append("lipids")
			self.features["lipid"] = 1

	def allenes(self):

		"""
		4 - Allenes: Compounds in which one carbon atom has double bonds with each of its two adjacent carbon centers.
		The definition includes both the hydrocarbon molecules and their derivatives obtained by substitution.
		"""
		split_bounding = self.smile.split("=C=")
		if len(split_bounding) == 1:
			return
		if (split_bounding[0][-1] in ["C",")"]) and (split_bounding[1][0] in ["C","("]):
			self.superclasses.append("allenes")
			self.features["allene"] = 1

	def phenylpropanoids_and_polyketides(self):

		"""
		6 - Phenylpropanoids and polyketides: Organic compounds that are synthesized either from the amino acid phenylalanine (phenylpropanoids) 
		or the decarboxylative condensation of malonyl-CoA (polyketides). Phenylpropanoids are aromatic compounds based on the phenylpropane skeleton. 
		Polyketides usually consist of alternating carbonyl and methylene groups (beta-polyketones), biogenetically derived from repeated condensation of
		 acetyl coenzyme A (via malonyl coenzyme A), and usually the compounds derived from them by further condensations.

		phenilpropane: CCCC1=CC=CC=C1
		carbonyl: C=O
		methylene: [CH2]
		"""
		if ("benzenoids" in self.superclasses) and ("C"*3 in self.no_groups_arbitrary_rings_smile):
			self.superclasses.append("phenylpropanoids_and_polyketides")
			self.features["phenylpropanoids_and_polyketides"] = 1
			return

		if ("C=O[CH2]" in self.smile) or ("[CH2]C=O" in self.smile) or \
			("CC=O" in self.smile) or ("C=OC" in self.smile):
			self.superclasses.append("phenylpropanoids_and_polyketides")
			self.features["phenylpropanoids_and_polyketides"] = 1
			return

	def organic_acids(self):

		"""
		7 - Organic acids: Compounds an organic acid or a derivative thereof.
		https://www.sciencedirect.com/topics/neuroscience/organic-acids
		Organic acids without carboxyl group:	https://pubs.acs.org/doi/pdf/10.1021/ed077p910

		"""
		if (self.carboxyl == True):
			#Carboxyl acids
			self.superclasses.append("organic_acids")
			self.features["carboxyl"] = 1
			self.features["organic_acid"] = 1


	def alkaloids(self):

		"""
		8 - Alkaloids and derivatives: Naturally occurring chemical compounds that contain mostly basic nitrogen atoms. 
		This group also includes some related compounds with neutral and even weakly acidic properties. 
		Also some synthetic compounds of similar structure are attributed to alkaloids. In addition to carbon, hydrogen and nitrogen, 
		alkaloids may also contain oxygen, sulfur and more rarely other elements such as chlorine, bromine, and phosphorus.
		"""
		if ("N" in self.only_atoms_smile) and \
				((int(self.smile.count("-")) - int(self.smile.count("+"))) > 0):
			self.superclasses.append("alkaloids")
			self.features["alkaloid"] = 1

	def organic_salts(self):

		"""
		9 - Organic salts: Organic compounds consisting of an assembly of cations and anions.
		"""
		split_by_sep = self.smile.split("{")
		if len(split_by_sep) > 1:
			split_block = split_by_sep[1]
		else:
			split_by_sep = self.smile.split("[")
			if len(split_by_sep) > 1:
				split_block = split_by_sep[1]
			else:
				split_block = ""
		if ("+" in split_block) or ("-" in split_block):
			self.superclasses.append("organic_salts")
			self.features["organic_salt"] = 1

	def organic_count(self, atom_type = "N"):

		"""
		12 - Organic nitrogen compounds: Organic compounds containing a nitrogen atom.
		14 - Organic oxygen compounds: Organic compounds that contain one or more oxygen atoms.
		"""
		if atom_type == "N":
			superclass_name = "organic_nitrogen"
		if atom_type == "O":
			superclass_name = "organic_oxygen"
		if atom_type == "P":
			superclass_name = "organophosphorus"

		if atom_type in self.only_atoms_smile:
			self.superclasses.append(superclass_name)
			self.features[superclass_name] = 1

	def identify_rings(self):

		"""
		Identify, if possible, aromatic rings, retrieve a dictionary with the list of rings
		"""
		count = 1
		self.rings = {}
		possible_ring = True
		while possible_ring == True:
			current_position_ring = str(count)
			try:
				current_ring = self.smile.split(current_position_ring)[0][-1] + current_position_ring + \
						self.smile.split(current_position_ring)[1] + current_position_ring
				self.rings[count] = current_ring
				count += 1
			except:
				possible_ring = False

		if bool(self.rings) != False:
			self.features["aromatic_rings"] = count

	def identify_highlighted_molecules(self):

		"""
		Identify, if possible, highlighted molecules, retrieve a dictionary with the list of these molecules
		"""

		count = 1
		self.molecules = {}
		possible_molecule = True
		while possible_molecule == True:
			current_position_molecule = str(count)
			try:
				current_molecule = self.smile.split("[")[count]
				self.molecules[count] = current_molecule.split("]")[0]
				count += 1
			except:
				possible_molecule = False

	def bounds_to_single_atoms(self, input_mode = ""):

		"""
		10 - Organohalogen compounds: Organic compounds containing a bond between a carbon atom and a halogen atom (At, F, Cl, Br, I).
		11 - Organometallic compounds: Organic compounds containing a bond between a carbon atom and metal atom.
		"""
		if input_mode == "organohalogens":
			usable_atoms = HALOGENS

		if input_mode == "organometallics":
			usable_atoms = METALS

		registered_superclass = False
		for current_atom in usable_atoms:
			bond = "C" + current_atom
			bond_reverse = current_atom + "C"
			aromatic = "c" + current_atom
			aromatic_reverse = current_atom + "c"
			if (bond in self.only_atoms_smile) or (bond_reverse in self.only_atoms_smile) or \
				(aromatic in self.only_atoms_smile) or (aromatic_reverse in self.only_atoms_smile):
				if registered_superclass == False:
					self.superclasses.append(input_mode)
					registered_superclass = True
				self.features[input_mode] += 1

	def nucleotides_nucleosides(self):

		"""
		13 - Nucleosides and nucleotides: Compounds containing a nucleobase linked to a ribose or deoxyribose sugar 
		via a beta-glycosidic linkage. The ribose or deoxyribose moiety bears at least a phosphate group in case of nucleotides.
		cytosine: C1=C(NC(=O)N=C1)N
		adenine: C1=NC2=NC=NC(=C2N1)N
		guanine: C1=NC2=C(N1)C(=O)NC(=N2)N
		thymine: CC1=CNC(=O)NC1=O
		uracil: C1=CNC(=O)NC1=O
		ribose: C([C@H]([C@H]([C@H](C=O)O)O)O)O
		deoxyribose: C(C=O)[C@@H]([C@@H](CO)O)O
		"""
		import re
		nucleobases = ["C&=C(NC(=O)N=C&)N", "C&=NC&=NC=NC(=C&N&)N", "C&=NC&=C(N&)C(=O)NC(=N&)N", \
						"CC&=CNC(=O)NC&=O", "C&=CNC(=O)NC&=O"]
		sugars = ["C([C@H]([C@H]([C@H](C=O)O)O)O)O","C(C=O)[C@@H]([C@@H](CO)O)O"]
		first_condition, second_condition = False, False
		for current_nucleobase in nucleobases:
			if current_nucleobase in self.arbitrary_rings_smile:
				self.features["nucleobases"] += 1
				first_condition = True

		for current_sugar in sugars:
			if current_sugar in self.arbitrary_rings_smile:
				self.features["sugars"] += 1
				second_condition = True

		if first_condition == True and second_condition == True:
			self.superclasses.append("nucleosides_nucleotides")

	def homogenous_metal(self):

		"""
		1 - Homogenous metal compounds: Inorganic compounds that contain only metal atoms.
		"""
		new_string = self.only_atoms_smile
		for current_atom in METALS:
			new_string = new_string.replace(current_atom, "")
		if new_string == "":
			self.superclasses.append("homogenous_metal")
			self.features["homogenous_metal"] = 1

	def homogenous_non_metal(self):

		"""
		2 - Homogenous non-metal compounds: Inorganic compounds that contain only non-metal atoms.
		"""
		new_string = self.only_atoms_smile
		for current_atom in METALS:
			new_string = new_string.replace(current_atom, "")
		if new_string == self.only_atoms_smile:
			self.superclasses.append("homogenous_non_metal")
			self.features["homogenous_non_metal"] = 1

	def mixed_metal_non_metal(self):

		"""
		3 - Mixed metal/non-metal compounds: Inorganic compounds that contain both metal and non metal atoms.
		"""
		new_string = self.only_atoms_smile
		for current_atom in METALS:
			new_string = new_string.replace(current_atom, "")
		if (new_string != self.only_atoms_smile) and (new_string != ""):
			self.superclasses.append("mixed_metal_non_metal")
			self.features["mixed_metal_non_metal"] = 1

	def lignans_neolignans(self):

		"""
		16 - Lignans, neolignans and related compounds: Plant products of low molecular weight formed primarily from oxidative coupling of two \
		p-propylphenol moieties. They can also be described as micromolecules with two phenylpropanoid units coupled together. \
		They can be attached in various manners, like C5-C5', C8-C8'. Most known natural lignans are oxidized at C9 and C9´ and, \
		based upon the way in which oxygen is incorporated into the skeleton and on the cyclization patterns, a wide range of lignans \
		of very different structural types can be formed.
		p-propyphenol: CCCC1=CC=C(C=C1)O
		phenylpropane: CCCC1=CC=CC=C1
		"""
		if len(self.arbitrary_rings_smile.split("C&=CC=C(C=C&)O")) >= 2:
			self.superclasses.append("lignans_neolignans")
			self.features["lignans_neolignans"] = 1
			return
		elif len(self.arbitrary_rings_smile.split("CCCC&=CC=CC=C&")) >= 2:
			self.superclasses.append("lignans_neolignans")
			self.features["lignans_neolignans"] = 1
			return

	def organic_polymers(self, min_pattern_size = 5):

		"""
		17 - Organic Polymers: Organic compounds, generally large molecules or macromolecules, which are composed of many repeating units.
		"""
		import re
		regex_string = "\.{" + str(min_pattern_size) + "}"
		pattern = re.compile(regex_string)
		located_patterns = pattern.findall(self.smile)
		for current_pattern in located_patterns:
			if int(self.smile.count(current_pattern)) > 2:
				self.superclasses.append("organic_polymers")
				self.features["polymer_length"] += int(self.smile.count(current_pattern))
				return

	def hydrocarbon_derivatives(self):

		"""
		18 - Hydrocarbon derivatives: Derivatives of hydrocarbons obtained by substituting one or more carbon atoms by an heteroatom. They contain at least one carbon atom and heteroatom.
		19 - Hydrocarbons: Organic compounds made up only of carbon and hydrogen atoms.
		"""
		if self.only_atoms_smile.replace("C","") != "":
			self.superclasses.append("hydrocarbon_derivatives")
			self.features["hydrocarbon_derivative"] = 1
		else:
			self.superclasses.append("hydrocarbon")
			self.features["hydrocarbon"] = 1

	def charged_categories(self, target_charges = ""):

		"""
		20 - Organic Anions: Organic compounds that have a negative electric charge.
		21 - Organic Cations: Organic compounds with a positive electric charge.
		22 - Organic Zwitterions: Organic neutral compounds having formal unit electrical charges of opposite sign.
		"""
		self.charges_log_dictionary = {"+": int(self.smile.count("+")), "-": int(self.smile.count("-"))}
		self.features["positive"] = self.charges_log_dictionary["+"]
		self.features["negative"] = self.charges_log_dictionary["-"]
		if (self.charges_log_dictionary["+"] == 0) and (self.charges_log_dictionary["-"] == 0):
			return

		if (self.charges_log_dictionary["+"] - self.charges_log_dictionary["-"]) == 0:
			self.superclasses.append("organic_zwitterions")
			return

		if (self.charges_log_dictionary["+"] - self.charges_log_dictionary["-"]) > 0:
			self.superclasses.append("organic_cations")
			return

		if (self.charges_log_dictionary["+"] - self.charges_log_dictionary["-"]) < 0:
			self.superclasses.append("organic_anions")
			return

	def carbenes(self):

		"""
		23 - Carbenes: The electrically neutral species H2C and its derivatives, in which the carbon is covalently 
		bonded to two univalent groups of any kind or a divalent group and bears two non-bonding electrons, which may be 
		spin-paired (singlet state) or spin-non-paired (triplet state). Subclasses of carbenes include acyl carbenes, imidoyl carbenes, and vinyl carbenes.
		Note: In a SMILE, "." denotes unpaired electrons
		"""
		for current_carbene_config in ["[C-]", "[C--]","[C---]","[C-2]","[C-3]"]:
			if current_carbene_config in self.smile:
				self.superclasses.append("carbenes")
				self.features["carbene"] = 1
				return
				
		if "." in self.smile:
			split_segments = self.smile.split(".")
			for current_segment in split_segments[0:-1]:
				clean_segment = only_atoms_string(current_segment)
				if clean_segment[-1] == "C":
					self.superclasses.append("carbenes")
					self.features["carbene"] = 1
					return

	def organic_1_3_dipolar(self):

		"""
		24 - Organic 1,3-dipolar compounds: Electrically neutral organic molecules carrying a positive and a negative charge in one of their major canonical descriptions.
		 In most dipolar compounds the charges are delocalized; however the term is also applied to species where this is not the case. 
		 The term 1,3-dipolar compounds is used for those in which a significant canonical resonance form can be represented by a separation of charge over three atoms 
		 (in connection with 1,3-dipolar cycloadditions).
		 The dipole has at least one resonance structure with positive and negative charges having a 1,3 relationship which can generally be denoted as a+ - b - c-
		  where "a" may be a carbon, oxygen or nitrogen, "b" may be nitrogen or oxygen and "c" may be a carbon, oxygen or nitrogen
		"""
		if (self.charges_log_dictionary["+"] != 0) and (self.charges_log_dictionary["-"] != 0) and \
			((self.charges_log_dictionary["+"] - self.charges_log_dictionary["-"]) == 0):
			for a in ["C","O","N"]:
				for b in ["N","O"]:
					for c in ["C","O","N"]:
						current_combos = ["[" + a + "+]" + b + "[" + c + "-]", "["+ a + "-]" + b + "[" + c + "+]"]
						if current_combos[0] in self.smile:
							self.superclasses.append("organic_1_3_dipolar")
							self.features["organic_1_3_dipolar"] = 1
							return
						if current_combos[1] in self.smile:
							self.superclasses.append("organic_1_3_dipolar")
							self.features["organic_1_3_dipolar"] = 1
							return

	def inorganic_salts(self):
	
		"""
		5 - Inorganic salts: Inorganic compounds consisting of an assembly of cations and anions.
		"""
		split_smile_bracket = self.smile.split("]")
		if len(split_smile_bracket) == 1:
			return
		charges_log_dictionary = {"+": 0, "-": 0}
		for current_block in split_smile_bracket[0:-1]:
			current_charge = current_block[-1]
			if current_charge in charges_log_dictionary:
				charges_log_dictionary[current_charge] += 1
		if (charges_log_dictionary["+"] > 0) and (charges_log_dictionary["-"] > 0):
			self.superclasses.append("inorganic_salts")
			self.features["inorganic_salt"] = 1
			self.features["positive"] = charges_log_dictionary["+"]
			self.features["negative"] = charges_log_dictionary["-"]

	def organopnictogen(self):

		"""
		25 - Organopnictogen compounds: Compounds containing a bond between carbon and a pnictogen atom. Pnictogens are p-block element atoms that are in the group 15 of the periodic table.
		"""
		upper_only_atoms_smile = self.only_atoms_smile.upper()
		for current_pnictogen in GROUP_15:
			if (("C" + current_pnictogen) in upper_only_atoms_smile) or ((current_pnictogen + "C") in upper_only_atoms_smile) or \
					(("C=" + current_pnictogen) in upper_only_atoms_smile) or ((current_pnictogen + "=C") in upper_only_atoms_smile) or \
					(("C#" + current_pnictogen) in upper_only_atoms_smile) or ((current_pnictogen + "#C") in upper_only_atoms_smile) or \
					(("C+" + current_pnictogen) in upper_only_atoms_smile) or ((current_pnictogen + "*C") in upper_only_atoms_smile):
				self.superclasses.append("organopnictogens")
				self.features["organopnictogen"] = 1
				return

	def acetylide(self):

		"""
		26 - Acetylides: chemical compounds with the chemical formulas MC≡CH and MC≡CM, where M is a metal. 
		"""
		for current_metal in METALS:
			for current_metal_2 in METALS:
				formula_1_A = current_metal + "C#C"
				formula_1_B = "C#C" + current_metal
				formula_2 = current_metal + "C#C" + current_metal_2
				if (formula_1_A in self.smile) or (formula_1_B in self.smile) or (formula_2 in self.smile):
					self.superclasses.append("acetylides")
					self.features["acetylide"] = 1
					return

	def get_organic_superclasses(self):

		"""
		26 possible superclasses
		1 - Organoheterocyclic compounds: compounds containing a ring with least one carbon atom and one non-carbon atom.
		2 - Organosulfur compounds: compounds containing a carbon-sulfur bond.
		3 - Lipids and lipid-like molecules: Fatty acids and their derivatives, and substances related biosynthetically or functionally to these compounds.
		4 - Allenes: Compounds in which one carbon atom has double bonds with each of its two adjacent carbon centers. The definition includes both the hydrocarbon molecules and their derivatives obtained by substitution.
		5 - Benzenoids: Aromatic compounds containing one or more benzene rings.
		6 - Phenylpropanoids and polyketides: Organic compounds that are synthesized either from the amino acid phenylalanine (phenylpropanoids) 
		or the decarboxylative condensation of malonyl-CoA (polyketides). Phenylpropanoids are aromatic compounds based on the phenylpropane skeleton. 
		Polyketides usually consist of alternating carbonyl and methylene groups (beta-polyketones), biogenetically derived from repeated condensation of
		 acetyl coenzyme A (via malonyl coenzyme A), and usually the compounds derived from them by further condensations.
		7 - Organic acids and derivatives: Compounds an organic acid or a derivative thereof.
		8 - Alkaloids and derivatives: Naturally occurring chemical compounds that contain mostly basic nitrogen atoms. 
		This group also includes some related compounds with neutral and even weakly acidic properties. 
		Also some synthetic compounds of similar structure are attributed to alkaloids. In addition to carbon, hydrogen and nitrogen, 
		alkaloids may also contain oxygen, sulfur and more rarely other elements such as chlorine, bromine, and phosphorus.
		9 - Organic salts: Organic compounds consisting of an assembly of cations and anions.
		10 - Organohalogen compounds: Organic compounds containing a bond between a carbon atom and a halogen atom (At, F, Cl, Br, I).
		11 - Organometallic compounds: Organic compounds containing a bond between a carbon atom and metal atom.
		12 - Organic nitrogen compounds: Organic compounds containing a nitrogen atom.
		13 - Nucleosides, nucleotides, and analogues: Compounds containing a nucleobase linked to a ribose or deoxyribose sugar 
		via a beta-glycosidic linkage. The ribose or deoxyribose moiety bears at least a phosphate group in case of nucleotides.
		14 - Organic oxygen compounds: Organic compounds that contain one or more oxygen atoms.
		15 - Organophosphorus compounds: Organic compounds containing the phosphorus atom.
		16 - Lignans, neolignans and related compounds: Plant products of low molecular weight formed primarily from oxidative coupling of two \
		p-propylphenol moieties. They can also be described as micromolecules with two phenylpropanoid units coupled together. \
		They can be attached in various manners, like C5-C5', C8-C8'. Most known natural lignans are oxidized at C9 and C9´ and, \
		based upon the way in which oxygen is incorporated into the skeleton and on the cyclization patterns, a wide range of lignans \
		of very different structural types can be formed.
		17 - Organic Polymers: 	Organic compounds, generally large molecules or macromolecules, which are composed of many repeating units.
		18 - Hydrocarbon derivatives: Derivatives of hydrocarbons obtained by substituting one or more carbon atoms by an heteroatom. They contain at least one carbon atom and heteroatom.
		19 - Hydrocarbons: Organic compounds made up only of carbon and hydrogen atoms.
		20 - Organic Anions: Organic compounds that have a negative electric charge.
		21 - Organic Cations: Organic compounds with a positive electric charge.
		22 - Organic Zwitterions: Organic neutral compounds having formal unit electrical charges of opposite sign.
		23 - Carbenes: The electrically neutral species H2C and its derivatives, in which the carbon is covalently 
		bonded to two univalent groups of any kind or a divalent group and bears two non-bonding electrons, which may be 
		spin-paired (singlet state) or spin-non-paired (triplet state). Subclasses of carbenes include acyl carbenes, imidoyl carbenes, and vinyl carbenes.
		24 - Organic 1,3-dipolar compounds: Electrically neutral organic molecules carrying a positive and a negative charge in one of their major canonical descriptions.
		 In most dipolar compounds the charges are delocalized; however the term is also applied to species where this is not the case. 
		 The term 1,3-dipolar compounds is used for those in which a significant canonical resonance form can be represented by a separation of charge over three atoms 
		 (in connection with 1,3-dipolar cycloadditions).
		25 - Organopnictogen compounds: Compounds containing a bond between carbon a pnictogen atom. Pnictogens are p-block element atoms that are in the group 15 of the periodic table.
		26 - Acetylides: chemical compounds with the chemical formulas MC≡CH and MC≡CM, where M is a metal. 
		"""

		self.sulfur_bond()
		self.organo_hetero_cyclic()
		self.lipids()
		self.allenes()
		self.phenylpropanoids_and_polyketides()
		self.organic_acids()
		self.alkaloids()
		self.organic_salts()

		self.bounds_to_single_atoms(input_mode = "organohalogens")
		self.bounds_to_single_atoms(input_mode = "organometallics")

		self.organic_count(atom_type = "N")
		self.nucleotides_nucleosides()
		self.organic_count(atom_type = "O")
		self.organic_count(atom_type = "P")
		self.lignans_neolignans()
		self.organic_polymers()
		self.hydrocarbon_derivatives()
		self.charged_categories() #Three possible superclasses in this bad boy
		self.carbenes()
		self.organic_1_3_dipolar() #Requires charged_categories to define charge dictionary
		self.organopnictogen()
		self.acetylide()

	def get_inorganic_superclasses(self):

		"""
		5 possible superclasses
		1 - Homogenous metal compounds: Inorganic compounds that contain only metal atoms.
		2 - Homogenous non-metal compounds: Inorganic compounds that contain only non-metal atoms.
		3 - Mixed metal/non-metal compounds: Inorganic compounds that contains both metal and non metal atoms.
		4 - Miscellaneous inorganic compounds: Inorganic compounds of miscellaneous structure, which can be a mixture of chemical elements of different types.
		5 - Inorganic salts: Inorganic compounds consisting of an assembly of cations and anions.
		"""

		self.homogenous_metal()
		self.homogenous_non_metal()
		self.mixed_metal_non_metal()
		self.inorganic_salts()

		if len(self.superclasses) == 0:
			self.superclasses.append("miscellaneous_inorganic")
			self.features["miscellaneous_inorganic"] = 1

	def __repr__(self):

		return  "SMILE: " + self.smile + "\nKingdom: " + self.kingdom +"\nSuperclasses: " + "-".join(sorted(self.superclasses))