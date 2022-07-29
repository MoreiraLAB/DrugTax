# DrugTax
Categorize small ligands according to chemical properties. Derive simple and explainable features. Only requires SMILEs as inputs. For a more detailed description of the reasoning behind DrugTax address the scientific paper *DrugTax: package for drug taxonomy identification and explainable feature extraction*, Preto, AJ *et al.*, 2022.

# Installation

To install **DrugTax**, first make sure you have **python 3.9.x** installed. Then, run:

`pip install drugtax`

# Usage
Firstly, import DrugTax

`import drugtax`

The basic usage of DrugTax stems from the DrugTax class, which takes as input a single SMILE.

`molecule = drugtax.DrugTax("OC(=O)C1=C(C(O)=O)C(C(O)=O)=C(C(O)=O)C(C(O)=O)=C1C(O)=O")`

The `molecule` object now has a series of useful properties such as:

- `molecule.smile`: allows the user to check the SMILE at any time
- `molecule.superclasses`: displays all the superclasses to which the input SMILE belongs
- `molecule.features`: retrieves simple and explainable features that can be used on prediction tasks or dataset characterization
- `molecule.kingdom`: informs on whether the molecule is *organic* or *inorganic*

# Bulk analysis

For superclass computation, instead of directly invoking the **DrugTax** class, it is possible to use `retrieve_taxonomic_class` on several different inputs. The example below shows an example using only a SMILEs list, however, it is also possible to feed a drugs list - in which case DrugTax leverages pubchempy to retrieve the isomeric SMILEs - or a file. This function outputs a table with the SMILEs and their respective taxonomy, as well as a summary table detailing how many of which superclass combinations are present on the dataset.

`smiles_table, summary_table = drugtax.retrieve_taxonomic_class(["CCNO","CCC"], input_mode = "smiles_list", output_name = "testing", write_values = True)`

The `retrieve_taxonomic_class` function has different arguments that can be used to pick input and output information for bulk analysis submission. These are:
- `input_data`: this is the only mandatory argument, corresponding to either a smiles list, a drugs names list or a file.
- `input_mode`: depending on the `input_data` this argument needs to be changed. The default mode is `file`, which requires an input `.csv` file. This argument needs to be coupled with `target_column`, specifying the name of the column containing the input SMILEs. To input a list of SMILE, `input_mode` needs to be changed to `smiles_list`. To input a list of drug names from which isomeric SMILEs are to be retrieved with the aid of `pubchempy`, the user needs to change `input_mode` to `drugs_list`.
- `target_column`: this argument needs to be specified when using the `file` input mode.
- `output_name`: if the user wishes to save the file, he should specify this argument.
- `write_values`: with default `False`, the user needs to change this argument to `True` if he wishes to save the output tables.
- `input_sep`: if the input mode is `file`, the user can change this argument depending on the column separator on the input file. 

# Plotting

In order to visualize the data retrieved from bulk analysis, **DrugTax** leverages `UpSetPlot`, a package designed to allow the visualization of a large number of intersecting sets. This computation requires a file generated in the above Bulk Analysis sections. When writing the files, the summary table will be the one with the termination `*_assess,csv`, the beginning of the name depends on the users chosen `output_name`. This file is the one that can be fed to the `plot_categories` function.

`drugtax.plot_categories("testing_assess.csv", output_name = "plot")`

The `plot_categories` function has three arguments:
- `input_file`: the name of the `*_assess.csv` previously retrieved.
- `output_name`: a name for the output `*.png` file.
- `threshold`: with default 1, this function triggers an aggregation of low populated superclass combinations to their above counterpart. `threshold` is the minimum number of entries in the file for it to be aggregated. 