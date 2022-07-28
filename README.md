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
