![Eco](header.png)


Primary Devloper: John Docter
Please reach out to john.docter@colorado.edu with any questions 

Introduction
============
EcoGenoRisk a risk assessment tool designed to assess the risk a synthetic organism can introduce into the environment. 

EcoGenoRisk is composed of HazID, EnvCen, and RiskQ(under development).

**HazID**
=========
Utilizes DIAMOND sequence aligner to run a functional comparison of a synbio organism against all known organisms. 

_1_ through _5_: Outputs the most susceptible organism to competition. 

_6_Competitor_Find.py:  Outputs the InchiKeys that are the shared reactants for different enzymes, i.e. what reactants the topmatch organism and synbio organism will be competing for. 

_7_SPLENDA.py: Takes in the JSON file of BRENDA's inhibitor list per enzyme and converts it to a csv for use in _7_PoisInhibitor 
_7_PoisInhibitor.py: Takes in a list of all Inhibitors and their EC numbers and outputs a list of species and EC's that are inhibitied by a product of the synbio org and what inhibits them 

**EnCen**
=========
Biome specific functional analysis

Runs on DIAMOND and JGI to identify most susceptible biome to invasion across multiple biome inputs 

Allows for deep dive into bin lineage and enzymatic effect of synbio invasion on any given environment 

