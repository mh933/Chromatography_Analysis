# Chromatography_Analysis

This python/Jupyter Notebook program has been produced to create a pipeline in which input data, from the program Unicorn 7 (Cytivia), can be automatically analysed to a visually appealing figure. This pipeline will also be used to calculate the protein concentration in collected fractions to aid in analysis of purification success. Currently this data only works with specific files where fractions are labelled A1-A13 and with a 5 mL sample loop.

To run:

  1.) Download AdkAnalysis.py file. This is the raw python code.
  2.) Download AdkAnalysisNotebook.ipnb file. This is a Jupyter Notebook file, in which is a module for Adk Affinity Chromatography Analysis. This is more user-friendly than using       the raw python code.
  3.) Download the example data (Adk1 500mL Culture Data.csv and Adk4 Arctic Express 500mL Culture Data.csv) AND control data (Adk 4 Blank Data Example.csv).
  4.) Open AdkAnalysisNotebook.ipnb in Jupyter Notebook
  5.) SHIFT+ENTER on the first code cell: this will run the code and define a function Adk_aff_chrom_analysis()

Adk_aff_chrom_analysis will require 2 arguments to run: Adk_aff_chrom_analysis("Input_Data_File_Name", AU_Peak_Minumum)
  > The first argument is the input file: this needs to be in the format of .csv and must be enclosed by quotation marks/apostrophes e.g. "Adk1 500mL Culture Data.csv"
  > The second argument is the minimum AU value at which a peak will be detected. A smaller number will be more sensitive than a greater number
  
For example:
  
  Adk_aff_chrom_analysis("Adk1 500mL Culture Data.csv", 350)
  
  OR:
  
  Adk_aff_chrom_analysis("Adk4 Arctic Express 500mL Culture Data.csv", 250)
  
  (Both examples have been included at the bottom of the .ipnb file
  
