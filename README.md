Welcome to Cancer_ABM_Project,

Here we simultate HELA Cells (or any other kind of cells) with the use of Biodynamo ABM and we tune the models using Bayesian Optimization or Genetic Algorithm.

Code files included:

Bayesian optimaziation code generator python code 

Genetic algorithm code generator python code

Python code to generate the target CSV files from Data XLS files

Python code to generate graphs

Genreated Python model code for both Bayesian optimization and Genetic algorithm (This is just included for refrence only and will not work. Please generate your own code instead)

Model configuration:

Model configuration is included in JSON files with example models included in the repo

How to use the example configurations to generate models:

In the example configuration JSON files replace "/your/dir" with your directory path and also include a correct and existing path for the results storage directory "results_dir" as well as for the
directory that the target input files are stored ex. "Target_Treated" 

After finishing the configuration of your JSON file go to the python code generator file of your model that you wish to generate (Bayesian Optimization or Genetic Algorithm) and replace the path 
in this line with the location of your JSON file.

f=open('/your/dir/Your_Model_Example_Conf.json')


Congrats, the model is generated and ready to be used.



PS. Biodynamo instalation is required to run any model

Reference here: https://github.com/BioDynaMo/biodynamo
