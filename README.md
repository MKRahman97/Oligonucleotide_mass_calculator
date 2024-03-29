﻿# Oligonucleotide_mass_calculator

To run the script ensure that all the files are in the same directory. Define the sequence in myconfig.py, as well as the 3' and 5' ends. Modifications to the sugar ring and phosphate linker can also be defined here. A list of neutral losses is needed and then selected at the bottom if they are to be calculated. Mcluckey cleavages, base losses, and internal fragments are also defined here. 

In order to calculate the masses you need to execute the following command:
python Mass_calc.py

This will read the parameters defined in myconfig

Adjust the sequence as required in order to calculate the McLuckey cleavages. Modifications are defined in the Nucleic_acids.csv file, which can be customised as needed. Note the character 'd' is reserved for locked nucleic acid modifications only. So far, the common oligonucleotide elements have been noted in the Element_masses.csv file, which can be expanded when needed.
