# socorro2cif
Converts Socorro crystal and new_crystal files to CIF format.

## Instructions
Run socorro2cif.py in directory containing crystal file. It will write this crystal in cif format to crystal_0.cif. If there is a new_crystal file present, it will convert this into multiple cif files (crystal_1.cif, crystal_2.cif, ..., crystal_N.cif).

This has not been ported to python 3, although I suspect only the print statements need to be changed or removed.
