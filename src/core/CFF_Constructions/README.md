This folder is for individual constructions of CFFs (and the finite field wrapper)

To add a new construction, make a new .c and .h that implements:

-A CFF_Construction_And_Name_Functions

-A function to fill in the tables that will be called in the function makeTables() in CFF_Tables.c, and add the new function to makeTables()

(see Doubling_Construction for a good example)