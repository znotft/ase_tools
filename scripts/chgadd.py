#!/usr/bin/env python
"""
A script which reads CHGCAR files and adds the 
densities together.
Output is to a file called CHGADD.
Depends on ase.
"""

import os
import sys
import numpy
import time
from ase.calculators.vasp import VaspChargeDensity

starttime = time.clock() 
print "Starting calculation at",
print time.strftime("%H:%M:%S on %a %d %b %Y")

# Find out how many arguments were on the command line, 
# all but the first two should contain files with densities
# to subtract
nsubtract = len(sys.argv)-2
if not nsubtract >= 1:
    print "\n** ERROR: Must specify name of at least two files on command line."
    print "eg. chgadd.py CHGCAR1 CHGCAR2  [CHGCAR3 ...]"
    print "The reference density is taken from the first filename."                                            
    print "The densities in the files after this will be subtracted from the reference."                       
    sys.exit(0)

# Check that files exist
for name in sys.argv[1:]:
    if not os.path.isfile(name):
        print "\n** ERROR: Input file %s was not found." % name
        sys.exit(0)

# Read information from command line
# First specify location of CHGCAR file with reference density
CHGCARfile1 = sys.argv[1].lstrip()

# Open geometry and density class objects
#-----------------------------------------
print "Reading density data from file %s ..." % CHGCARfile1,
sys.stdout.flush()
vasp_charge1 = VaspChargeDensity(filename = CHGCARfile1)
chg1 = vasp_charge1.chg[-1]
atoms1 = vasp_charge1.atoms[-1]
del vasp_charge1
print "done." 

chgadd=chg1
for CHGCARfile2 in sys.argv[2:]:
    CHGCARfile2 = CHGCARfile2.strip()
    print "Reading density data from file %s ..." % CHGCARfile2,
    sys.stdout.flush()
    vasp_charge2 = VaspChargeDensity(filename = CHGCARfile2)
    chg2 = vasp_charge2.chg[-1]
    del vasp_charge2
    print "done." 

    # Make sure that the second data set is on the same grid
    #--------------------------------------------------------
    if chg2.shape != chg1.shape:
       print "\n**ERROR: Two sets of data are not on the same grid."
       print "Data from file %s on %dx%dx%d grid." % (CHGCARfile1,chg1.shape[0],chg1.shape[1],chg1.shape[2])
       print "Data from file %s on %dx%dx%d grid.\n" % (CHGCARfile2,chg2.shape[0],chg2.shape[1],chg2.shape[2])
       sys.exit(0)
    else:
       print "Adding data from file %s ..." % CHGCARfile2,
       sys.stdout.flush()

    # Add charge density
    #-----------------
    chgadd += chg2
    print "done." 

zero = raw_input("Set negative values of the added charge density to zero (Yes/No): ")

vasp_charge_add = VaspChargeDensity(filename=None)
vasp_charge_add.atoms=[atoms1,]
vasp_charge_add.chg=[chgadd,] 

# Print out charge density
#--------------------------
# Check whether CHGADD exists
if os.path.isfile("./CHGADD"):
   print "\n**WARNING: A file called CHGADD already exists in this directory."
   yesno=raw_input("Type y to continue and overwrite it, any other key to stop\n")
   if yesno!="y":
      sys.exit(0)


print "Writing added density data to file CHGADD ...", 
sys.stdout.flush()
vasp_charge_add.write(filename="CHGADD",format="chgcar")
print "done."

endtime = time.clock() 
runtime = endtime-starttime
print "\nEnd of calculation." 
print "Program was running for %.2f seconds." % runtime
