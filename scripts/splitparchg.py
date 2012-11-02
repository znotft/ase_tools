#!/usr/bin/env python
"""
A script which reads VASP charge density files from spin polarised calculations
and outputs three new files containing the majority spin data, the minority spin
data and the difference between the two.
Depends on ase
"""

import os
import sys
import time
import numpy
from ase.calculators.vasp import VaspChargeDensity

starttime = time.clock() 
print "Starting calculation at",
print time.strftime("%H:%M:%S on %a %d %b %Y")

if len(sys.argv) != 2:
    print "\n** ERROR: Must specify name of file on command line."
    print "eg. splitparchg.py PARCHG ."                                         
    sys.exit(0)

if not os.path.isfile(sys.argv[1]):
    print "\n** ERROR: Input file %s was not found." % sys.argv[1]
    sys.exit(0)

# Read information from command line
# First specify location of PARCHG 
PARCHGfile = sys.argv[1].lstrip()

# Open geometry and density class objects
#-----------------------------------------
print "Reading potential data from file %s ..." % PARCHGfile,
sys.stdout.flush()
vasp_charge_data = VaspChargeDensity(filename=PARCHGfile)
print "done." 
# Check data is spin polarised
if not vasp_charge_data.is_spin_polarized():
    print "\n** ERROR: Input file does not contain spin polarised data."
    sys.exit(0)
# Make Atoms object and arrays of density data
geomdata = vasp_charge_data.atoms[-1]
parchg_sum = vasp_charge_data.chg[-1]
parchg_diff = vasp_charge_data.chgdiff[-1]

# Read in potential data
#------------------------
ngridpts = numpy.array(parchg_sum.shape)
totgridpts = ngridpts.prod()
print "Potential stored on a %dx%dx%d grid" % (ngridpts[0],ngridpts[1],ngridpts[2])
print "Total number of points is %d" % totgridpts

# Calculate up spin density (sum+diff)/2
#----------------------------------------
parchg_up = (parchg_sum+parchg_diff)/2.

# Calculate down spin density (sum-diff)/2
#----------------------------------------
parchg_down = (parchg_sum-parchg_diff)/2.

# Write out parchg files
#------------------------
parchgfile = PARCHGfile + ".UP"
print "Writing up spin data to file %s..." % parchgfile,
sys.stdout.flush()
output_charge_data = VaspChargeDensity(filename=None)
output_charge_data.atoms=[geomdata,]
output_charge_data.chg=[parchg_up,]
output_charge_data.write(filename=parchgfile,format="chgcar")
print "done."

parchgfile = PARCHGfile + ".DOWN"
print "Writing down spin data to file %s..." % parchgfile,
sys.stdout.flush()
output_charge_data = VaspChargeDensity(filename=None)
output_charge_data.atoms=[geomdata,]
output_charge_data.chg=[parchg_down,]
output_charge_data.write(filename=parchgfile,format="chgcar")
print "done."

parchgfile = PARCHGfile + ".DIFF"
print "Writing up-down spin data to file %s..." % parchgfile,
sys.stdout.flush()
output_charge_data = VaspChargeDensity(filename=None)
output_charge_data.atoms=[geomdata,]
output_charge_data.chg=[parchg_diff,]
output_charge_data.write(filename=parchgfile,format="chgcar")
print "done."

endtime = time.clock() 
runtime = endtime-starttime
print "\nEnd of calculation." 
print "Program was running for %.2f seconds." % runtime
