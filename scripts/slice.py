#!/usr/bin/env python

########################################################################
#                      slice.py - version 0.3                          #
#----------------------------------------------------------------------#
# Copyright (C) 2008, 2009 Matthew Dyer and Jonas Bjork                #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program.  If not, see http://www.gnu.org/licenses/.  #
#----------------------------------------------------------------------#
# A script to import volume data from a file with CHGCAR format and to #
# either interpolate the data onto a user specified 2D grid in an      #
# arbitrary plane or use the actual grid points.                       #
# The 2D grid data is then written to file.                            #
# Interpolation is performed using Fourier Transforms.                 #
# Alternatively 2D grid data can be imported from a file.              #
# Following either process a 2D contour plot can be made using pylab   #
# and matplotlib.                                                      #
#----------------------------------------------------------------------#
# The script can be run as an executable and it prompts the user for   #
# input in an interactive way. Alternatively input can be piped in to  #
# standard input using a command like:                                 #
#    script.py < input_file                                            #
#----------------------------------------------------------------------#
# The script depends on the following python modules and packages:     #
#    numpy - http://numpy.scipy.org/                                   #
#    ase - https://wiki.fysik.dtu.dk/ase/                              #
#    matplotlib - http://matplotlib.sourceforge.net/                   #
#----------------------------------------------------------------------#
# Authors: Matthew Dyer (msd30@liv.ac.uk), Jonas Bjork                 #
# Date of last revision: 19/01/09                                      #
########################################################################

import sys
import numpy
import time
from ase.calculators.vasp import VaspChargeDensity

print "----------------------------------------------------------------"
print "                     slice.py, version 0.3                      "
print "----------------------------------------------------------------"
print "Copyright (C) 2008, 2009 Matthew Dyer and Jonas Bjork."
print "This program is distributed under the GNU General Public Licence"
print "and comes with ABSOLUTELY NO WARRANTY; for details see the file"
print "GPL_LICENCE.txt or go to http://www.gnu.org/licenses/"
print "----------------------------------------------------------------"
starttime = time.clock() 
print "Starting calculation at",
print time.strftime("%H:%M:%S on %a %d %b %Y")
print "----------------------------------------------------------------"

# Should we import the volume data from a CHGCAR type of file or import 
# data on a 2D grid?
inputstr=raw_input("Import data from:\n %s %s %s" % ("1: A CHGCAR format file and interpolate\n",
                                                     "2: A CHGCAR format file and use actual grid\n",
                                                     "3: Data from a previous calculation\n"))
option=inputstr.strip()

#First option is to read data from a CHGCAR file and then interpolate density onto
#an arbitrary plane
if option=="1":
    # Ask which file to plot
    inputstr=raw_input("Enter filename of CHGCAR format file containing data:\n")
    
    print "Reading charge/potential data from file...",
    sys.stdout.flush()
    # Read density and atoms
    vasp_charge = VaspChargeDensity(filename = inputstr)
    density = vasp_charge.chg[-1]
    atoms = vasp_charge.atoms[-1]
    del vasp_charge
    print "Done.\n"

    # Read size of grid 
    ngridpts = numpy.array(density.shape)

    # Read total number of grid points
    totgridpts = ngridpts.prod()

    # Read scaling factor and unit cell
    unit_cell=atoms.get_cell()

    # Take Fourier transform of density
    print "Performing Fourier transform...",
    sys.stdout.flush()
    fdensity = numpy.fft.fftn(density)
    print "done."

    #Remove real space density from memory
    del density

    #Define the plane that we wish to get data for
    inputstr=raw_input("Enter origin in Cartesian coordinates for the plotting plane:\n")
    if len(inputstr.split())!=3:
        print "Syntax error. Quiting program."
        sys.exit(0)
    origin=[]
    for i in range(3):
        try:
            origin.append(float(inputstr.split()[i]))
        except:
            print "Syntax error. Quiting program."
            sys.exit(0)
    origin=numpy.array(origin)
    print "Origin of the plot is at %5f %5f %5f" % (origin[0],origin[1],origin[2])

    inputstr=raw_input("Enter the vector along the bottom of the plotting plane.\n Use Angstroms and the Cartesian axis system:\n") 
    if len(inputstr.split())!=3:
        print "Syntax error. Quiting program."
        sys.exit(0)
    xvector=[]
    for i in range(3):
        try:
            xvector.append(float(inputstr.split()[i]))
        except:
            print "Syntax error. Quiting program."
            sys.exit(0)
    xvector=numpy.array(xvector)
    xunit=xvector/numpy.sqrt(numpy.dot(xvector,xvector))
    print "Vector along bottom of plot is %5f %5f %5f" % (xvector[0],xvector[1],xvector[2])
    
    inputstr=raw_input("Enter the vector along the side of the plotting plane.\n Use Angstroms and the Cartesian axis system:\n") 
    if len(inputstr.split())!=3:
        print "Syntax error. Quiting program."
        sys.exit(0)
    yvector=[]
    for i in range(3):
        try:
            yvector.append(float(inputstr.split()[i]))
        except:
            print "Syntax error. Quiting program."
            sys.exit(0)
    yvector=numpy.array(yvector)
    yunit=yvector/numpy.sqrt(numpy.dot(yvector,yvector))
    print "Vector along side of plot is %5f %5f %5f" % (yvector[0],yvector[1],yvector[2])

    #Define the spacing of grid points
    inputstr=raw_input("Enter the number of grid points along bottom and side of plot:\n") 
    if len(inputstr.split())!=2:
        print "Syntax error. Quiting program."
        sys.exit(0)
    nplotgrid=[]
    for i in range(2):
        try:
            nplotgrid.append(int(inputstr.split()[i])+1)
        except:
            print "Syntax error. Quiting program."
            sys.exit(0)
    nplotgrid=numpy.array(nplotgrid)
    print "Number of grid points for plot is %dx%d" % (nplotgrid[0],nplotgrid[1])

    #Build array with Cartesian positions of grid points
    cartgrid=numpy.zeros((nplotgrid[0],nplotgrid[1],3),numpy.float)
    for j in range(nplotgrid[1]):
        for i in range(nplotgrid[0]):
            cartgrid[i][j]=origin+float(i)/(nplotgrid[0]-1)*xvector+float(j)/(nplotgrid[1]-1)*yvector

    #Perform Fourier Transform using array methods
    #Convert into scaled coordinates and move inside unit cell
    #Then multiply by 2*pi*i
    scaledgrid=[]
    for j in range(nplotgrid[1]):
        for i in range(nplotgrid[0]):
            scaledgrid.append(numpy.dot(cartgrid[i][j],numpy.linalg.inv(unit_cell))%1)
    scaledgrid=numpy.array(scaledgrid)
    scaledgrid=2.*numpy.pi*1.j*scaledgrid

    #Then build an array of reciprocal lattice vectors and
    #make a 1D array of the Fourier coefficients in fdensity
    gvectors=[]
    temp=[]
    for i in range(-ngridpts[0]/2+1,ngridpts[0]/2+1):
        for j in range(-ngridpts[1]/2+1,ngridpts[1]/2+1):
            for k in range(-ngridpts[2]/2+1,ngridpts[2]/2+1):
                gvectors.append([float(i),float(j),float(k)])
                temp.append(fdensity[i][j][k])
    gvectors=numpy.array(gvectors)
    gvectors=gvectors.T
    fdensity=numpy.array(temp)
    del temp

    density2D=numpy.zeros((nplotgrid[0]*nplotgrid[1]),numpy.float)
    #Loop over grid points to reduce memory use (probably slower than not doing)
    for i in range(nplotgrid[0]*nplotgrid[1]):
        #Build an array of the dot products between the grid points 
        #and the reciprocal lattice vectors and take the exponential
        temp=numpy.dot(scaledgrid[i],gvectors)
        temp=numpy.exp(temp)
        #Multiply by Fourier coefficients and sum
        density2D[i]=numpy.dot(fdensity,temp.T).real
    #Delete gvectors and fdensity to save space
    del temp
    del gvectors
    del fdensity

    #Take real part and divide by total number of grid points
    density2D=density2D/totgridpts

    #Reshape into 2D array for plotting
    density2D=density2D.reshape(nplotgrid[1],nplotgrid[0])

    #Work out x and y arrays for plotting - still needs some work!!
    #xvector will be plotted as the x-axis
    #Must be same dimensions as density2D
    #Find length of xvector and yvector
    xlength=numpy.sqrt(numpy.dot(xvector,xvector.T))
    ylength=numpy.sqrt(numpy.dot(yvector,yvector.T))
    #Find projection of yvector onto xvector
    yontox=numpy.dot(xvector,yvector.T)/xlength
    #Find component of yvector perpendicular to xvector
    ynormal=numpy.cross(xvector,yvector.T)/xlength
    ynormal=numpy.sqrt(numpy.dot(ynormal,ynormal.T))
    #Make arrays containing x and y values for each point
    xarray=numpy.zeros((nplotgrid[1],nplotgrid[0]),numpy.float)
    yarray=numpy.zeros((nplotgrid[1],nplotgrid[0]),numpy.float)
    for j in range(nplotgrid[1]):
        for i in range(nplotgrid[0]):
            xarray[j][i]=float(i)/float(nplotgrid[0]-1)*xlength+float(j)/float(nplotgrid[1]-1)*yontox
            yarray[j][i]=float(j)/float(nplotgrid[1]-1)*ynormal

    inputstr = raw_input("Enter name of file to write plot data to:\n")

    output_file=open(inputstr,"w")
    for j in range(nplotgrid[1]):
        for i in range(nplotgrid[0]):
            output_file.write("%20.10F %20.10F %20.10E\n" % (xarray[j][i],yarray[j][i],density2D[j][i]))
        output_file.write("\n")
    output_file.close()

    endtime = time.clock() 
    runtime = endtime-starttime
    print "\nEnd of calculation." 
    print "Program was running for %.2f seconds." % runtime
# Second option is to read density for a CHGCAR file and then to take the 
# density in a plane on the original grid
elif option=="2":
    # Ask which file to plot
    inputstr=raw_input("Enter filename of CHGCAR format file containing data:\n")
    # Read density and atoms
    print "Reading charge/potential data from file...",
    sys.stdout.flush()
    # Define the vasp charge density object
    # Read density and atoms
    vasp_charge = VaspChargeDensity(inputstr)
    density = vasp_charge.chg[-1]
    atoms = vasp_charge.atoms[-1]
    del vasp_charge
    print "Done.\n"
    # Read size of grid
    ngridpts = numpy.array(density.shape)
    print ngridpts

    # Read total number of grid points
    totgridpts = ngridpts.prod()

    # Read scaling factor and unit cell
    unit_cell=atoms.get_cell()

    #Define the plane that we wish to get data for
    inputstr=raw_input("Enter lattice vector not in the plane to plot (a, b or c):\n")
    normal=inputstr.strip().lower()[0]
    inormal = 'abc'.find(normal)
    if inormal==0:
       iplane1 = 1
       iplane2 = 2
    elif inormal==1:
       iplane1 = 0
       iplane2 = 2
    elif inormal==2:
       iplane1 = 0
       iplane2 = 1
    else:
       raise SyntaxError('Lattice vector must be either a, b, or c.')
#    if inormal==-1:
#        raise SyntaxError('Lattice vector must be either a, b, or c.')
#    iplane1 = (inormal+1)%3
#    iplane2 = (inormal+2)%3
    print "Vector not in plotting plane is %s" % normal

    inputstr=raw_input("Enter the distance in Angstroms along this vector to make the cut.\n")
    try:
        distance=float(inputstr.strip())
    except:
        print "Syntax error. Quiting program."
        sys.exit(0)
    print "Attempting to find a plane at a distance of %5f Angs" % distance

    #Find nearest plane
    #First calculate length of cell vectors
    cell_lengths=numpy.sqrt(numpy.dot(unit_cell,unit_cell.transpose()).diagonal())
    #Then find integer corresponding to closest plane on grid
    plane_index=int(round(ngridpts[inormal]*distance/cell_lengths[inormal]))%ngridpts[inormal]
    #Write out which distance we are actually using
    print "Using index %d which corresponds to a distance of %f Angstroms.\n" % (plane_index,float(plane_index)/float(ngridpts[inormal])*cell_lengths[inormal])

    #Cut out plane from 3D real space density
    if inormal==0:
       density2D=density[plane_index,:,:]
    elif inormal==1:
       density2D=density[:,plane_index,:].T
       density2D=density[:,plane_index,:]
    else:
       density2D=density[:,:,plane_index]

    #Make arrays of x and y values
    #First vector will be plotted as the x-axis
    #Must be same dimensions as density2D
    #Find projection of second vector onto first
    yontox=numpy.dot(unit_cell[iplane1],unit_cell[iplane2].T)/cell_lengths[iplane1]
    #Find component of yvector perpendicular to xvector
    ynormal=numpy.cross(unit_cell[iplane1],unit_cell[iplane2].T)/cell_lengths[iplane1]
    ynormal=numpy.sqrt(numpy.dot(ynormal,ynormal.T))
    #Make arrays containing x and y values for each point
    xarray=numpy.zeros((ngridpts[iplane1],ngridpts[iplane2]),numpy.float)
    yarray=numpy.zeros((ngridpts[iplane1],ngridpts[iplane2]),numpy.float)
    for i in range(ngridpts[iplane1]):
        for j in range(ngridpts[iplane2]):
            xarray[i][j]=float(i)/float(ngridpts[iplane1])*cell_lengths[iplane1]+float(j)/float(ngridpts[iplane2])*yontox
            yarray[i][j]=float(j)/float(ngridpts[iplane2])*ynormal

    #Write data to output file
    inputstr=raw_input("Enter name of file to write plot data to:\n") 

    output_file=open(inputstr,"w")
    for j in range(ngridpts[iplane1]):
        for i in range(ngridpts[iplane2]):
            output_file.write("%20.10F %20.10F %20.10E\n" % (xarray[j][i],yarray[j][i],density2D[j][i]))
        output_file.write("\n")
    output_file.close()

    endtime = time.clock() 
    runtime = endtime-starttime
    print "\nEnd of calculation." 
    print "Program was running for %.2f seconds." % runtime

elif option=="3":
    #Open input file
    inputstr=raw_input("Enter name of file containing plot data:\n")
    try:
       input_file=open(inputstr,"r")
    except:
       print "Error opening file.\n"
       sys.exit(0)

    #Read x-points and work out size of grid
    data=input_file.readlines()
    i=0
    while data[i]!="\n":
        i=i+1
    nplotgrid=[]
    nplotgrid.append(i)
    nplotgrid.append(len(data)/(i+1))
    nplotgrid=numpy.array(nplotgrid)

    #Read in y-points and density
    density2D=numpy.zeros((nplotgrid[1],nplotgrid[0]),numpy.float)
    xarray=numpy.zeros((nplotgrid[1],nplotgrid[0]),numpy.float)
    yarray=numpy.zeros((nplotgrid[1],nplotgrid[0]),numpy.float)
    k=0
    for j in range(nplotgrid[1]):
        for i in range(nplotgrid[0]):
            xarray[j][i]=float(data[k].split()[0])
            yarray[j][i]=float(data[k].split()[1])
            density2D[j][i]=float(data[k].split()[2])
            k=k+1
        k=k+1

    endtime = time.clock() 
    runtime = endtime-starttime
    print "\nEnd of calculation." 
    print "Program was running for %.2f seconds." % runtime
else:
    print "Syntax error. Quiting program."
    sys.exit(0)

#Have the option to plot in matplotlib
inputstr=raw_input("Would you like to plot the density using matplotlib? Yes / No:\n") 
if inputstr.strip().lower()=="yes":
    import pylab
    inputstr=raw_input("Enter number of contours for plot:\n") 
    try:
        ncontours=int(inputstr.strip())
    except:
        print "Syntax error. Quiting program."
        exit(0)
    inputstr=raw_input("Would you like to fill the contours in the plot? Yes / No:\n") 
    if inputstr.strip().lower()=="yes":
       density_plot = pylab.contourf(xarray,yarray,density2D,ncontours,antialiased=False)
    else:
       density_plot = pylab.contour(xarray,yarray,density2D,ncontours,antialiased=False)
    inputstr=raw_input("Would you like a colorbar in the plot? Yes / No:\n") 
    if inputstr.strip().lower()=="yes":
       pylab.colorbar(extend="both",format="%.3E")
    pylab.axis('image')
    pylab.show()
