# -*- coding: utf-8 -*-
"""
Created on Thu Dec 16 12:39:23 2021

@author: hbalsche

Manual: set file name ..................... in line 26
        set guessed diameter .............. in line 60
        set rotation about cylinder axis .. in line 61
        set limit value of color bands .... in line 72
        call desired plotting fuctions .... in line 74 to 82 
"""

import subMods.plotting as pl
import subMods.zylFit as zylFit
import pathlib 


##############################################################################
#                                                                            #
# file name and paths                                                        #
#                                                                            #
##############################################################################

# name of stl file in folder 01_STL_Daten
stlFilename = 'MeshSample.stl'

# define directory names:
currentDir = pathlib.Path.cwd()
home = currentDir.parents[0]
stlDir = home/'01_STL_Daten'
imgDir = home/'03_Diagramme'


##############################################################################
#                                                                            #
# run assessment                                                             #
#                                                                            #
##############################################################################

stleFilePath = stlDir/stlFilename
zylFitter = zylFit.fitter()   
zylFitter.readStl(stleFilePath)

# guess axes position by center of bounding box
x_max, x_min = max(zylFitter.data[0]), min(zylFitter.data[0])
y_max, y_min = max(zylFitter.data[1]), min(zylFitter.data[1])
z_max, z_min = max(zylFitter.data[2]), min(zylFitter.data[2]) 
 
x0 = 0.5*(x_min+x_max)
y0 = 0.5*(y_min+y_max)
z0 = 0.5*(z_min+z_max)
 
dx = (x_min-x_max)
dy = (y_min-y_max)
dz = (z_min-z_max)

# compute best for cylinder 

diam_guess = 200
rotation   = 155
zylFitter.computeBestFitZylinder(diam_guess, [x0,y0,z0], [dx,dy,dz], 
                                 rot_z_deg = rotation)#, reverse = True)

#zylFitter.removePointsBetween(1500, 2500)

##############################################################################
#                                                                            #
# draw sketches                                                              #
#                                                                            #
##############################################################################

delBot = -5  # bottom limit of color band 
delTop =  5  # top limit of color band
sketcher = pl.plotter(zylFitter.getPlottingInputs())


sketcher.plotFitting(base = 'guess', minTol = delBot, maxTol=delTop, nStep=10, 
                     fPath = imgDir)
sketcher.plotUnwrapped(base = 'guess', minTol = delBot, maxTol=delTop, 
                       fPath = imgDir, nStep = 10, reverseSamples=True)
sketcher.SurfPlotUnwrapped(nPhi = 200, nz = 50)