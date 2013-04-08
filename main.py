'''
Created on Mar 19, 2013

@author: tristan
'''

import sys
import numpy as np
from meshBuilder import *
from spacialDiscretization import *
from modelChecker import *
from timestepper import *
from fluxCalculations import *
from sourceCalculations import *
from spacialDiscretization import *
from boundaryConditions import *
from dataSaver import *

np.seterr(all='raise')

# Directory for saving fort.14 and fort.63 files
saveOutput = True
workingDir = '/home/tristan/Desktop/'

# Define all constants here
gridSize = 20
gridPadding = 2
cellWidth = 1.0  # meters
cellHeight = 1.0  # meters
floorElevation = 1.0  # meters above 0.0
waterDepth = 3.0  # meters
manningsN = 0.03

runTime = 10.0  # simulated run time in seconds
currTime = 0.0  # starting time
iterationCount = 0  # cumulative number of timesteps


if __name__ == "__main__":

    #################################################
    # Build mesh -- Create global constant variables as well
    #               as the arrays that will be used to store U
    #               data
    # ------------------------------------------------------------
    # meshCoordinates - [i][j][bottomLeft(x, y, z)]    (m+1)x(n+1)
    # meshBottomIntegrationPoints - [i][j][zS, zW]     (m+1)x(n+1)
    # meshBottomSlopes - [i][j][dB/dx, dB/dy]          (m)x(n)
    # meshBottomCenters - [i][j][z]                    (m)x(n)
    #
    # meshU - [i][j][w, hu, hv]
    # meshUstar - [i][j][w, hu, hv]
    # meshUnext - [i][j][w, hu, hv]
    #################################################


    ###### Flat bottom ######
    # meshCoordinates, meshBottomIntegrationPoints, meshBottomSlopes, meshBottomCenters = buildTestMesh(gridSize, gridPadding, cellWidth, cellHeight, floorElevation)
    # Constant free surface
    # meshU = buildConstantTestU(meshCoordinates, meshBottomCenters, floorElevation + waterDepth)
    # Pyramid free surface
    # meshU = buildPyramidTestU(meshCoordinates, meshBottomCenters, floorElevation + waterDepth, floorElevation + 1.0, 8)
    #########################

    ###### Sloping bottom with beach ######
    # meshCoordinates, meshBottomIntegrationPoints, meshBottomSlopes, meshBottomCenters = buildSlopingTestMesh(gridSize, gridPadding, cellWidth, cellHeight, 1.0, 5.0)
    # meshU = buildConstantTestU(meshCoordinates, meshBottomCenters, 10.0)
    #######################################

    ###### Partial sloping bottom #########
    meshCoordinates, meshBottomIntegrationPoints, meshBottomSlopes, meshBottomCenters = buildDoubleSlopingTestMesh(gridSize, gridPadding, cellWidth, cellHeight, 1.0, 0.1, 12)
    # meshU = buildConstantTestU(meshCoordinates, meshBottomCenters, 15.0)
    meshU = buildPyramidTestU(meshCoordinates, meshBottomCenters, 5.0, 4.5, 8)
    #######################################

    meshUstar = np.zeros((meshU.shape[0], meshU.shape[1], 3))
    meshUnext = np.zeros((meshU.shape[0], meshU.shape[1], 3))

    if saveOutput:
        writeCustomFort14(workingDir, meshCoordinates)
        fort63 = createFort63(workingDir, meshCoordinates)


    # print("Starting serial timestepping...\n")

    # Timestepping begins here
    while (currTime < runTime):

        # Save timestep data to fort.63
        if saveOutput:
            writeCustomTimestep(fort63, meshU)
            # writeTimestep(fort63, meshU, meshBottomCenters)

        # printMatrix(meshU)
        print meshU[0][4], meshU[1][4], meshU[2][4], meshU[3][4]

        ##### Building of meshUstar begins here #####
        # First calculate the values of U at integration points and the one-sided propagation speeds
        meshUIntegrationPoints, meshOneSidedPropagationSpeeds = reconstructFreeSurface(meshCoordinates, meshU, meshBottomIntegrationPoints, gridPadding, cellWidth, cellHeight)
        # print "a+: " + str(meshOneSidedPropagationSpeeds[1][27][2]) + "\ta-: " + str(meshOneSidedPropagationSpeeds[2][27][3])
        # meshUIntegrationPoints, meshOneSidedPropagationSpeeds = newReconstructFreeSurface(meshU, meshBottomCenters, meshBottomIntegrationPoints, cellWidth, cellHeight)

        # printMatrix(meshOneSidedPropagationSpeeds)

        # Next calculate fluxes and bed slope source terms
        meshFluxes = newFluxSolver(meshUIntegrationPoints, meshBottomIntegrationPoints, meshOneSidedPropagationSpeeds)
        meshSlopeSource = bedSlopeSourceSolver(meshUIntegrationPoints, meshBottomIntegrationPoints, cellWidth, cellHeight)

        # Calculate dt
        dt = calculateTimestep(meshOneSidedPropagationSpeeds, cellWidth, cellHeight)
        sys.stdout.write("Iteration " + str(iterationCount) + "\tSimulated Time: " + "%.6f" % currTime + "\tdt: " + \
        	"%.6f" % dt + "\tSimulated Time Remaining: " + "%.6f" % (runTime - currTime) + "\r")
        sys.stdout.flush()

        # Calculate R and bed shear source terms
        meshR = buildRValues(meshSlopeSource, meshFluxes)
        # meshShearSource = bedShearSourceSolver(meshU, meshBottomCenters, manningsN)
        meshShearSource = newBedShearSourceSolver(meshU, meshBottomCenters, manningsN, cellWidth, cellHeight)

        # Assemble meshUstar
        meshUstar = buildU_star(meshUstar, meshU, dt, meshR, meshShearSource)




        ##### Apply boundary conditions
        meshUstar = applyWallBoundaries(meshUstar, gridSize, gridPadding)
        # print 'Ustar:'
        # printMatrix(meshUstar)
        # print meshUstar[50][26], meshUstar[51][26], meshUstar[52][26], meshUstar[53][26]


        ##### Building of meshUnext begins here
        # First calculate values of U at integration points and the one-sided propagation speeds
        meshUIntegrationPoints, meshOneSidedPropagationSpeeds = reconstructFreeSurface(meshCoordinates, meshUstar, meshBottomIntegrationPoints, gridPadding, cellWidth, cellHeight)
        # meshUIntegrationPoints, meshOneSidedPropagationSpeeds = newReconstructFreeSurface(meshUstar, meshBottomCenters, meshBottomIntegrationPoints, cellWidth, cellHeight)

        # Next calculate fluxes and bed slope source terms
        meshFluxes = newFluxSolver(meshUIntegrationPoints, meshBottomIntegrationPoints, meshOneSidedPropagationSpeeds)
        meshSlopeSource = bedSlopeSourceSolver(meshUIntegrationPoints, meshBottomIntegrationPoints, cellWidth, cellHeight)

        # Calculate R and bed shear source terms
        meshR = buildRValues(meshSlopeSource, meshFluxes)
        # meshShearSource = bedShearSourceSolver(meshUstar, meshBottomCenters, manningsN)
        meshShearSource = newBedShearSourceSolver(meshUstar, meshBottomCenters, manningsN, cellWidth, cellHeight)

        # Assemble meshUnext
        meshUnext = buildU_next(meshUnext, meshU, meshUstar, dt, meshR, meshShearSource)
        # print "Unext:"
        # printMatrix(meshUnext)

        ##### Set meshU equal to the next timestep U with boundaries applied
        meshU = applyWallBoundaries(meshUnext, gridSize, gridPadding)



        # Print results
        # printMatrix(meshUnext)

        ##### Increment the timestep
        currTime += dt
        iterationCount += 1


    if saveOutput:
        # closeFort63(workingDir, fort63, meshCoordinates, iterationCount)
        closeCustomFort63(workingDir, fort63, meshU, iterationCount)
        print "Finished. " + str(iterationCount) + " timesteps recorded"
    else:
        print "Finished."
