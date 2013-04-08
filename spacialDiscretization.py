'''
Created on Mar 20, 2013

@author: tristan
'''
import numpy as np
from mesh.modelChecker import print4dMatrix, printMatrix
import sys


def newReconstructFreeSurface(meshU, meshBottomCenters, meshBottomIntegrationPoints, cellWidth, cellHeight):

    # These sizes include all padding
    m = meshU.shape[0]
    n = meshU.shape[1]

    # Constants
    theta = 1.3
    g = 9.81
    sqrt2 = np.sqrt(2.0)
    K0 = 0.01
    Kappa = K0 * max(1.0, min(cellWidth, cellHeight))

    Q = np.zeros((m, n, 2, 3))  # Slopes at cell centers
    meshUIntegrationPoints = np.zeros((m, n, 4, 3))  # U at the integration points inside the cell
    huv = np.zeros((m, n, 4, 3))  # h, u, v at integration points inside the cell
    meshOneSidedPropagationSpeeds = np.zeros((m, n, 4))  # N, S, E, W propagation speeds inside each cell


    for i in range(m):
        for j in range(n):
            if (meshU[i][j][0] <= meshBottomCenters[i][j]):
                meshU[i][j][0] = meshBottomCenters[i][j]
                meshU[i][j][1] = 0.0
                meshU[i][j][2] = 0.0
                meshUIntegrationPoints[i][j][0] = np.array([meshBottomIntegrationPoints[i][j + 1][0], 0.0, 0.0])
                meshUIntegrationPoints[i][j][1] = np.array([meshBottomIntegrationPoints[i][j][0], 0.0, 0.0])
                meshUIntegrationPoints[i][j][2] = np.array([meshBottomIntegrationPoints[i + 1][j][1], 0.0, 0.0])
                meshUIntegrationPoints[i][j][3] = np.array([meshBottomIntegrationPoints[i][j][1], 0.0, 0.0])
                huv[i][j] = np.zeros((4, 3))


    for i in range(1, m - 1):
        for j in range(1, n - 1):
            # Calculate free surface slopes for wet cells only
            if (meshU[i][j][0] > meshBottomCenters[i][j]):

                # ## Calculate d/dx
                # dw/dx
                forwardX = (meshU[i + 1][j][0] - meshU[i][j][0]) / cellWidth
                centralX = (meshU[i + 1][j][0] - meshU[i - 1][j][0]) / (2 * cellWidth)
                backwardX = (meshU[i][j][0] - meshU[i - 1][j][0]) / cellWidth
                Q[i][j][0][0] = minmod(theta * forwardX, centralX, theta * backwardX)
                # dhu/dx
                forwardX = (meshU[i + 1][j][1] - meshU[i][j][1]) / cellWidth
                centralX = (meshU[i + 1][j][1] - meshU[i - 1][j][1]) / (2 * cellWidth)
                backwardX = (meshU[i][j][1] - meshU[i - 1][j][1]) / cellWidth
                Q[i][j][0][1] = minmod(theta * forwardX, centralX, theta * backwardX)
                # dhv/dx
                forwardX = (meshU[i + 1][j][2] - meshU[i][j][2]) / cellWidth
                centralX = (meshU[i + 1][j][2] - meshU[i - 1][j][2]) / (2 * cellWidth)
                backwardX = (meshU[i][j][2] - meshU[i - 1][j][2]) / cellWidth
                Q[i][j][0][2] = minmod(theta * forwardX, centralX, theta * backwardX)

                # ## Calculate d/dy
                # dw/dy
                forwardY = (meshU[i][j + 1][0] - meshU[i][j][0]) / cellHeight
                centralY = (meshU[i][j + 1][0] - meshU[i][j - 1][0]) / (2 * cellHeight)
                backwardY = (meshU[i][j][0] - meshU[i][j - 1][0]) / cellHeight
                Q[i][j][1][0] = minmod(theta * forwardY, centralY, theta * backwardY)
                # dhu/dy
                forwardY = (meshU[i][j + 1][1] - meshU[i][j][1]) / cellHeight
                centralY = (meshU[i][j + 1][1] - meshU[i][j - 1][1]) / (2 * cellHeight)
                backwardY = (meshU[i][j][1] - meshU[i][j - 1][1]) / cellHeight
                Q[i][j][1][1] = minmod(theta * forwardY, centralY, theta * backwardY)
                # dhv/dy
                forwardY = (meshU[i][j + 1][2] - meshU[i][j][2]) / cellHeight
                centralY = (meshU[i][j + 1][2] - meshU[i][j - 1][2]) / (2 * cellHeight)
                backwardY = (meshU[i][j][2] - meshU[i][j - 1][2]) / cellHeight
                Q[i][j][1][2] = minmod(theta * forwardY, centralY, theta * backwardY)


    for i in range(1, m - 1):
        for j in range(1, n - 1):
            # Calculate U at integration points for wet cells
            if (meshU[i][j][0] > meshBottomCenters[i][j]):
                meshUIntegrationPoints[i][j][0] = meshU[i][j] + (cellHeight / 2.0) * Q[i][j][1]  # North = center + dist*d/dy
                meshUIntegrationPoints[i][j][1] = meshU[i][j] - (cellHeight / 2.0) * Q[i][j][1]  # South = center - dist*d/dy
                meshUIntegrationPoints[i][j][2] = meshU[i][j] + (cellWidth / 2.0) * Q[i][j][0]  # East = center + dist*d/dx
                meshUIntegrationPoints[i][j][3] = meshU[i][j] - (cellWidth / 2.0) * Q[i][j][0]  # West = center - dist*d/dx

                huv[i][j][0][0] = meshUIntegrationPoints[i][j][0][0] - meshBottomIntegrationPoints[i][j + 1][0]  # h - North
                huv[i][j][1][0] = meshUIntegrationPoints[i][j][1][0] - meshBottomIntegrationPoints[i][j][0]  # h - South
                huv[i][j][2][0] = meshUIntegrationPoints[i][j][2][0] - meshBottomIntegrationPoints[i + 1][j][1]  # h - East
                huv[i][j][3][0] = meshUIntegrationPoints[i][j][3][0] - meshBottomIntegrationPoints[i][j][1]  # h - West
                huv[i][j][0][1] = (sqrt2 * huv[i][j][0][0] * meshUIntegrationPoints[i][j][0][1]) / (np.sqrt(huv[i][j][0][0] ** 4 + max(huv[i][j][0][0] ** 4, Kappa)))  # u - North
                huv[i][j][1][1] = (sqrt2 * huv[i][j][1][0] * meshUIntegrationPoints[i][j][1][1]) / (np.sqrt(huv[i][j][1][0] ** 4 + max(huv[i][j][1][0] ** 4, Kappa)))  # u - South
                huv[i][j][2][1] = (sqrt2 * huv[i][j][2][0] * meshUIntegrationPoints[i][j][2][1]) / (np.sqrt(huv[i][j][2][0] ** 4 + max(huv[i][j][2][0] ** 4, Kappa)))  # u - East
                huv[i][j][3][1] = (sqrt2 * huv[i][j][3][0] * meshUIntegrationPoints[i][j][3][1]) / (np.sqrt(huv[i][j][3][0] ** 4 + max(huv[i][j][3][0] ** 4, Kappa)))  # u - West
                huv[i][j][0][2] = (sqrt2 * huv[i][j][0][0] * meshUIntegrationPoints[i][j][0][2]) / (np.sqrt(huv[i][j][0][0] ** 4 + max(huv[i][j][0][0] ** 4, Kappa)))  # v - North
                huv[i][j][1][2] = (sqrt2 * huv[i][j][1][0] * meshUIntegrationPoints[i][j][1][2]) / (np.sqrt(huv[i][j][1][0] ** 4 + max(huv[i][j][1][0] ** 4, Kappa)))  # v - South
                huv[i][j][2][2] = (sqrt2 * huv[i][j][2][0] * meshUIntegrationPoints[i][j][2][2]) / (np.sqrt(huv[i][j][2][0] ** 4 + max(huv[i][j][2][0] ** 4, Kappa)))  # v - East
                huv[i][j][3][2] = (sqrt2 * huv[i][j][3][0] * meshUIntegrationPoints[i][j][3][2]) / (np.sqrt(huv[i][j][3][0] ** 4 + max(huv[i][j][3][0] ** 4, Kappa)))  # v - West

                meshUIntegrationPoints[i][j][:][1] = huv[i][j][:][0] * huv[i][j][:][1]  # Calculate new values of hu for N, S, E, W
                meshUIntegrationPoints[i][j][:][2] = huv[i][j][:][0] * huv[i][j][:][2]  # Calculate new values of hv for N, S, E, W

    print "huv:"
    print4dMatrix(huv, 2, 2)
    for i in range(1, m - 1):
        for j in range(1, n - 1):
            # Calculate North propagation speed of current element
            meshOneSidedPropagationSpeeds[i][j][0] = min(huv[i][j][0][2] - np.sqrt(g * huv[i][j][0][0]),
                                                         huv[i][j + 1][1][2] - np.sqrt(g * huv[i][j + 1][1][0]),
                                                         0.0)

            # Calculate South propagation speed of element above current one
            meshOneSidedPropagationSpeeds[i][j + 1][1] = max(huv[i][j][0][2] + np.sqrt(g * huv[i][j][0][0]),
                                                           huv[i][j + 1][1][2] + np.sqrt(g * huv[i][j + 1][1][0]),
                                                           0.0)

            # Calculate East propagation speed of current element
            meshOneSidedPropagationSpeeds[i][j][2] = min(huv[i][j][2][1] - np.sqrt(g * huv[i][j][2][0]),
                                                         huv[i + 1][j][3][1] - np.sqrt(g * huv[i + 1][j][3][0]),
                                                         0.0)

            # Calculate West propagation speed of the element to the right
            meshOneSidedPropagationSpeeds[i + 1][j][3] = max(huv[i][j][2][1] + np.sqrt(g * huv[i][j][2][0]),
                                                           huv[i + 1][j][3][1] + np.sqrt(g * huv[i + 1][j][3][0]),
                                                           0.0)

    return meshUIntegrationPoints, meshOneSidedPropagationSpeeds



def reconstructFreeSurface(meshCoordinates, meshU, meshBottomIntegrationPoints, paddingSize, cellWidth, cellHeight):

    m = meshCoordinates.shape[0]
    n = meshCoordinates.shape[1]

    # Constants
    theta = 1.3
    g = 9.81
    sqrt2 = np.sqrt(2.0)
    K0 = 0.01
    Kappa = K0 * max(1.0, min(cellWidth, cellHeight))

    Q = np.zeros((m, n, 6))  # [i][j][dw/dx, dw/dy, dhu/dx, dhu/dy, dhv/dx, dhv/dy]
    forwardX = 0.0  # forward approximation of slope in x-direction
    forwardY = 0.0  # forward approximation of slope in y-direction
    centralX = 0.0  # central approximation of slope in x-direction
    centralY = 0.0  # central approximation of slope in y-direction
    backwardX = 0.0  # backward approximation of slope in x-direction
    backwardY = 0.0  # backward approximation of slope in y-direction

    for i in range(1, m - 2):
        for j in range(1, n - 2):

            # Calculate free surface slopes
            forwardX = (meshU[i + 1][j][0] - meshU[i][j][0]) / cellWidth
            centralX = (meshU[i + 1][j][0] - meshU[i - 1][j][0]) / (2 * cellWidth)
            backwardX = (meshU[i][j][0] - meshU[i - 1][j][0]) / cellWidth

            forwardY = (meshU[i][j + 1][0] - meshU[i][j][0]) / cellHeight
            centralY = (meshU[i][j + 1][0] - meshU[i][j - 1][0]) / (2 * cellHeight)
            backwardY = (meshU[i][j][0] - meshU[i][j - 1][0]) / cellHeight

            Q[i][j][0] = minmod(theta * forwardX, centralX, theta * backwardX)
            Q[i][j][1] = minmod(theta * forwardY, centralY, theta * backwardY)

            # Calculate hu slopes
            forwardX = (meshU[i + 1][j][1] - meshU[i][j][1]) / cellWidth
            centralX = (meshU[i + 1][j][1] - meshU[i - 1][j][1]) / (2 * cellWidth)
            backwardX = (meshU[i][j][1] - meshU[i - 1][j][1]) / cellWidth

            forwardY = (meshU[i][j + 1][1] - meshU[i][j][1]) / cellHeight
            centralY = (meshU[i][j + 1][1] - meshU[i][j - 1][1]) / (2 * cellHeight)
            backwardY = (meshU[i][j][1] - meshU[i][j - 1][1]) / cellHeight

            Q[i][j][2] = minmod(theta * forwardX, centralX, theta * backwardX)
            Q[i][j][3] = minmod(theta * forwardY, centralY, theta * backwardY)

            # Calculate hv slopes
            forwardX = (meshU[i + 1][j][2] - meshU[i][j][2]) / cellWidth
            centralX = (meshU[i + 1][j][2] - meshU[i - 1][j][2]) / (2 * cellWidth)
            backwardX = (meshU[i][j][2] - meshU[i - 1][j][2]) / cellWidth

            forwardY = (meshU[i][j + 1][2] - meshU[i][j][2]) / cellHeight
            centralY = (meshU[i][j + 1][2] - meshU[i][j - 1][2]) / (2 * cellHeight)
            backwardY = (meshU[i][j][2] - meshU[i][j - 1][2]) / cellHeight

            Q[i][j][4] = minmod(theta * forwardX, centralX, theta * backwardX)
            Q[i][j][5] = minmod(theta * forwardY, centralY, theta * backwardY)

            # if (i == 1 and j == 26):
            #    print "Q= " + str(Q[i][j])

    meshUIntegrationPoints = np.zeros((m, n, 4, 3))  # [i][j][N, S, E, W][w, hu, hv]
    huv = np.zeros((m, n, 4, 3))  # [i][j][N, S, E, W][h, u, v]


    for i in range(m - 1):
        for j in range(n - 1):

            # Calculate U at all four integration points
            meshUIntegrationPoints[i][j][0][0] = meshU[i][j][0] + (cellHeight / 2.0) * Q[i][j][1]  # w - North
            meshUIntegrationPoints[i][j][1][0] = meshU[i][j][0] - (cellHeight / 2.0) * Q[i][j][1]  # w - South
            meshUIntegrationPoints[i][j][2][0] = meshU[i][j][0] + (cellWidth / 2.0) * Q[i][j][0]  # w - East
            meshUIntegrationPoints[i][j][3][0] = meshU[i][j][0] - (cellWidth / 2.0) * Q[i][j][0]  # w - West
            meshUIntegrationPoints[i][j][0][1] = meshU[i][j][1] + (cellHeight / 2.0) * Q[i][j][3]  # hu - North
            meshUIntegrationPoints[i][j][1][1] = meshU[i][j][1] - (cellHeight / 2.0) * Q[i][j][3]  # hu - South
            meshUIntegrationPoints[i][j][2][1] = meshU[i][j][1] + (cellWidth / 2.0) * Q[i][j][2]  # hu - East
            meshUIntegrationPoints[i][j][3][1] = meshU[i][j][1] - (cellWidth / 2.0) * Q[i][j][2]  # hu - West
            meshUIntegrationPoints[i][j][0][2] = meshU[i][j][2] + (cellHeight / 2.0) * Q[i][j][5]  # hv - North
            meshUIntegrationPoints[i][j][1][2] = meshU[i][j][2] - (cellHeight / 2.0) * Q[i][j][5]  # hv - South
            meshUIntegrationPoints[i][j][2][2] = meshU[i][j][2] + (cellWidth / 2.0) * Q[i][j][4]  # hv - East
            meshUIntegrationPoints[i][j][3][2] = meshU[i][j][2] - (cellWidth / 2.0) * Q[i][j][4]  # hv - West

            # Check for water elevations below ground elevation and fix slopes
            # accordingly (see Kurganov 3.20-3.23)
            if (meshUIntegrationPoints[i][j][0][0] < meshBottomIntegrationPoints[i][j + 1][0]):  # N < N?
                # Check to see if the South integration point is also below the ground
                # if (meshUIntegrationPoints[i][j][1][0] < meshBottomIntegrationPoints[i][j][0]):
                #    meshUIntegrationPoints[i][j][0][0] = meshBottomIntegrationPoints[i][j + 1][0]
                #    meshUIntegrationPoints[i][j][1][0] = meshBottomIntegrationPoints[i][j][0]
                # else:
                meshUIntegrationPoints[i][j][0][0] = meshBottomIntegrationPoints[i][j + 1][0]
                meshUIntegrationPoints[i][j][1][0] = 2 * meshU[i][j][0] - meshBottomIntegrationPoints[i][j + 1][0]
            elif (meshUIntegrationPoints[i][j][1][0] < meshBottomIntegrationPoints[i][j][0]):  # S < S?
                meshUIntegrationPoints[i][j][1][0] = meshBottomIntegrationPoints[i][j][0]
                meshUIntegrationPoints[i][j][0][0] = 2 * meshU[i][j][0] - meshBottomIntegrationPoints[i][j][0]
            if (meshUIntegrationPoints[i][j][2][0] < meshBottomIntegrationPoints[i + 1][j][1]):  # E < E?
                # Check to see if the West integration point is also below the ground
                # if (meshUIntegrationPoints[i][j][3][0] < meshBottomIntegrationPoints[i][j][1]):
                #    meshUIntegrationPoints[i][j][2][0] = meshBottomIntegrationPoints[i + 1][j][1]
                #    meshUIntegrationPoints[i][j][3][0] = meshBottomIntegrationPoints[i][j][1]
                # else:
                meshUIntegrationPoints[i][j][2][0] = meshBottomIntegrationPoints[i + 1][j][1]
                meshUIntegrationPoints[i][j][3][0] = 2 * meshU[i][j][0] - meshBottomIntegrationPoints[i + 1][j][1]
            elif (meshUIntegrationPoints[i][j][3][0] < meshBottomIntegrationPoints[i][j][1]):  # W < W?
                meshUIntegrationPoints[i][j][3][0] = meshBottomIntegrationPoints[i][j][1]
                meshUIntegrationPoints[i][j][2][0] = 2 * meshU[i][j][0] - meshBottomIntegrationPoints[i][j][1]

            # Calculate the values of h and (u or v) for each integration point
            # These values are used to calculate the one-sided propagation speeds
            huv[i][j][0][0] = meshUIntegrationPoints[i][j][0][0] - meshBottomIntegrationPoints[i][j + 1][0]  # h - North
            huv[i][j][1][0] = meshUIntegrationPoints[i][j][1][0] - meshBottomIntegrationPoints[i][j][0]  # h - South
            huv[i][j][2][0] = meshUIntegrationPoints[i][j][2][0] - meshBottomIntegrationPoints[i + 1][j][1]  # h - East
            huv[i][j][3][0] = meshUIntegrationPoints[i][j][3][0] - meshBottomIntegrationPoints[i][j][1]  # h - West
            if huv[i][j][0][0] > 0.0:
                huv[i][j][0][1] = (sqrt2 * huv[i][j][0][0] * meshUIntegrationPoints[i][j][0][1]) / (np.sqrt(huv[i][j][0][0] ** 4 + max(huv[i][j][0][0] ** 4, Kappa)))  # u - North
                huv[i][j][1][1] = (sqrt2 * huv[i][j][1][0] * meshUIntegrationPoints[i][j][1][1]) / (np.sqrt(huv[i][j][1][0] ** 4 + max(huv[i][j][1][0] ** 4, Kappa)))  # u - South
                huv[i][j][2][1] = (sqrt2 * huv[i][j][2][0] * meshUIntegrationPoints[i][j][2][1]) / (np.sqrt(huv[i][j][2][0] ** 4 + max(huv[i][j][2][0] ** 4, Kappa)))  # u - East
                huv[i][j][3][1] = (sqrt2 * huv[i][j][3][0] * meshUIntegrationPoints[i][j][3][1]) / (np.sqrt(huv[i][j][3][0] ** 4 + max(huv[i][j][3][0] ** 4, Kappa)))  # u - West
                huv[i][j][0][2] = (sqrt2 * huv[i][j][0][0] * meshUIntegrationPoints[i][j][0][2]) / (np.sqrt(huv[i][j][0][0] ** 4 + max(huv[i][j][0][0] ** 4, Kappa)))  # v - North
                huv[i][j][1][2] = (sqrt2 * huv[i][j][1][0] * meshUIntegrationPoints[i][j][1][2]) / (np.sqrt(huv[i][j][1][0] ** 4 + max(huv[i][j][1][0] ** 4, Kappa)))  # v - South
                huv[i][j][2][2] = (sqrt2 * huv[i][j][2][0] * meshUIntegrationPoints[i][j][2][2]) / (np.sqrt(huv[i][j][2][0] ** 4 + max(huv[i][j][2][0] ** 4, Kappa)))  # v - East
                huv[i][j][3][2] = (sqrt2 * huv[i][j][3][0] * meshUIntegrationPoints[i][j][3][2]) / (np.sqrt(huv[i][j][3][0] ** 4 + max(huv[i][j][3][0] ** 4, Kappa)))  # v - West
            else:
                for k in range(4):
                    huv[i][j][k][1] = 0.0
                    huv[i][j][k][2] = 0.0

            # if (i == 1 and j == 27):
            #    print "huv_east = " + str(huv[i][j][2])
            # if (i == 2 and j == 27):
            #    print "huv_west = " + str(huv[i][j][3])

            # Change hu and hv in U at integration points to reflect the newly calculated velocities
            meshUIntegrationPoints[i][j][0][1] = huv[i][j][0][0] * huv[i][j][0][1]  # hu - North
            meshUIntegrationPoints[i][j][1][1] = huv[i][j][1][0] * huv[i][j][1][1]  # hu - South
            meshUIntegrationPoints[i][j][2][1] = huv[i][j][2][0] * huv[i][j][2][1]  # hu - East
            meshUIntegrationPoints[i][j][3][1] = huv[i][j][3][0] * huv[i][j][3][1]  # hu - West
            meshUIntegrationPoints[i][j][0][2] = huv[i][j][0][0] * huv[i][j][0][2]  # hv - North
            meshUIntegrationPoints[i][j][1][2] = huv[i][j][1][0] * huv[i][j][1][2]  # hv - South
            meshUIntegrationPoints[i][j][2][2] = huv[i][j][2][0] * huv[i][j][2][2]  # hv - East
            meshUIntegrationPoints[i][j][3][2] = huv[i][j][3][0] * huv[i][j][3][2]  # hv - West

    # huv[huv < 0.00000001] = 0.0
    # print4dMatrix(huv, 2, 1)
    # Calculate the one-sided propagation speeds
    # See Kurganov 3.26 - 3.29
    meshOneSidedPropagationSpeeds = np.zeros((m, n, 4))

    for i in range(m - 1):
        for j in range(n - 1):

            try:
                # Calculate North propagation speed of current element
                meshOneSidedPropagationSpeeds[i][j][0] = min(huv[i][j][0][2] - np.sqrt(g * huv[i][j][0][0]),
                                                             huv[i][j + 1][1][2] - np.sqrt(g * huv[i][j + 1][1][0]),
                                                             0.0)

                # Calculate South propagation speed of element above current one
                meshOneSidedPropagationSpeeds[i][j + 1][1] = max(huv[i][j][0][2] + np.sqrt(g * huv[i][j][0][0]),
                                                               huv[i][j + 1][1][2] + np.sqrt(g * huv[i][j + 1][1][0]),
                                                               0.0)

                # Calculate East propagation speed of current element
                meshOneSidedPropagationSpeeds[i][j][2] = min(huv[i][j][2][1] - np.sqrt(g * huv[i][j][2][0]),
                                                             huv[i + 1][j][3][1] - np.sqrt(g * huv[i + 1][j][3][0]),
                                                             0.0)

                # Calculate West propagation speed of the element to the right
                meshOneSidedPropagationSpeeds[i + 1][j][3] = max(huv[i][j][2][1] + np.sqrt(g * huv[i][j][2][0]),
                                                               huv[i + 1][j][3][1] + np.sqrt(g * huv[i + 1][j][3][0]),
                                                               0.0)
            except:
                print "Error: i=" + str(i) + "\tj=" + str(j) + "\nhuv[i][j][0] = " + str(huv[i][j][0]) + "\nhuv[i][j+1][1] = " + str(huv[i][j + 1][1]) + "\thuv[i][j][2] = " + str(huv[i][j][2]) + "\thuv[i+1][j][3] = " + str(huv[i + 1][j][3])
                sys.exit()

    # print "a+: " + str(meshOneSidedPropagationSpeeds[1][27][2]) + "\ta-: " + str(meshOneSidedPropagationSpeeds[2][27][3])

    return meshUIntegrationPoints, meshOneSidedPropagationSpeeds



def minmod(a, b, c):
    if (a > 0.0 and b > 0.0 and c > 0.0):
        return min(a, b, c)
    elif (a < 0.0 and b < 0.0 and c < 0.0):
        return max(a, b, c)
    else:
        return 0.0

def buildRValues(meshH_b, meshFluxes):

    m = meshH_b.shape[0]
    n = meshH_b.shape[1]

    meshR = np.zeros((m - 2, n - 2, 3))

    for i in range(m - 2):
        i += 1
        for j in range(n - 2):
            j += 1
            meshR[i - 1][j - 1] = meshH_b[i][j] - (meshFluxes[i][j][1] - meshFluxes[i - 1][j][1]) - (meshFluxes[i][j][0] - meshFluxes[i][j - 1][0])

    return meshR














































