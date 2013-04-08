'''
Created on Mar 25, 2013

@author: tristan
'''

import numpy as np
from math import sqrt


g = 9.81


def bedSlopeSourceSolver(meshUIntegrationPoints, meshBottomIntegrationPoints, cellWidth, cellHeight):

    m = meshUIntegrationPoints.shape[0]
    n = meshUIntegrationPoints.shape[1]

    # Bed slope source terms (cell-centered): [i][j][[0, xDir, yDir]]
    meshH_b = np.zeros((m, n, 3))

    for i in range(m - 1):
        for j in range(n - 1):
            hX = ((meshUIntegrationPoints[i][j][2][0] - meshBottomIntegrationPoints[i + 1][j][1]) + (meshUIntegrationPoints[i][j][3][0] - meshBottomIntegrationPoints[i][j][1])) / 2.0
            hY = ((meshUIntegrationPoints[i][j][0][0] - meshBottomIntegrationPoints[i][j + 1][0]) + (meshUIntegrationPoints[i][j][1][0] - meshBottomIntegrationPoints[i][j][0])) / 2.0
            slopeX = (meshBottomIntegrationPoints[i + 1][j][1] - meshBottomIntegrationPoints[i][j][1]) / cellWidth
            slopeY = (meshBottomIntegrationPoints[i][j + 1][0] - meshBottomIntegrationPoints[i][j][0]) / cellHeight

            meshH_b[i][j][1] = -g * slopeX * hX
            meshH_b[i][j][2] = -g * slopeY * hY

    return meshH_b


# def bedShearSourceSolver(meshU, meshBottomIntegrationPoints):
#
#     m = meshU.shape[0]
#     n = meshU.shape[1]
#
#     meshH_f = np.zeros((m, n, 3))
#
#     for i in range(m):
#         for j in range(n):
#             bottomElevation = (meshBottomIntegrationPoints[i + 1][j][1] - meshBottomIntegrationPoints[i][j][1]) / 2.0
#             meshH_f[i][j][1] = -g * sqrt((meshU[i][j][1] ** 2 + meshU[i][j][2]))



def bedShearSourceSolver(meshU, meshBottomCenters, manningsN):

    m = meshU.shape[0]
    n = meshU.shape[1]

    meshShearSource = np.zeros((m, n, 3))

    for i in range(m):
        for j in range(n):
            h = meshU[i][j][0] - meshBottomCenters[i][j]
            sol = (-g * sqrt(meshU[i][j][1] ** 2 + meshU[i][j][2] ** 2)) / (h ** 2 * (h ** (1.0 / 6.0) / manningsN) ** 2)
            meshShearSource[i][j][1] = sol
            meshShearSource[i][j][2] = sol

    return meshShearSource


def newBedShearSourceSolver(meshU, meshBottomCenters, manningsN, cellWidth, cellHeight):

    m = meshU.shape[0]
    n = meshU.shape[1]

    sqrt2 = np.sqrt(2.0)
    K0 = 0.01
    Kappa = K0 * max(1.0, min(cellWidth, cellHeight))

    meshShearSource = np.zeros((m, n, 3))

    for i in range(m):
        for j in range(n):
            if (meshU[i][j][0] > meshBottomCenters[i][j]):
                h = meshU[i][j][0] - meshBottomCenters[i][j]
                u = (sqrt2 * h * meshU[i][j][1]) / (np.sqrt(h ** 4 + max(h ** 4, Kappa)))
                v = (sqrt2 * h * meshU[i][j][2]) / (np.sqrt(h ** 4 + max(h ** 4, Kappa)))
                Cz = (h ** (1.0 / 6.0)) / manningsN
                solution = (-g * np.sqrt(u ** 2 + v ** 2)) / (h * Cz ** 2)
                meshShearSource[i][j][1] = solution
                meshShearSource[i][j][2] = solution
            else:
                meshShearSource[i][j] = np.zeros(3)

    return meshShearSource



