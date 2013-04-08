'''
Created on Mar 21, 2013

@author: tristan
'''
import numpy as np

g = 9.81

def newFluxSolver(meshUIntegrationPoints, meshBottomIntegrationPoints, meshOneSidedPropagationSpeeds):

    m = meshUIntegrationPoints.shape[0]
    n = meshUIntegrationPoints.shape[1]

    # Each cell calculates and stores only the fluxes across its
    # north [0] and east [1] boundaries
    fluxes = np.zeros((m, n, 2, 3))

    for i in range(m - 2):
        for j in range(n - 2):

            ############### NORTH FLUX STUFF
            b_plus = meshOneSidedPropagationSpeeds[i][j + 1][1]
            b_minus = meshOneSidedPropagationSpeeds[i][j][0]
            U_north = meshUIntegrationPoints[i][j][0]
            U_south = meshUIntegrationPoints[i][j + 1][1]
            B = meshBottomIntegrationPoints[i][j + 1][0]
            G_north = G(U_north, B)
            G_south = G(U_south, B)

            if (b_plus - b_minus != 0.0):
                fluxes[i][j][0] = (b_plus * G_north - b_minus * G_south) / (b_plus - b_minus) + ((b_plus * b_minus) / (b_plus - b_minus)) * (U_south - U_north)
            else:
                fluxes[i][j][0] = np.array([0.0, 0.0, 0.0])

            ############### EAST FLUX STUFF
            a_plus = meshOneSidedPropagationSpeeds[i + 1][j][3]
            a_minus = meshOneSidedPropagationSpeeds[i][j][2]
            U_east = meshUIntegrationPoints[i][j][2]
            U_west = meshUIntegrationPoints[i + 1][j][3]
            B = meshBottomIntegrationPoints[i + 1][j][1]
            F_east = F(U_east, B)
            F_west = F(U_west, B)

            if (a_plus - a_minus != 0.0):
                fluxes[i][j][1] = (a_plus * F_east - a_minus * F_west) / (a_plus - a_minus) + ((a_plus * a_minus) / (a_plus - a_minus)) * (U_west - U_east)
            else:
                fluxes[i][j][1] = np.array([0.0, 0.0, 0.0])

    return fluxes



def fluxSolver(meshUIntegrationPoints, meshBottomIntegrationPoints, meshOneSidedPropagationSpeeds):

    m = meshUIntegrationPoints.shape[0]
    n = meshUIntegrationPoints.shape[1]

    # Each cell calculates and stores only the fluxes across its
    # north [0] and east [1] boundaries
    fluxes = np.zeros((m, n, 2, 3))

    for i in range(1, m - 2):
        for j in range(1, n - 2):

            # Calculate North flux for the current cell

#             b_plus = meshOneSidedPropagationSpeeds[i][j + 1][1]
#             b_minus = meshOneSidedPropagationSpeeds[i][j][0]
#             G_north = G(meshUIntegrationPoints[i][j][0], meshBottomIntegrationPoints[i][j + 1][0])
#             G_south = G(meshUIntegrationPoints[i][j + 1][1], meshBottomIntegrationPoints[i][j + 1][0])
#             top = b_plus * G_north - b_minus * G_south
#             bottom = b_plus - b_minus
#             ans = top / bottom
#             print meshUIntegrationPoints[i][j + 1][1] - meshUIntegrationPoints[i][j][0]

            if (i == 1 and j == 27):
                print "a+: " + str(meshOneSidedPropagationSpeeds[2][27][3]) + "\ta-: " + str(meshOneSidedPropagationSpeeds[1][27][2])
                print "Ue = " + str(meshUIntegrationPoints[i][j][2]) + "Uw = " + str(meshUIntegrationPoints[i + 1][j][3])

            if (meshOneSidedPropagationSpeeds[i][j][0] != 0.0 and meshOneSidedPropagationSpeeds[i][j + 1][1] != 0.0):
                fluxes[i][j][0] = ((meshOneSidedPropagationSpeeds[i][j + 1][1] * G(meshUIntegrationPoints[i][j][0], meshBottomIntegrationPoints[i][j + 1][0]) -
                               meshOneSidedPropagationSpeeds[i][j][0] * G(meshUIntegrationPoints[i][j + 1][1], meshBottomIntegrationPoints[i][j + 1][0])) /
                               (meshOneSidedPropagationSpeeds[i][j + 1][1] - meshOneSidedPropagationSpeeds[i][j][0])) + (
                               (meshOneSidedPropagationSpeeds[i][j + 1][1] * meshOneSidedPropagationSpeeds[i][j][0]) / (meshOneSidedPropagationSpeeds[i][j + 1][1] - meshOneSidedPropagationSpeeds[i][j][0])) * (
                                meshUIntegrationPoints[i][j + 1][1] - meshUIntegrationPoints[i][j][0])
            else:
                fluxes[i][j][0] = 0.0


            # Calculate East flux for the current cell
            if (meshOneSidedPropagationSpeeds[i][j][2] != 0.0 and meshOneSidedPropagationSpeeds[i + 1][j][3] != 0.0):
                fluxes[i][j][1] = ((meshOneSidedPropagationSpeeds[i + 1][j][3] * F(meshUIntegrationPoints[i][j][2], meshBottomIntegrationPoints[i + 1][j][1]) -
                               meshOneSidedPropagationSpeeds[i][j][2] * F(meshUIntegrationPoints[i + 1][j][3], meshBottomIntegrationPoints[i + 1][j][1])) /
                               (meshOneSidedPropagationSpeeds[i + 1][j][3] - meshOneSidedPropagationSpeeds[i][j][2])) + (
                               (meshOneSidedPropagationSpeeds[i + 1][j][3] * meshOneSidedPropagationSpeeds[i][j][2]) / (meshOneSidedPropagationSpeeds[i + 1][j][3] - meshOneSidedPropagationSpeeds[i][j][2])) * (
                                meshUIntegrationPoints[i + 1][j][3] - meshUIntegrationPoints[i][j][2])
            else:
                fluxes[i][j][1] = 0.0

            if (i == 1 and j == 27):
                print "fluxes: " + str(fluxes[i][j][1])

    return fluxes


def F(U, z):

    # U = [w, hu, hv]
    # z = bottom elevation at integration point

    h = U[0] - z
    return np.array([U[1],
                     (U[1] ** 2) / h + 0.5 * g * h ** 2,
                     U[1] * U[2] / h])

def G(U, z):

    # U = [w, hu, hv]
    # z = bottom elevation at integration point

    h = U[0] - z
    return np.array([U[2],
                     U[1] * U[2] / h,
                     (U[2] ** 2) / h + 0.5 * g * h ** 2])



