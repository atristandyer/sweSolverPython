'''
Created on Mar 19, 2013

@author: tristan
'''

def calculateTimestep(meshOneSidedPropagationSpeeds, dx, dy):
    amax = -99999.0
    bmax = -99999.0

    m = meshOneSidedPropagationSpeeds.shape[0]
    n = meshOneSidedPropagationSpeeds.shape[1]

    for i in range(2, m - 2):
        for j in range(2, n - 2):
            if abs(meshOneSidedPropagationSpeeds[i][j][0]) > bmax:
                bmax = abs(meshOneSidedPropagationSpeeds[i][j][0])
            if abs(meshOneSidedPropagationSpeeds[i][j][1]) > bmax:
                bmax = abs(meshOneSidedPropagationSpeeds[i][j][1])
            if abs(meshOneSidedPropagationSpeeds[i][j][2]) > amax:
                amax = abs(meshOneSidedPropagationSpeeds[i][j][2])
            if abs(meshOneSidedPropagationSpeeds[i][j][3]) > amax:
                amax = abs(meshOneSidedPropagationSpeeds[i][j][3])

    # print amax, bmax
    dt = min(dx / (4.0 * amax), dy / (4.0 * bmax))

    return dt


def buildU_star(meshUstar, meshU, dt, meshR, meshShearSource):

    m = meshR.shape[0]
    n = meshR.shape[1]

    for i in range(m):
        for j in range(n):
            meshUstar[i + 1][j + 1] = (meshU[i + 1][j + 1] + dt * meshR[i][j]) / (1 + dt * meshShearSource[i + 1][j + 1])

    return meshUstar



def buildU_next(meshUnext, meshU, meshUstar, dt, meshRstar, meshShearSourcestar):

    m = meshRstar.shape[0]
    n = meshRstar.shape[1]

    for i in range(m):
        for j in range(n):
            meshUnext[i + 1][j + 1] = (0.5 * meshU[i + 1][j + 1] + 0.5 * (meshUstar[i + 1][j + 1] + dt * meshRstar[i][j])) / (1 + 0.5 * dt * meshShearSourcestar[i + 1][j + 1])

    return meshUnext

























