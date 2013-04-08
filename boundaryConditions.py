'''
Created on Mar 26, 2013

@author: tristan
'''

# Grid padding must be 2 for this boundary condition
def applyWallBoundaries(meshU, gridSize, gridPadding):

    m = meshU.shape[0]
    n = meshU.shape[1]

    # Left and right side ghost cells
    for j in range(n):
        meshU[0][j] = meshU[3][j]
        meshU[0][j][1] = -meshU[0][j][1]
        meshU[1][j] = meshU[2][j]
        meshU[1][j][1] = -meshU[1][j][1]
        meshU[m - 1][j] = meshU[m - 4][j]
        meshU[m - 1][j][1] = -meshU[m - 1][j][1]
        meshU[m - 2][j] = meshU[m - 3][j]
        meshU[m - 2][j][1] = -meshU[m - 2][j][1]

    # Bottom and top ghost cells
    for i in range(m):
        meshU[i][0] = meshU[i][3]
        meshU[i][0][2] = -meshU[i][0][2]
        meshU[i][1] = meshU[i][2]
        meshU[i][1][2] = -meshU[i][1][2]
        meshU[i][n - 1] = meshU[i][n - 4]
        meshU[i][n - 1][2] = -meshU[i][n - 1][2]
        meshU[i][n - 2] = meshU[i][n - 3]
        meshU[i][n - 3][2] = -meshU[i][n - 3][2]


    return meshU

def applyFreeBoundaries(meshU, gridSize, gridPadding):

    m = meshU.shape[0]
    n = meshU.shape[1]

    # Left side ghost cells
    for i in range(gridPadding):
        for j in range(gridPadding, n - gridPadding):
            meshU[i][j] = meshU[gridPadding][j]

    # Bottom ghost cells
    for j in range(gridPadding):
        for i in range(gridPadding, m - gridPadding):
            meshU[i][j] = meshU[i][gridPadding]

    # Top ghost cells
    for j in range(n - gridPadding, n):
        for i in range(gridPadding, m - gridPadding):
            meshU[i][j] = meshU[i][n - gridPadding - 1]

    # Right ghost cells
    for i in range(m - gridPadding, m):
        for j in range(gridPadding, n - gridPadding):
            meshU[i][j] = meshU[m - gridPadding - 1][j]

    # Bottom left corner ghost cells
    for i in range(gridPadding):
        for j in range(gridPadding):
            meshU[i][j] = meshU[gridPadding][gridPadding]

    # Top left corner ghost cells
    for i in range(gridPadding):
        for j in range(n - gridPadding, n):
            meshU[i][j] = meshU[gridPadding][n - gridPadding - 1]

    # Bottom right corner ghost cells
    for i in range(m - gridPadding, m):
        for j in range(gridPadding):
            meshU[i][j] = meshU[m - gridPadding - 1][gridPadding]

    # Top right corner ghost cell
    for i in range(m - gridPadding, m):
        for j in range(n - gridPadding, n):
            meshU[i][j] = meshU[m - gridPadding - 1][n - gridPadding - 1]

    return meshU
