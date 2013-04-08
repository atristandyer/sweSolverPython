'''
Created on Mar 21, 2013

@author: tristan
'''


def printMatrix(matrix, padding=2):
    dim = matrix.ndim
    m = matrix.shape[0]
    n = matrix.shape[1]

    if dim == 3:
        o = matrix.shape[2]
        for i in reversed(range(m)):
            line = ""
            for j in range(n):
                if (i >= padding and i <= m - padding - 1 and j >= padding and j <= n - padding - 1):
                        line += "{"
                else:
                    line += "["
                for k in range(o):
                    line += "%.2f" % matrix[j][i][k]
                    if k == o - 1:
                        if (i >= padding and i <= m - padding - 1 and j >= padding and j <= n - padding - 1):
                            line += "}"
                        else:
                            line += "]"
                    else:
                        line += ", "
                line += "\t"
            print line

def print4dMatrix(matrix, padding=2, dim=0):

    m = matrix.shape[0]
    n = matrix.shape[1]
    o = matrix.shape[3]
    for i in reversed(range(m)):
        line = ""
        for j in range(n):
            if (i >= padding and i <= m - padding - 1 and j >= padding and j <= n - padding - 1):
                    line += "{"
            else:
                line += "["
            for k in range(o):
                line += "%.2f" % matrix[j][i][dim][k]
                if k == o - 1:
                    if (i >= padding and i <= m - padding - 1 and j >= padding and j <= n - padding - 1):
                        line += "}"
                    else:
                        line += "]"
                else:
                    line += ", "
            line += "\t"
        print line


