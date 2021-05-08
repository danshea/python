import numpy as np
def levenshtein(a,b):
    if len(a) == 0:
        return len(b)
    if len(b) == 0:
        return len(a)
    I,J = len(a)+1,len(b)+1
    M = np.zeros((I,J))
    for i in range(I):
        M[i,0] = i
    for j in range(J):
        M[0,j] = j
    # loop over b
    for j, lb in enumerate(b,1):
        for i, la in enumerate(a,1):
            if la == lb:
                cost = 0
            else:
                cost = 1
            M[i,j] = min([M[i,j-1]+1, M[i-1,j]+1, M[i-1,j-1]+cost])
    return {'matrix': M, 'score': M[-1,-1]}