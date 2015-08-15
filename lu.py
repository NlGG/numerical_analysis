# coding: utf-8

import numpy as np

def lu(matrix):
    A = matrix
    N = len(A)
    u_t = np.matrix(np.zeros((N, N), dtype=float))
    l = np.matrix(np.zeros((N, N), dtype=float))
    l_t = np.transpose(l)

    for i in range(0, N):
        A_t = np.transpose(A)
        u_t[i] = A[i]
        l_t[i] = A_t[i]/A[i, i]
        l = np.transpose(l_t)
        A = A - l[:, i]*u_t[i]
    
    print "L =", l
    print "U =", u_t
