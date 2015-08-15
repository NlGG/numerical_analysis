# coding: utf-8

import numpy as np
from numpy import linalg as la
import math 

class QR():
    def __init__(self, matrix):
        self.matrix  = matrix
        self.N = len(self.matrix)
        
    def hessenberg_matrix(self):
        A = self.matrix
        N = len(A)
        I = np.identity((N), dtype = float)  
        
        for i in range(N-2):
            A_t = np.transpose(A)
            sum_x_j = 0
            x_t = A_t[i]          
            for k in range(i+1, N):
                sum_x_j += pow(x_t[0,k], 2)                
            y_t = np.matrix(np.zeros(N))
            for n in range(i+1):
                y_t[0, n] = x_t[0, n]
            y_t[0, i+1] = (-np.sign(x_t[0, i+1]))*math.sqrt(sum_x_j)
            u_t = x_t - y_t
            u = np.transpose(u_t)
            H = I - 2*u*u_t/(pow(np.linalg.norm(u_t[0]),2))                             
            A = H*A
            A = A*H
            
        return A
        

    def qr(self, A=None):
        if A is None:
            A = self.hessenberg_matrix()
        N = len(A)
        I = np.identity((N), dtype = float)
        
        for i in range(N-1):
            A_t = np.transpose(A)
            sum_x_j = 0
            x_t = A_t[i]
            for k in range(i, N):
                sum_x_j += pow(x_t[0,k], 2)
            y_t = np.matrix(np.zeros(N))              
            for n in range(i):
                y_t[0, n] = x_t[0, n]
            y_t[0, i] = (-np.sign(x_t[0, i]))*math.sqrt(sum_x_j)
            u_t = x_t - y_t
            u = np.transpose(u_t)
            H = I - 2*u*u_t/(pow(np.linalg.norm(u_t[0]),2))                    
            A = H*A
            if i == 0:
                H_p = H
            else: 
                H_p = H * H_p
            qr =[[],[]]
            qr[0] = H_p
            qr[1] = A               
        return qr
                
    def qr_method(self, A=1, t=10):
        if A == 0:
            A = self.matrix
        elif A == 1:
            A = self.hessenberg_matrix()
        for i in range(t):
            qr = self.qr(A)
            Q = qr[0]
            R = qr[1]
            A = R*Q         
        return A
