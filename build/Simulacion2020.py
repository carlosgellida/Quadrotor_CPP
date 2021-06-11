#####################################################################
# Program to calculate dynamic of a quadrotor with 3 DoF robot amr  #
#####################################################################

import numpy as np 
import math 
from scipy.linalg import expm
from math import sin, cos

def inv_xi (Chi): 
    # Process to efficiently invert an EHTM 
    R = Chi[:3,:3]
    PxR = Chi[:3,3:]
    R_T = np.transpose(R)
    PxR_T = -np.matmul(np.matmul(R_T, PxR),R_T) 
    Chi_inv_U = np.hstack((R_T, PxR_T))
    Chi_inv_L = np.hstack((np.zeros((3,3)), R_T))
    return np.vstack((Chi_inv_U, Chi_inv_L))
    
def skew3 (a):
    # Skew symmetric for R3 vectors
    A = np.array([[0, -a[2], a[1]], 
                  [a[2], 0, -a[0]],
                  [-a[1], a[0], 0]])
    return A 

def skew31_ext (S):
    # Skew symmetric for R3 vectors
    A = np.zeros((S.shape[0], 2, 3, 3))
    A[:, 0, 0, 1] = -S[:, 2]
    A[:, 0, 1, 0] = S[:, 2]
    A[:, 0, 0, 2] = S[:, 1]
    A[:, 0, 2, 0] = -S[:, 1]
    A[:, 0, 2, 1] = -S[:, 0]
    A[:, 0, 1, 2] = -S[:, 0]
    
    A[:, 1, 0, 1] = -S[:, 5]
    A[:, 1, 1, 0] = S[:, 5]
    A[:, 1, 0, 2] = S[:, 4]
    A[:, 1, 2, 0] = -S[:, 4]
    A[:, 1, 1, 2] = -S[:, 3]
    A[:, 1, 2, 1] = S[:, 3]
    return A

def skew32_ext (P):
    # Skew symmetric for R3 vectors
    A = np.zeros((P.shape[0], 3, 3))
    A[:, 0, 1] = -P[:, 2]
    A[:, 1, 0] = P[:, 2]
    A[:, 0, 2] = P[:, 1]
    A[:, 2, 0] = -P[:, 1]
    A[:, 2, 1] = -P[:, 0]
    A[:, 1, 2] = P[:, 0]
    return A

def skew61 (b): 
    # First extended skew symmetric
    B = np.vstack((np.hstack((skew3(b[3:6]), skew3(b[:3]))),
                   np.hstack((np.zeros((3,3)), skew3(b[3:6]))) ))
    return B 

def skew61_ext (S): 
    # First extended skew symmetric
    A = skew31_ext(S)
    B = np.zeros((S.shape[0], 6, 6))
    B[:, :3, :3] = A[:, 1, :, :]
    B[:, :3, 3:] = A[:, 0, :, :]
    B[:, 3:, 3:] = A[:, 1, :, :]
    
    return B 

def skew62 (b): 
    # Second extended skew symmetric
    B = np.vstack((np.hstack(( np.zeros((3,3)), skew3(b[:3]) )) ,
                   np.hstack(( skew3(b[:3]), skew3(b[3:6])    )) ))
    return B  

def skew62_ext (V):
    A = skew31_ext(V)
    B = np.zeros((V.shape[0], 6, 6))
    B[:, :3, 3:] = A[:, 0, :, :]
    B[:, 3:, :3] = A[:, 0, :, :]
    B[:, 3:, 3:] = A[:, 1, :, :]
    return B 

def Exponential (S,q):
    if sum(S[:3]) == 0: 
        Exp = np.identity((6)) + math.sin(q) * skew61(S) + (1 - math.cos(q)) * np.matmul(skew61(S),skew61(S))
    else:
        Exp = np.identity((6)) + q*skew61(S)
    return Exp

def Chi_zero_short (P):
    Chi_0 = np.vstack((np.hstack((np.identity(3), skew3(P) )),
                           np.hstack((np.zeros((3,3)), np.identity(3) )) ))
    return Chi_0

def Exponential_ext (S, q, Nj):
    Iext = np.zeros((Nj, 6, 6))
    Iext[:, 0, 0] = 1
    Iext[:, 1, 1] = 1
    Iext[:, 2, 2] = 1
    Iext[:, 3, 3] = 1
    Iext[:, 4, 4] = 1
    Iext[:, 5, 5] = 1
    K = skew61_ext(S)
    Exp = Iext + np.reshape(np.sin(q), (Nj, 1, 1)) * K + np.reshape((1 - np.cos(q)), (Nj, 1, 1)) * np.matmul(K, K)
    return Exp 

def Chi_zero_ext (P):
    Chi_0 = np.zeros((P.shape[0], 6, 6))
    A = skew32_ext(P)
    Chi_0[:, :3, :3] = np.identity(3)
    Chi_0[:, :3, 3:] = A
    Chi_0[:, 3:, 3:] = np.identity(3)
    return Chi_0

def matmul2(b, a):
    ''' Make an special product defined by:
    np.matmul(b[i, :], a[i, :, :]) but avoiding for cycle '''
    lenght = b.shape[0]*b.shape[1] 
    width = int((a.shape[0]*a.shape[1]*a.shape[2])/lenght)
    b2 = b.reshape(1, lenght)
    a2 = np.transpose(a, axes=((2,0,1))).reshape(width, lenght)
    result1 = np.reshape(np.transpose(a2*b2, axes=((1,0))), (a.shape[0], a.shape[1], a.shape[2]) )
    result = np.sum(result1, axis=1)
    return result

def EIM (I_info1, I_info2):
    
    I = np.zeros((np.shape(I_info1)[0], 3, 3))
    R = np.zeros((np.shape(I_info1)[0], 3, 3))
    Id = np.zeros((np.shape(I_info1)[0], 3, 3))
    Id = np.identity(3)
    Px = np.zeros((np.shape(I_info1)[0], 3, 3))
    #M = np.zeros((np.shape(I_info1)[0], 6, 6))
    M = np.zeros((int(I_info2[0, -1]) + 1, 6, 6))
    # Agregar una función que sume las EIM de los cuerpos que están 
    # unidos rígidamente, usando el máximo índice de I_info2[0, i]

    I[:, 0, 0] = I_info1[:, 0, 0]
    I[:, 1, 1] = I_info1[:, 0, 1]
    I[:, 2, 2] = I_info1[:, 0, 2]
    I[:, 0, 1] = I_info1[:, 1, 0]
    I[:, 1, 0] = I_info1[:, 1, 0]
    I[:, 0, 2] = I_info1[:, 1, 1]
    I[:, 2, 0] = I_info1[:, 1, 1]
    I[:, 1, 2] = I_info1[:, 1, 2]
    I[:, 2, 1] = I_info1[:, 1, 2]
    
    theta = np.linalg.norm(I_info1[:, 3, :], axis=1)
    K = np.zeros((np.shape(I_info1)[0], 3, 3))
    j = 0  
    
    for i in range(np.shape(I_info1)[0]):
        Px[i, :, :] = skew3(I_info1[i, 2, :])
        if np.allclose(I_info1[i, 3, :], 0): 
            R[i, :, :] = np.identity(3)
        else:
            K[i, :, :] = skew3(I_info1[i, 3, :]/ theta[i])
            R[i, :, :] = np.identity(3) + math.sin(theta[i])*K[i, :, :] + (1 - math.cos(theta[i]))*np.matmul(K[i, :, :], K[i, :, :])
        
        
        if I_info2[0, i] == I_info2[0, i-1]:
            j = j - 1 
            M[j, 3:, 3:] += -I_info2[1, i]*np.matmul(Px[i, :, :], Px[i, :, :]) + np.matmul(R[i, :, :], np.matmul(I[i, :, :], np.transpose(R[i, :, :])))
            M[j, :3, :3] += I_info2[1, i]*np.identity(3)
            M[j, 3:, :3] += -I_info2[1, i]*Px[i, :, :]
            M[j, :3, 3:] += I_info2[1, i]*Px[i, :, :]  
        else:
            M[j, 3:, 3:] = -I_info2[1, i]*np.matmul(Px[i, :, :], Px[i, :, :]) + np.matmul(R[i, :, :], np.matmul(I[i, :, :], np.transpose(R[i, :, :])))
            M[j, :3, :3] = I_info2[1, i]*np.identity(3)
            M[j, 3:, :3] = -I_info2[1, i]*Px[i, :, :]
            M[j, :3, 3:] = I_info2[1, i]*Px[i, :, :]        
            
        
        j = j + 1
        
    return M, I


def Inverse_dyn_imp(Chi_base, dot_q, q, P, S_0, Nj, Vp, M):
    
    J = np.zeros((Nj+1, 6, Nj+6)) # Screw axis expressed at descendants coordinate systems
    dotJ = np.zeros((Nj+1, 6, Nj+6))
    V = np.zeros((Nj+1, 6)) # Twis of each body w.r.t. inertial frame in local coordinates
    V_frame0 = np.zeros((Nj+1, 6)) # Twist of each body in inertial frame coordinates
    Hi = np.zeros((Nj+1, Nj+6, Nj+6))
    Mu = np.zeros((Nj+1, 6, 6)) #Body momentum 
    Ci1 = np.zeros((Nj+6, Nj+6, Nj+6))
    Ci2 = np.zeros((Nj+6, Nj+6, Nj+6))
    Hi = np.zeros((Nj+6, Nj+6, Nj+6))
    
    R_storage = np.zeros((Nj+1, 6, 6))
    R_storage[0, 3:, 3:] = np.copy(Chi_base[3:, 3:])
    R_storage[0, :3, :3] = np.copy(Chi_base[3:, 3:])

    J[0, :6, :6] = np.identity(6)
    V[0, :] = J[0, :, :] @ dot_q

    V_frame0[0, :] = R_storage[0, :, :] @V[0, :]
    
    for i in range(1, Nj+1):
        S = S_0[i-1]
        Exp = Exponential(-S, q[i+5])
        Chi = Exp @ Chi_zero_short(-P[i-1]) 
        R = np.copy(Chi).T
        R[3:, :3] = np.zeros((3, 3))
        #Chi_inverse = inv_xi(Chi)

        J[i, :, i+5] = S
        J[i, :, :] = J[i, :, :] + Chi @ J[Vp[i], :, :] 
        dotJ[i, :, :] =  (Chi @ dotJ[Vp[i], :, :]) - (skew61(S*dot_q[i+5]) @ Chi @ J[Vp[i], :, :] )
        R_storage[i, :, :] = R[:, :] @ R_storage[Vp[i], :, :] 
        V[i, :] = J[i, :, :]@dot_q
        V_frame0[i, :] = R_storage[i, :, :] @ V[i, :]
        Mu[i, :, :] = skew62( M[i, :, :] @ V[i, :] )
        Hi[i, :, :] = (J[i, :, :].T @ M[i, :, :]) @ J[i, :, :]
        Ci1[i, :, :] = (J[i, :, :].T @ M[i, :, :]) @ dotJ[i, :, :]
        Ci2[i, :, :] =  (J[i, :, :].T @ Mu[i, :, :]) @ J[i, :, :]

    Mu[0, :, :] = skew62( M[0, :, :] @ dot_q[:6] )
    
    Hi[0, :, :] = (J[0, :, :].T @ M[0, :, :]) @ J[0, :, :]
    Ci1[0, :, :] = (J[0, :, :].T @ M[0, :, :]) @ dotJ[0, :, :]
    Ci2[0, :, :] =  (J[0, :, :].T @ Mu[0, :, :]) @ J[0, :, :]

    H = np.sum(Hi, axis=0)
    
    Ci = Ci1 - Ci2
    C = np.sum(Ci, axis=0)
    
       
    return H, C, V_frame0, Chi 

def gs_cofficient(v1, v2):
    return np.dot(v2, v1) / np.dot(v1, v1)

def multiply(cofficient, v):
    return map((lambda x : x * cofficient), v)




def direct_dynamic(dot_q, Chi_base, q):
    H, C, V_frame0, Chi = Inverse_dyn_imp(Chi_base, dot_q, q, P, S_0, Nj, Vp, M)
    b = (tau - C@dot_q)
    #b = tau
    L = np.linalg.cholesky(H)
    x = np.linalg.solve(L, b)
    ddot_q = np.linalg.solve(L.transpose(), x)
    return ddot_q, V_frame0 


'''Introduction of robot configuration:'''


'''Rows of matrix S contains the axis of action of
    each joint'''
    
S_0 = np.zeros((7, 6))
S_0[0, 5] = 1
S_0[1, 5] = 1
S_0[2, 5] = 1
S_0[3, 5] = 1
S_0[4, 5] = 1
S_0[5, 3] = 1
S_0[6, 3] = 1

Nb = 8  # Number of bodies including base one and auxiliar ones
Nj = Nb-1  # Number of joints including auxiliars

''' Rows of matrix P contain the position
        vectors of frame i w.r.t. frame i-1'''
        
P = np.array([[0.25, -0.25, 0.1],
              [0.25, 0.25, 0.1],
              [-0.25, 0.25, 0.1],
              [-0.25, -0.25, 0.1],
              [0, 0, -0.1],
              [0, 0, -0.05],
              [0, 0.2*sin(0.7854), -0.2*cos(0.7854)]])

Vp = np.array([-1, 0, 0, 0, 0, 0, 5, 6]) # Parent's vector 


''' Introduction of inertia and center of mass '''

I_info1 = np.zeros((10, 4, 3))

I_info1[0, 0, :] = np.array([0.0008125, 0.0008125, 0.001125])

I_info1[1, 0, :] = np.array([0.66667e-6, 0.00444417, 0.00444417])
I_info1[1, 3, :] = np.array([0, 0, 0.7854])

I_info1[2, 0, :] = np.array([0.66667e-6, 0.00444417, 0.00444417])
I_info1[2, 3, :] = np.array([0, 0, -0.7854])

I_info1[3, 0, :] = np.array([1.36561e-05, 1.36771e-05, 2.13783e-08])
I_info1[3, 3, :] = np.array([-1.5708, 0, 0]) #-90 degrees rotation over x axis

I_info1[4, 0, :] = np.array([1.36561e-05, 1.36771e-05, 2.13783e-08])
I_info1[4, 3, :] = np.array([-1.5708, 0, 0]) #-90 degrees rotation over x axis

I_info1[5, 0, :] = np.array([1.36561e-05, 1.36771e-05, 2.13783e-08])
I_info1[5, 3, :] = np.array([-1.5708, 0, 0]) #-90 degrees rotation over x axis

I_info1[6, 0, :] = np.array([1.36561e-05, 1.36771e-05, 2.13783e-08])
I_info1[6, 3, :] = np.array([-1.5708, 0, 0]) #-90 degrees rotation over x axis

I_info1[7, 0, :] = np.array([8.41667e-05, 8.41667e-05, 1.66667e-06])
#I_info1[7, 2, :] = np.array([0, 0, 0.05])    #traslation

I_info1[8, 0, :] = np.array([6.68333e-04, 6.68333e-04, 3.33333e-06])
I_info1[8, 2, :] = np.array([0, -0.1*sin(0.7854), 0.1*cos(0.7854)])    #traslation
I_info1[8, 3, :] = np.array([0.7854, 0, 0]) #Rotation 

I_info1[9, 0, :] = np.array([6.68333e-04, 6.68333e-04, 3.33333e-06])
I_info1[9, 2, :] = np.array([0, -0.1, 0])
I_info1[9, 3, :] = np.array([1.5708, 0, 0])

I_info2 = np.array([[0, 0, 0, 1, 2, 3, 4, 5, 6, 7],
                   [0.3, 0.1, 0.1, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 0.1, 0.2, 0.2]])

''' Initial conditions are declared here: '''

dot_q_0 = np.zeros(13)
dot_q_0[10] = 0.1
#dot_q_0[6] = 300
#dot_q_0[7] = 300
#dot_q_0[8] = 300
#dot_q_0[9] = 300
#dot_q_0[12] = 0.05

q_part_0 = np.zeros(7)  # It doesn't include the attitude of base body

R_base_0 = np.identity(3) 
P_base_0 = np.ones((3, 1))

Chi_base_0 = np.zeros((6, 6))

Chi_base_0[:3, :3] = R_base_0
Chi_base_0[3:6, 3:6] = R_base_0
Chi_base_0[:3, 3:6] = skew3(P_base_0)

''' Beginning of algorithm '''

dot_q = dot_q_0
Chi_base = Chi_base_0
#q_part = q_part_0
q_part = np.zeros((Nj+6))
q_part[6:] = q_part_0
#q_part2[Nj] = Chi_base[2, 4]
#q_part2[Nj+1] = Chi_base[1, 3]
#q_part2[Nj+2] = Chi_base[0, 5]


M, I = EIM(I_info1, I_info2)

tau = np.zeros((dot_q_0.shape))
#tau[3] = -0.01

''' Integration of Runge Kutta, Bogacki–Shampine Method '''

dot_q_save = np.zeros(dot_q.shape)
q = np.zeros((Nj+6))
q_save = np.zeros((Nj+6))

t0 = 0
tend = 10
h = 5e-3
n = int((tend-t0)/h)


V_save = np.zeros((Nj+1, 6))
V_0_save = np.zeros((6))
ddot_q2, V_frame0 = direct_dynamic(dot_q, Chi_base, q)
V_save = V_frame0
R_base = np.zeros((6, 6))


for i in range(0,n):
    t = i*h

    k1, V_frame0 = direct_dynamic(dot_q, Chi_base, q_part)
    k2, Vfake = direct_dynamic(dot_q + 0.5*k1*h, Chi_base, q_part)
    k3, Vfake = direct_dynamic(dot_q + 0.75*k2*h, Chi_base, q_part)   
    
    
    
    dot_q = dot_q + (2/9)*k1*h + (1/3)*k2*h + (4/9)*k3*h
    
    k1b = dot_q
    k2b = dot_q + 0.5*k1b*h
    k3b = dot_q + 0.75*k2b*h
    
    q_part_increment = (2/9)*k1b*h + (1/3)*k2b*h + (4/9)*k3b*h
    q_part = q_part + q_part_increment 
        
    q_save = np.vstack((q_save, q_part))
    dot_q_save = np.vstack((dot_q_save, dot_q))
    V_save = np.vstack((V_save, V_frame0))

    Chi_base = Chi_base@expm(skew61(q_part_increment[:6]))
    print(t*100/tend, "%") 

    
 
np.savetxt("q.txt", q_save)
np.savetxt("dotq.txt", dot_q_save)
np.savetxt("V.txt", V_save)
