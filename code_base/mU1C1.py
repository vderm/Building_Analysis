#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 12:15:50 2015

@author: Vasken

Variable Names
Uin: Conductance matrix input by user, upper triangle only, (nN x nN) (W/K)
U: Conductance matrix (symmetrical) with added capacitance for diagonal term, (nN x nN) (W/K)
C: Capacitance vector, (nN x 1) (J/K)
F: Conductance matrix of nodes connected to a known temperature source, (nN x nM) (W/K)
T: Temperature vector per timestep, (nT x nN) (degC)
TK: Temperature vector of known temperatures per timestep, (nT x nM) (degC)
Qin: Heat flow, only external sources, (nN x 1) (W)
Q: Heat flow vector + external sources + capacitance from previous timestep (implicit only), (nN x 1) (W)

nN: Number of nodes
nM: Number of nodes with known temperatures
nT: Number of timesteps

Node Number: Object
0: effective room node, connected to capacitor and T_ambient (in the Ms)

Node Number with known temperatures: Object
0: ambient air
"""

def mU1C1(U_in, C_in, dt):
    # Load dependencies
    from numpy import zeros
    from numpy import sum as npsum
    from numpy.linalg import inv

    # #### Control
    nN = 1          # number of nodes
    nM = 1          # number of nodes with known temperatures

    #%% Nodal Connections
    # Declare variables
    Uin = zeros((nN,nN))     # K/W
    F = zeros((nN,nM))       # K/W
    C = zeros((nN,1))        # J/K

    # How are the nodes connected?
    # Uin[0,1] = (1/R + U + dx/kA)**-1

    # Connected to temperature sources
    F[0,0] = (1/U_in)**-1

    # Nodes with capacitance
    C[0] = C_in

    #%% U-matrix completion, and its inverse
    U = -Uin - Uin.T  # U is symmetrical, non-diagonals are -ve
    s = -npsum(U,1)
    for i in range(0,nN):
        U[i,i] = s[i] + npsum(F[i,]) + C[i]/dt
    U_inv = inv(U)

    #%% Ship it
    return (U_inv, F, C, nN, nM)
