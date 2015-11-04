'''
============================================
Atmos_tools
--------------------------------------------
A collection of functions in Atmospheric Sciences
                ----- Author: Yingkai Sha
============================================
2014/12/21 File created
'''

import numpy as np

def central_diff(T)
    '''
    % ==========================================================================
    % Calculate the derivative of variable on Earth by central difference method
    %   build-in function
    % ==========================================================================
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/2/25
    % ==========================================================================
        Moved from MATLAB to Python on 2014/12/21
    ----------------------------------------------------------------------------
        dx, dy = central_diff(T)
    ============================================================================
    '''
    M, N=np.size(T);
    dx=np.zeros([M, N]);
    dy=np.zeros([M, N]);
    # dx
    dx[:, 0]=T[:, 1]-T[:, 0];
    dx[:, N-1]=T[:, N-1]-T(:, N-2];    
    for j in range(1, N-2):
        dx[:, j]=(T[:, j+1]-T[:, j-1])/2;
    # dy
    dy[0, :]=T[1, :]-T[0, :];
    dy[M-1, :]=T[M, :]-T[M-1, :];
    for i in range(1, M-2):
        dy[i, :]=(T[i+1, :]-T[i-1, :])/2;
    return dx, dy
    
def EOF(H, nmode=10, ndim=3, reverse=1):
    '''
    Converted from MATLAB to Python 2.7 code @ 2015/06/15 - YKS
     + ndim: [LAT, LON, TIME] data (=3) or [MAP, TIME] data (=2)
     + reverse: normalized spatial pattern (=0), normalized PC (=1)
    % ======================================================================= %
    % Input
    %   H: Variable required for EOF comutation, H(LAT, LON, Time) 
    %       or H(Space, Time) is accepted.
    %   nmode: Number of modes output
    % Output
    %   EOFs: EOF Spatial Pattern
    %   PC: Timeseries cooresponding to each EOFs
    %   expvar: Explained variance
    % ======================================================================= %
    % Author
    %   Yingkai Sha
    %       yingkaisha@gmail.com
    % 2014/3/18
    % ======================================================================= %
    '''
    ##import scipy.linalg.eig as eig
    # Get the size of array
    if ndim == 3:
        LAT, LON, T = H.shape
    elif ndim == 2:
        LON, T = H.shape
        LAT = 1
    # Covarience
    H = np.reshape(H, [LAT*LON, T]).T
    R=np.dot(H, H.T); N = np.size(R, 0)
    # Allocation
    PC     = np.zeros([nmode, T]);
    expvar = np.zeros([nmode]);
    eof    = np.zeros([N, LAT*LON]);
    EOFs   = np.zeros([LAT, LON, nmode]);
    # Eigvector analysis
    L, E = np.linalg.eig(R)
    # Get modes
    E    = np.dot(H.T, E)
    #sq   = (np.sqrt(np.diag(L))).T
    #sq   = sq[0, :]
    sq = np.sqrt(L)
    E    = E/sq
    Z    = np.dot(E.T, H.T)
    for i in range(nmode):
        eof[i, :] = np.squeeze(E[:, i]).T
        PC[i, :]  = np.squeeze(Z[i, :])
    # Get expvar
    L = np.abs(L)
    dsum = np.sum(np.abs(L))
    # Output
    for i in range(nmode):
        expvar[i] = L[i]/dsum
        EOFs[:, :, i] = np.reshape(eof[i, :], [LAT, LON])
    if reverse==1:
        EOFs, PC = reverse_std(EOFs, PC, nmode)
    return EOFs, PC, expvar

def reverse_std(EOFs, PC, nmode):
    for i in range(nmode):
        STD = np.nanstd(PC[i, :])
        PC[i, :] = PC[i, :]/STD
        EOFs[:, :, i] = EOFs[:, :, i]*STD
    return EOFs, PC
    
def seasonal_decomp(data, method=0):
    '''
    =======================================================================
    Remove the seasonal cycle from 1D data
                            ----- created on 2015/06/15, Yingkai (Kyle) Sha
    -----------------------------------------------------------------------
        data = seasonal_decomp(...)
    -----------------------------------------------------------------------
    Input:
            data
            method: removal done by anomaly (=0) or normalize (=1)
    ======================================================================= 
    '''
    for mon in range(12):
        temp_data = data[mon:len(data):12]
        if method == 0:
            data[mon:len(data):12] = data[mon:len(data):12]-np.nanmean(temp_data)
        elif method == 1:
            data[mon:len(data):12] = (data[mon:len(data):12]-np.nanmean(temp_data))/np.nanstd(temp_data)
    return data
    
def seasonal_decomp3d(data, method=0):
    '''
    =======================================================================
    Remove the seasonal cycle from 3D data
                            ----- created on 2015/06/15, Yingkai (Kyle) Sha
    -----------------------------------------------------------------------
        data = seasonal_decomp(...)
    -----------------------------------------------------------------------
    Input:
            data
            method: removal done by anomaly (=0) or normalize (=1)
    ======================================================================= 
    '''
    for i in range(np.size(data, 0)):
        for j in range(np.size(data, 1)):
            data[i, j, :] = seasonal_decomp(data[i, j, :], method=0)
    return data