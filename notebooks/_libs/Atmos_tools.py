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
