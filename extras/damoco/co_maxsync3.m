function[M_SyncIn, maxind, n_theta1, m_theta2, l_theta3]=...
    co_maxsync3(theta1, theta2, theta3, or)
% DAMOCO Toolbox, function CO_MAXSYNC3, version 06.03.14
%
% This function computes a matrix of  n:m:l synchronization indices, 
%           where n=0,...,or and m,l=-or,...,or.
% It also determines the maximal index and corresponding values of n,m,l
% 
% Form of call: [M_SyncIn,maxind,n_theta1,m_theta2,l_theta3]=co_maxsync3(theta1,theta2,theta3,or)
% INPUT     theta1:   (Proto)phase of the 1st oscillator
%           theta2:   (Proto)phase of the 2nd oscillator
%           theta3:   (Proto)phase of the 3rd oscillator
%           or:       Maximal order to be considered 
%   
% OUTPUT    M_SyncIn: 3-dimensional array of  n:m:l synchronization indices, i.e. 
%                     M_SynIn(n,m,l) is the synchronization index of order n:m:l
%           maxind:   maximal synchronization index
%           n_theta1: the value of n, corresponding to the maximal 
%                     synchronization index
%           m_theta2: the value of m, corresponding to the maximal 
%                     synchronization index
%           l_theta3: the value of l, corresponding to the maximal 
%                     synchronization index
M_SyncIn = zeros(or,or,or); maxind=0;
for n = 0 : or        % Computing the 3D-array of synchronizations indices
    for m = -or : or
        for l = -or : or;
            index=co_sync3(theta1, theta2, theta3, n, m, l);  
            
            if index > maxind && (n~=0 || m~=0 || l~=0)
                maxind=index; n_theta1=n; m_theta2=m; l_theta3=l;
            end
            M_SyncIn(n+1, m+or+1, l+or+1) = index;
        end;
    end
end
end