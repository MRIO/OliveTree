function [Qcoef1, Qcoef2, q1, q2]=co_fcplfct2(phi1,phi2,dphi1,dphi2,N,ngrid)
% DAMOCO Toolbox, function CO_FCPLFCT2, version 03.03.14
%
% Given two protophases and their derivatives, the function yields 
% the coupling functions q1(phi1,phi2) and q2(phi2,phi1) 
% via fitting a Fourier series. 
%
% Form of call: 
%          [Qcoef1,Qcoef2,q1,q2]=co_fcplfct2(phi1,phi2,dphi1,dphi2,N,ngrid)
%          [Qcoef1,Qcoef2]=co_fcplfct2(phi1,phi2,dphi1,dphi2,N,ngrid)
%          [Qcoef1,Qcoef2,q1,q2]=co_fcplfct2(phi1,phi2,dphi1,dphi2,N)
% Input:   phi1:        phase of the 1st system,
%          phi2:        phase of the 2nd system ('external'),
%          dphi1:       derivative of the phase of the 1st system
%          dphi2:       derivative of the phase of the 2nd system
%          N:           maximal order of Fourier expansion,
%          ngrid:       size of the grid for function computation, 
%                       by default ngrid =  100
% Output:  Qcoef1,Qcoef2 are Fourier coefficients of the coupling functions
%          q1,q2 are the functions, computed on a grid                
%          
phi1 = unwrap(phi1);    phi2 = unwrap(phi2);
A = zeros(4*N+1, 4*N+1);                           % This matrix contains the coefficients A(n+k),(m+l) 
                                                   % for the linear system of equations for the coefficients Qn,m. 
or2=2*N;  or21=or2+1;  or1=N+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Computing the coefficients of the matrix A using symmetries of the coefficients 
for n = -or2 : or2
    for m = -or2 : n
        A(n+or21, m+or21) =  mean(exp(1i*(n*phi1 + m*phi2) ));
        A(-n+or21, -m+or21)=conj(A(n+or21, m+or21));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Computing the coefficients of the matrices Bnm
B1 = zeros(or21*or21);    % This vector contains the coefficients B1n,m for the linear equation system for
                          % the coefficients Q1n,m of the phase phi1
B2 = B1;                  % This vector contains the coefficients B2n,m for the linear equation system for
                          % the coefficients Q2n,m of the phase phi2                                        
C = B1;                   % The elements of the matrix A are reorganized in C to match the requirements of the 
                          % MATLAB function to solve systems of linear equations
ind=1; 
for n = -N : N
    i1_1=(n+N)*or21; 
    for m = -N : N
        i1=i1_1+m+or1; i4=m+or21;
        tmp=exp(-1i*( n*phi1 + m*phi2) );
        B1(ind)= mean(dphi1.* tmp);
        B2(ind)= mean(dphi2.* tmp);
        ind=ind+1;
        for r = -N : N
            i3=(r+N)*or21 +or1; i2=(n-r)+or21;
            for s = -N : N;     % Elements of the matrix A are reorganized in C 
                C(i1,i3 + s) = A(i2,i4-s);                                                      
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%       Solving the system of linear equations to obtain the coefficients qc_nm
qc1 = conj(C) \ B1;    
qc2 = conj(C) \ B2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%         Reorganizing the Fourier cofficients qc in the matrix Qcoef
Qcoef1 = zeros(or21,or21); Qcoef2 = Qcoef1;
for n = 1 : or21
    k=(n-1)*or21;
    for m = 1 : or21
        Qcoef1(n, m)=qc1(k+m);    
        Qcoef2(n, m)=qc2(k+m);   
    end
end
Qcoef2 = Qcoef2.';      % Reorganizing the matrix Qcoef2 to match our convention
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%       Computing the coupling functions q1,q2 on a grid, if required
if nargout == 4;
    if nargin == 5; ngrid=100; end;     %Default value
    [Y,X]=meshgrid(2*pi*(0:ngrid-1)/(ngrid-1),2*pi*(0:ngrid-1)/(ngrid-1));    
    q1 = zeros(ngrid,ngrid); q2=q1;
    for n = -N : N
        for m = -N : N
            tmp=exp(1i*n*X + 1i*m*Y);
            q1 = q1 + Qcoef1(n+or1, m+or1) * tmp;    
            q2 = q2 + Qcoef2(n+or1, m+or1) * tmp; 
        end
    end
    q1 = real(q1);           % Eliminating imaginary rests (numerical)  
    q2 = real(q2);           
end
end
