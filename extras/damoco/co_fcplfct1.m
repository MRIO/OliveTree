function [Qcoef, q]=co_fcplfct1(phi1,phi2,dphi1,N,ngrid)
% DAMOCO Toolbox, function CO_FCPLFCT1, version 03.03.14
%
% Given two protophases and the derivatives of the first, the function yields 
% the coupling functions q1(phi1,phi2)
% via fitting a Fourier series. 
%
% Form of call: 
%          [Qcoef, q]   =   co_fcplfct1(phi1,phi2,dphi1,N,ngrid)
%          [Qcoef]      =   co_fcplfct1(phi1,phi2,dphi1,N,ngrid)
%          [Qcoef,q]    =   co_fcplfct1(phi1,phi2,dphi1,N)
% Input:   phi1:        phase of the 1st system,
%          phi2:        phase of the 2nd system ('external'),
%          dphi1:       derivative of the phase of the 1st system
%          N:           maximal order of Fourier expansion,
%          ngrid:       size of the grid for function computation, 
%                       by default ngrid =  100
% Output:  Qcoef are Fourier coefficients of the coupling functions
%          q is the functions, computed on a grid                
%          
phi1 = unwrap(phi1);    phi2 = unwrap(phi2);
A = zeros(4*N+1, 4*N+1);                           % This matrix contains the coefficients A(n+k),(m+l) 
                                                   % for the linear system of equations for the coefficients Qn,m. 
or2=2*N;  or21=or2+1;  or1=N+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Computing the coefficients of the matrix A using symmetries of the cofficients 
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
C = B1;                   % The elements of the matrix A are reorganized in C to match the requirements of the 
                          % MATLAB function to solve systems of linear equations
ind=1; 
for n = -N : N
    i1_1=(n+N)*or21; 
    for m = -N : N
        i1=i1_1+m+or1; i4=m+or21;
        tmp=exp(-1i*( n*phi1 + m*phi2) );
        B1(ind)= mean(dphi1.* tmp);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%         Reorganizing the Fourier cofficients qc in the matrix Qcoef
Qcoef = zeros(or21,or21); 
for n = 1 : or21
    k=(n-1)*or21;
    for m = 1 : or21
        Qcoef(n, m)=qc1(k+m);      
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%       Computing the coupling functions q on a grid, if required
if nargout == 2;
    if nargin == 4; ngrid=100; end;     %Default value
    [Y,X]=meshgrid(2*pi*(0:ngrid-1)/(ngrid-1),2*pi*(0:ngrid-1)/(ngrid-1));    
    q = zeros(ngrid,ngrid);
    for n = -N : N
        for m = -N : N
            tmp=exp(1i*n*X + 1i*m*Y);
            q = q + Qcoef(n+or1, m+or1) * tmp;          
        end
    end
    q = real(q);           % Eliminating imaginary rests (numerical)            
end
end


