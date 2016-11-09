function[Qcoef1, Qcoef2, Qcoef3]=co_fcpltri(phi1, phi2, phi3, Dphi1, Dphi2, Dphi3, N)
% DAMOCO Toolbox, function CO_FCPLTRI, version 06.03.14
% Given phases of 3 oscillators and their derivatives, the function yields the coupling functions 
%       q1(phi1, phi2, phi3), 
%       q2(phi2, phi3, phi1), and
%       q3(phi3, phi1, phi2)
% via fitting a Fourier series: exp( i(n*phi1 + m*phi2 +k*phi3) ). 
% The results for each oscillator are given by 3-dimensional matrix of coefficients
%       Qcoef1(n,m,k)
%       Qcoef2(m,k,n)
%       Qcoef3(k,n,m)
% where the indices are associated with the phases as indicated above:
%       n to ph1, 
%       m to phi2
%       k to phi3
% Hence, the first dimension of each matrix always corresponds to the phase, which dynamics
% is represented by the coefficients.
% Form of call: 
%          [Qcoef1, Qcoef2, Qcoef3] = co_fcpltri(phi1, phi2, phi3, or, safr)
%
% Input:   phi1:    phase of the 1st oscillator,
%          phi2:    phase of the 2nd oscillator,
%          phi3:    phase of the 3rd oscillator,
%          Dphi1:   derivative of phi1
%          Dphi2:   derivative of phi2
%          Dphi3:   derivative of phi3
%          N:       order of the Fourier series
% Output:  Qcoef1, Qcoef2, Qcoef3 are Fourier coefficients of the coupling functions
phi1 = unwrap(phi1); % unwrapped phases are needed here
phi2 = unwrap(phi2);
phi3 = unwrap(phi3);
 
A=zeros(4*N+1,4*N+1,4*N+1);  % This matrix contains the coefficients A(n+p),(m+q),(k+r) 
                                % for the linear system of equations to
                                % obtain the coefficients Qcoef(n,m,k)
ncf=2*N+1; ncf2=ncf*ncf;       % number of coefficients in each dimension

for n = -2*N : 2*N;
    for m = - 2*N : 2*N;
        for k = -2*N : m;       
            A( n+ncf, m+ncf, k+ncf ) = mean( exp( 1i*( n*phi1 + m*phi2 + k*phi3 ) ) );
            A( -n+ncf, -m+ncf, -k+ncf ) = conj (A( n+ncf, m+ncf, k+ncf ));
        end;
    end;
end;

B1 = zeros(ncf^3, 1);          % This vector contains the coefficients B1(n,m,k) for the linear equation system to
                               % obtain the coefficients Qcoef1(n,m,k) of the phase dynamics of phi1
B2 = B1;                       % This vector contains the coefficients B2(n,m,k) for the linear equation system to
                               % obtain the coefficients Qcoef2(m,k,n) of
                               % the phase dynamics of phi2
B3 = B1;                       % This vector contains the coefficients B3(n,m,k) for the linear equation system to
                               % obtain the coefficients Qcoef3(k,n,m) of the phase dynamics of phi2
C = zeros(ncf^3,ncf^3);        % The elements of the matrix A are reorganized in C to match the requirements of the 
                               % MATLAB function to solve systems of linear equations
ind = 1;
for r = -N : N;
    for s = -N : N;
        for q = -N : N;
            EXP=exp( -1i*( r*phi1 + s*phi2 + q*phi3) );
            B1(ind) = mean( Dphi1 .* EXP );
            B2(ind) = mean( Dphi2 .* EXP );
            B3(ind) = mean( Dphi3 .* EXP );
            clear EXP;
            ind = ind+1;
            for n = -N : N;
                for m = -N : N;
                    for k = -N : N;
                        C((r+N)*ncf2 + (s+N)*ncf + (q+N+1) , (n+N)*ncf2 + (m+N)*ncf + (k+N+1) ) = A( (n-r)+ncf, (m-s)+ncf, (k-q)+ncf );
                        % Elements of the matrix A are reorganized in C 
                    end;
                end;
            end;
        end;
    end;
end;
clear A phi* Dphi*
coeff1 = C\B1; % solving linear system to obtain coefficients
coeff2 = C\B2; % solving linear system to obtain coefficients
coeff3 = C\B3; % solving linear system to obtain coefficients
clear C B*
%%% rewriting the coefficients in NxMxK matrix-form
Qcoef1 = zeros(ncf,ncf,ncf);  Qcoef2 = Qcoef1;  Qcoef3 = Qcoef1; 
for n = 1 : ncf;
    for m = 1 : ncf;
        for k = 1 : ncf;
            ind=(n-1)*ncf2 + (m-1)*ncf + k;
            Qcoef1(n, m, k) = coeff1(ind); 
            Qcoef2(m, k, n) = coeff2(ind); 
            Qcoef3(k, n, m) = coeff3(ind); 
        end;
    end;
end;
end