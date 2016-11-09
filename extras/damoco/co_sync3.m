function SyncIn=co_sync3(P1,P2,P3,n,m,l)
% DAMOCO Toolbox, function CO_SYNC3, version 06.03.14
%
% This function computes the n:m:l synchronization index from 
% (proto)phases P1, P2 and P3
%
% Form of call SyncIn=co_sync3(P1,P2,P3,n,m,l)
% 
SyncIn = abs(mean(exp(1i*( n*P1 + m*P2 + l*P3) )));
end