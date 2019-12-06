function [OP, ssv] = find_orthogonal_projection(meg,iR,Bt)
% [OP, ssv] = find_orthogonal_projection(meg,iR,Bt)
% OP size(Nchan x Nchan), 
% meg is meg data (Ntime x Nchan)
% iR is the time markers for R waves
% Bt is threshold to stop (optional)
% Fs is the system sample rate

% ssv is the set of vectors that were found.

if ~exist('Bt','var')
    Bt = 1e-13;
end

[~, Nc] = size(meg);
ssv = zeros(Nc,15);
mcg = time_locked_avg(meg,iR);
%Fs=312.5; [mcg, ~] = mcg_mode(mcg,Fs);   % remove statistical mode
OP = eye(length(meg(1,:)));
for n = 1:15
    mcgx = mcg*OP;
    [mx,ix] = max(max(abs(mcgx')));
    A = mcgx(ix,:)';
    ssv(:,n) = A;
    if mx < Bt ; 
        disp(['Number of vectors ' num2str(n-1)])
        %figure(666);plot(mcgx)
        break;
    end
    P = A*inv(A'*A)*A';
    NV = eye(size(P))-P;
    OP = OP*NV;
end
ssv = ssv(:,1:(n-1));