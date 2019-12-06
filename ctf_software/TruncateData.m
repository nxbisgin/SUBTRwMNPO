function [stTrunc,endTrunc]=TruncateData(L,sf,fres)
%    L: data length
% fres: frequency resolution
%   sf: sampling freq 
% stTrunc: sample to truncate from beginning of data
% endTrunc: sample to truncate from end of data
% usage: data(stTrunc:endTrunc,ch) 
if nargin==2
    fres=0.4; %frequency resolution
end
Ttr =1/fres;                                %time resolution
stTrunc= round((Ttr/2)*sf);                 %sample to truncate from beginning of data
endTrunc=L-round(Ttr*sf);                   %sample to truncate from end of data