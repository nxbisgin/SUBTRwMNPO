function D = high_pass_filter(cf1,D,sf)
% cf1: cut frequency  
% D: data [smaples X channels]
% sf: samplig freq
%
% 2013 DEscalona

[b,a] = butter(4,cf1/(sf/2),'high');
for i = 1 : size(D,2)
    D(:,i) = D(:,i) - mean(D(:,i));
end
D = filtfilt(b,a,D);

for i = 1 : size(D,2)
    D(:,i) = D(:,i) - mean(D(:,i));
end

return;