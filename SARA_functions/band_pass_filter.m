function D = band_pass_filter(cf1,cf2,D,sf)
% D = band_pass_filter(cf1,cf2,D,sf)
% cf1: low cut frequency  
% cf2: high cut freq
% D: data [smaples X channels]
% sf: samplig freq
%
% 2013 Diana Escalona

[b,a] = butter(4,[cf1/(sf/2) cf2/(sf/2)]);
for i = 1 : size(D,2)
    D(:,i) = D(:,i) - mean(D(:,i));
end
D = filtfilt(b,a,D);

for i = 1 : size(D,2)
    D(:,i) = D(:,i) - mean(D(:,i));
end

return;