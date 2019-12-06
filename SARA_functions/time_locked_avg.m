function [avg_sig,span] = time_locked_avg(data, index, span_range)
% avg_sig = time_locked_avg(data, index, span_range)
%
% data = input waveform(s). Vector time series or MxN matrix
%   M is time series, N is channel number
%
% index = time index (R markers)
%
% span_range = [lead lag] factors applied to meadian RR interval
%   span_range is optional, default is [0.4 0.6] and places R wave at 40%
% 
% If length(span_range)>2, then span is taken as explicit channel 
%   relative to the peak indicated by 'index'.
%  [-2 -1 0 1 2] indicates five channels.
% avg_sig = averaged waveform of size KxN
% K set by median RR interval, N is set by N of data size
%
% JD Wilson  9-3-2008

error(nargchk(2,3,nargin)); %check number of inputs inputs
if(nargin == 2)
    span_range = [0.4 0.6];
end
switch length(span_range)
    case 2
        lead = span_range(1);
        lag = span_range(2);
        span = floor(median(diff(index)));
        span = floor(-span*lead):floor(span*lag);
    otherwise
        span = span_range;
end

[~, col] = size(data);
avg_sig = zeros(length(span),col);
Nd = length(data(:,1)); % length of data
while((index(1) + span(1)) < 1)
    index = index(2:end);   % fix index less than 1
end
while((index(end) + span(end)) > Nd)
    index = index(1:end-1); % fix index greater than length of data
end
N = length(index);
for nc = 1:col
    for n = 1:N
        avg_sig(:,nc) = avg_sig(:,nc) + data(index(n)+span,nc);
    end
end
avg_sig = avg_sig/N;    % return average waveform(s)