function [dat,baseline] = ft_preproc_baselinecorrect(dat, begsample, endsample)

% FT_PREPROC_BASELINECORRECT performs a baseline correction, e.g. using the
% prestimulus interval of the data or using the complete data
%
% Use as
%   [dat] = ft_preproc_baselinecorrect(dat, begin, end)
% where
%   dat        data matrix (Nchans X Ntime)
%   begsample  index of the begin sample for the baseline estimate
%   endsample  index of the end sample for the baseline estimate
%
% If no begin and end sample are specified for the baseline estimate, it
% will be estimated on the complete data.
%
if nargin<2 || isempty(begsample)
  begsample = 1;
end
if nargin<3 || isempty(endsample)
  endsample = size(dat,2);
end

[dat,baseline] = ft_preproc_polyremoval(dat, 0, begsample, endsample);

% save this old code because it looks non-trivially optimized
%
% persistent hasbsxfun
% if isempty(hasbsxfun)
%   hasbsxfun = exist('bsxfun', 'builtin')==5;
% end
%
% % determine the size of the data
% [Nchans, Nsamples] = size(dat);
%
% % determine the interval to use for baseline correction
% if nargin<2 || isempty(begsample)
%   begsample = 1;
% end
% if nargin<3 || isempty(endsample)
%   endsample = Nsamples;
% end
%
% % estimate the baseline and subtract it
% baseline = mean(dat(:,begsample:endsample), 2);
%
% % ensure that the data is not represented as integer, otherwise "minus" fails
% dat = double(dat);
%
% if hasbsxfun
%   % it is even faster to do this
%   dat = bsxfun(@minus,dat,baseline);
% else
%   % it is faster to loop over samples than over channels due to the internal memory representation of Matlab
%   % for chan=1:Nchans
%   %  dat(chan,:) = dat(chan,:) - baseline(chan);
%   % end
%
%   for sample=1:Nsamples
%     dat(:,sample) = dat(:,sample) - baseline;
%   end
% end