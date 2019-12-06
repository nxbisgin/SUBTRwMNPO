function [dat,beta,x] = ft_preproc_polyremoval(dat, order, begsample, endsample, flag)

% FT_PREPROC_POLYREMOVAL removed an Nth order polynomal from the data
%
% Use as
%   dat = ft_preproc_polyremoval(dat, order, begsample, endsample, flag)
% where
%   dat        data matrix (Nchans X Ntime)
%   order      the order of the polynomial
%   begsample  index of the begin sample for the estimate of the polynomial
%   endsample  index of the end sample for the estimate of the polynomial
%   flag       optional boolean to specify whether the first order basis
%              vector will zscored prior to computing higher order basis
%              vectors from the first-order basis vector (and the beta
%              weights). This is to avoid numerical problems with the
%              inversion of the covariance when the polynomial is of high
%              order/number of samples is large
%
% If begsample and endsample are not specified, it will use the whole
% window to estimate the polynomial.
%
% For example
%   ft_preproc_polyremoval(dat, 0)
% removes the mean value from each channel and
%   ft_preproc_polyremoval(dat, 1)
% removes the mean and the linear trend.
%

% take the whole segment if begsample and endsample are not specified
if nargin<3 || isempty(begsample)
  begsample = 1;
end
if nargin<4 || isempty(endsample)
  endsample = size(dat,2);
end
if nargin<5 || isempty(order)
  flag = 0;
end

% this does not work on integer data
typ = class(dat);
dat = cast(dat, 'double');

% construct a "time" axis
nsamples = size(dat,2);
basis    = (1:nsamples)-1;
if flag
  basis = standardise(basis);
end

% create a set of basis functions that will be fitted to the data
x = zeros(order+1,nsamples);
for i = 0:order
  x(i+1,:) = basis.^(i);
end

% estimate the contribution of the basis functions
% beta = dat(:,begsample:endsample)/x(:,begsample:endsample); <-this leads to numerical issues, even in simple examples
invxcov = inv(x(:,begsample:endsample)*x(:,begsample:endsample)');
beta    = dat(:,begsample:endsample)*x(:,begsample:endsample)'*invxcov;

% remove the estimated basis functions
dat = dat - beta*x;

% convert back to the original input type
dat = cast(dat, typ);