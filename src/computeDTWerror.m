function error = computeDTWerror( Aerr, u, lag0 )
%
% Compute the accumulated error along the warping path for Dynamic
% Time Warping.
%
% USAGE: function error = computeDTWerror( Aerr, u, lag0 )
%
% INPUT:
%   Aerr = error MATRIX from err_diw (equation 13 in Hale, 2013)
%   u    = warping function (samples) VECTOR
%   lag0 = value of maximum lag (samples) SCALAR
%
%   See also ERR_DIW.
%
% Written by Dylan Mikesell
% Last modified: 21 August 2014

npts = numel(u);

if size(Aerr,1) ~= npts
    disp('Funny things with dimensions of error matrix: check inputs.');
    Aerr = transpose(Aerr);
end

error = 0; % initialize

% accumulate error
for ii = 1:npts
    error = error + Aerr( ii, u(ii) + (lag0+1) );
end

% NOTE: Aerr is computed as the squared difference between points in
% the trace. See err_diw.m for more on this error function

return