function error = computeDTWerror( Aerr, u, lag0 )
%
% Compute the accumulated error along the warping path for Dynamic
% Time Warping.
%
% USAGE: function error = computeDTWerror( Aerr, u, lag0 )
%
% INPUT:
%   Aerr = error MATRIX (equation 13 in Hale, 2013)
%   u    = warping function (samples) VECTOR
%   lag0 = value of maximum lag (samples) SCALAR
%
% Written by Dylan Mikesell
% Last modified: 25 February 2015

npts = numel(u);

if size(Aerr,1) ~= npts
    disp('Funny things with dimensions of error matrix: check inputs.');
    Aerr = transpose(Aerr);
end

error = 0; % initialize

% accumulate error
for ii = 1:npts
    idx = lag0 + 1 + u(ii); % index of lag
    error = error + Aerr( ii, idx ); % sum error
end

return