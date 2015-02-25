function d = accumulateErrorFunction( dir, err, nSample, lag, b )
%
% USAGE: d = accumulation_diw_mod( dir, err, nSample, lag, b )
%
% INPUT:
%   dir = accumulation direction ( dir > 0 = forward in time, dir <= 0 = backward in time)
%   err = the 2D error function; size = (nsamp,2*lag+1)
%   nSample = numer of points to compare in the traces
%   lag = maximum lag in sample number to search
%   b = strain limit (integer value >= 1)
% OUTPUT:
%    d = the 2D distance function; size = (nsamp,2*lag+1)
% 
% The function is equation 6 in Hale, 2013.
%
% Original by Di Yang
% Last modified by Dylan Mikesell (19 Dec. 2014)

nLag = ( 2 * lag ) + 1; % number of lags from [ -lag : +lag ]

% allocate distance matrix
d = zeros( nSample, nLag ); 

%--------------------------------------------------------------------------
% Setup indices based on forward or backward accumulation direction
%--------------------------------------------------------------------------
if dir > 0            % FORWARD
    iBegin = 1;       % start index
    iEnd   = nSample; % end index
    iInc   = 1;       % increment
else                  % BACKWARD
    iBegin = nSample; % start index
    iEnd   = 1;       % stop index
    iInc   = -1;      % increment
end
%--------------------------------------------------------------------------
% Loop through all times ii in forward or backward direction
for ii = iBegin : iInc : iEnd
    
    ji = max( 1, min( nSample, ii - iInc ) ); % current index i
    
    jb = max( 1, min( nSample, ii - iInc * b ) ); % current lag is*b
    
    % loop through all lags l
    for ll = 1 : nLag
        
        % -----------------------------------------------------------------
        % check limits on lag indices
        lMinus1 = ll - 1; % lag at l-1
        
        if lMinus1 < 1  % check lag index is greater than 1
            lMinus1 = 1; % make lag = first lag
        end
        
        lPlus1 = ll + 1; % lag at l+1
        
        if lPlus1 > nLag % check lag index less than max lag
%             lPlus1 = nLag - 1; % D.Y. version
            lPlus1 = nLag; % D.M. version
        end
        % -----------------------------------------------------------------

        % get distance at lags (ll-1, ll, ll+1)
        distLminus1 = d( jb, lMinus1 ); % minus
        distL       = d( ji, ll );      % actual
        distLplus1  = d( jb, lPlus1 );  % plus
        
        if (ji ~= jb)
            for kb = ji : -iInc : jb + iInc     
                distLminus1 = distLminus1 + err( kb, lMinus1 );
                distLplus1  = distLplus1  + err( kb, lPlus1  );
            end
        end
        
        % equation 6 in Hale (2013) after treating boundaries
        d( ii, ll ) = err( ii, ll ) + min( [ distLminus1, distL, distLplus1 ] ); 
    end
end

return