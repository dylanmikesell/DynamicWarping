function stbar = backtrackDistanceFunction( dir, d, err, lmin, b )
%
% USAGE: stbar = backtrackDistanceFunction( dir, d, err, lmin, b )
%
% INPUT:
%   dir   = side to start minimization ( dir > 0 = front, dir <= 0 =  back)
%   d     = the 2D distance function; size = (nsamp,2*lag+1)
%   err   = the 2D error function; size = (nsamp,2*lag+1)
%   lmin  = minimum lag to search over
%   b     = strain limit (integer value >= 1)
% OUTPUT:
%   stbar = vector of integer shifts subject to |u(i)-u(i-1)| <= 1/b
%
% The function is equation 2 in Hale, 2013.
%
% Original by Di Yang
% Last modified by Dylan Mikesell (19 Dec. 2014)

nSample = size(d,1); % number of samples
nLag    = size(d,2); % number of lags
stbar   = zeros(1,nSample); % allocate

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
% start from the end (front or back)
[~, ll ]      = min( d( iBegin, : ) ); % find minimum accumulated distance at front or back depending on 'dir'
stbar(iBegin) = ll + lmin - 1; % absolute value of integer shift
%--------------------------------------------------------------------------
% move through all time samples in forward or backward direction
ii = iBegin;
while (ii ~= iEnd)
    
    % min/max for edges/boundaries
    ji = max( 1, min( [ nSample, ii + iInc ] ) );
    jb = max( 1, min( [ nSample, ii + iInc * b ] ) );
    
    % -----------------------------------------------------------------
    % check limits on lag indices
    lMinus1 = ll - 1; % lag at l-1
    
    if lMinus1 < 1 % check lag index is greater than 1
        lMinus1 = 1; % make lag = first lag
    end
    
    lPlus1 = ll + 1; % lag at l+1
    
    if lPlus1 > nLag % check lag index less than max lag
        lPlus1 = nLag; % D.M. and D.Y. version
    end
    % -----------------------------------------------------------------
    
    % get distance at lags (ll-1, ll, ll+1)
    distLminus1 = d( jb, lMinus1 ); % minus:  d( i-b, j-1 )
    distL       = d( ji, ll );      % actual: d( i-1, j   )
    distLplus1  = d( jb, lPlus1 );  % plus:   d( i-b, j+1 )
    
    if (ji ~= jb) % equation 10 in Hale (2013)
        for kb = ji : iInc : jb - iInc % sum errors over i-1:i-b+1
            distLminus1 = distLminus1 + err( kb, lMinus1 );
            distLplus1  = distLplus1  + err( kb, lPlus1  );
        end
    end
    
    dl = min( [ distLminus1, distL, distLplus1 ] ); % update minimum distance to previous sample
    
    if ( dl ~= distL ) % then ll ~= ll and we check forward and backward
        if ( dl == distLminus1 )
            ll = lMinus1;
        else % ( dl == lPlus1 )
            ll = lPlus1;
        end
    end
    
    % assume ii = ii - 1
    ii = ii + iInc; % previous time sample
    
    stbar(ii) = ll + lmin - 1; % absolute integer of lag
    
    % now move to correct time index, if smoothing difference over many
    % time samples using 'b'
    if ( ll == lMinus1 || ll == lPlus1 ) % check edges to see about b values
        if ( ji ~= jb ) % if b>1 then need to move more steps
            for kb = ji : iInc : jb - iInc 
                ii = ii + iInc; % move from i-1:i-b-1
                stbar(ii) = ll + lmin - 1; % constant lag over that time
            end
        end
    end
end

return
