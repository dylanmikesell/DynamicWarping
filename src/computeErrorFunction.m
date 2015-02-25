function err = computeErrorFunction( u1, u0, nSample, lag, norm )
%
% USAGE: err = computeErrorFunction( u1, u0, nSample, lag )
%
% INPUT:
%   u1      = trace that we want to warp; size = (nsamp,1)
%   u0      = reference trace to compare with: size = (nsamp,1)
%   nSample = numer of points to compare in the traces
%   lag     = maximum lag in sample number to search
%   norm    = 'L2' or 'L1' (default is 'L2')
% OUTPUT:
%    err = the 2D error function; size = (nsamp,2*lag+1)
%
% The error function is equation 1 in Hale, 2013. You could umcomment the
% L1 norm and comment the L2 norm if you want on Line 29
%
% Original by Di Yang
% Last modified by Dylan Mikesell (25 Feb. 2015)

if nargin < 5
    norm = 'L2'; % set default
end 

if lag >= nSample
    error('computeErrorFunction:lagProblem','lag must be smaller than nSample');
end

% Allocate error function variable
err = zeros( nSample, 2 * lag + 1 );

%--------------------------------------------------------------------------
% initial error calculation
for ll = -lag : lag % loop over lags
    
    thisLag = ll + lag + 1;
    
    for ii = 1 : nSample % loop over samples
        
        if ( ii + ll >= 1 && ii + ll <= nSample ) % skip corners for now, we will come back to these
            
            diff = u1( ii ) - u0( ii + ll ); % sample difference
            
            switch norm
                case 'L2'
                    err( ii, thisLag ) = ( diff )^2; % difference squared error
                case 'L1'
                    err( ii, thisLag ) = abs( diff ); % absolute value error
            end
            
        end
        
    end
    
end

%--------------------------------------------------------------------------
% Now fix corners with constant extrapolation
for ll = -lag : lag % loop over lags
    
    thisLag = ll + lag + 1;
    
    for ii = 1 : nSample % loop over samples
        
        if ( ii + ll < 1 ) % lower left corner (negative lag, early time)
            
            err( ii, thisLag ) = err( -ll + 1, thisLag );
            
        elseif ( ii + ll > nSample ) % upper right corner (positive lag, late time)
            
            err( ii, thisLag ) = err( nSample - ll, thisLag );
            
        end
        
    end
    
end

return