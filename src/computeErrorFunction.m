function err = computeErrorFunction( u1, u0, nSample, lag )
%
% USAGE: err = computeErrorFunction( u1, u0, nSample, lag )
%
% INPUT:
%   u1      = trace that we want to warp; size = (nsamp,1)
%   u0      = reference trace to compare with: size = (nsamp,1)
%   nSample = numer of points to compare in the traces
%   lag     = maximum lag in sample number to search
% OUTPUT:
%    err = the 2D error function; size = (nsamp,2*lag+1)
%
% The error function is equation 1 in Hale, 2013. You could umcomment the
% L1 norm and comment the L2 norm if you want on Line 29
%
% Original by Di Yang
% Last modified by Dylan Mikesell (19 Dec. 2014)

% Allocate error function variable
err = zeros( nSample, 2 * lag + 1 );

% initial error calculation

for ll = -lag : lag % loop over lags
    
    for ii = 1 : nSample % loop over samples
        
        if ( ii + ll >= 1 && ii + ll <= nSample ) % skip corners for now, we will come back to these
            err( ii, ll + lag + 1 ) = ( u1( ii ) - u0( ii + ll ) )^2; % difference squared error
            % err(ii,ll+lag+1) = abs( u1(ii) - u0(ii+ll) ); % absolute value error
        end
        
    end
    
end

% Now fix corners with constant extrapolation

for ll = -lag : lag % loop over lags
    
    for ii = 1 : nSample % loop over samples
        
        if ( ii + ll < 1 ) % lower left corner (negative lag, early time)
            
            err( ii, ll + lag + 1 ) = err( -ll + 1, ll + lag + 1 );
            
        elseif ( ii + ll > nSample ) % upper right corner (positive lag, late time)
            
            err( ii, ll + lag + 1 ) = err( nSample - ll, ll + lag + 1 );
            
        end
        
    end
    
end

return