clear all
close all
clc

% This script makes the trace data to compare a step in the tine shifts

%% make Ricker wavelet

dt     = 0.004; % s
f      = 20; % Hz
[w,tw] = ricker(f,dt);

%% make random reflectivity sequence and convolve with Ricker

npts = 500; % length of time series
f    = randn(1,npts);
u0   = conv(f,w,'same'); % make waveform time series

%% make time varying shifts as sine wave

amp = 0.2 * dt;
dp  = 2 * pi / npts;
p   = 0 : dp : 2 * pi - dp;
st   = amp .* sin(p) ./ dt;

tvec  = ( 0 : npts - 1 ) .* dt;
tvec2 = tvec + st;

u1 = interp1( tvec, u0, tvec2);

st = st./dt;

save( 'sineShiftData.mat', 'dt', 'u0', 'u1', 'st' );

