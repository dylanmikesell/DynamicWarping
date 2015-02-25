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

%% make linear time shifts

shift1 = 15; % first shift
shift2 = 10; % second shift

u1           = u0 .* 0; % allocate shifted trace
st           = u1; % allocate a shift trace to compare with estimated shifts from DTW

u1(1:250)    = shiftr( u0(1:250), 1, shift1 ); % shift first part of trace
st(1:250)    = -shift1; % shift right is a delay --> so negative shift

u1(251:npts) = shiftr( u0(251:npts), 1, shift2 ); % shift second part of trace
st(251:npts) = -shift2;

save( 'stepShiftData.mat', 'dt', 'u0', 'u1', 'st', 'shift1', 'shift2' );

