clear all
close all
clc

% This script runs an example of dynamic time warping as an introduction to
% the method.
example = 1 ; % 1=sine wave, 2=square wave shift functions

% you can try two test cases below. Watch the 'lag' parameter. The two st
% functions have different magnitudes so if you don't go to large enough
% lags in the sine model, you won't converge to the correction solution.

switch example
    case 1
load('exampleData/stepShiftData.mat'); % data with constant shifts and step in the middle
lag  = 40; % max nuber of points to search forward and backward (can be npts, just takes longer and is unrealistic)
    case 2
load('exampleData/sineShiftData.mat'); % data with sine wave shifts of one cycle
lag  = 80; % max nuber of points to search forward and backward (can be npts, just takes longer and is unrealistic)
% in the case of the sine wave, you must have b=1 because the strains are
% large compared to dt!! Play with 'b' to see how you get stuck in the
% wrong minimum because you're not allowing large enough jumps.
end

idx = 1:500;

u0 = u0(idx);
u1 = u1(idx);

npts = numel(u0); 
lag = npts-1;

figure;
plot(u0); hold on; plot(u1); grid on;
xlabel('Sample No.'); ylabel('Amplitude [a.u.]'); title('Traces');
legend('u_{0}(t)','u(t)'); legend boxoff;

%% setup warping preliminaries

lvec = -lag:lag; % for plotting below

err = computeErrorFunction2( u1, u0, npts, lag ); % cpmpute error function over lages
% note the error function is independent of strain limit 'b'.

figure;
imagesc(err'); colormap('gray'); c = colorbar;
xlabel('Sample No.'); ylabel('Lag'); title('Error function');

%% Apply dynamic time warping in the forward direction

b = 1;  % impose a strain limit: 1 dt is b=1, half dt is b=2, etc.
% if you mess with b you will see the results in the real shifts compared
% to the esimated plots in the subplot below!

dist  = accumulateErrorFunction2( -1, err, npts, lag, b ); % forward accumulation to make distance function

%% backtrack

stbar = backtrackDistanceFunction( 1, dist, err, -lag, b ); % find shifts
% it is instructive to flip the sign of +/-1 here to see how the function
% changes as we start the backtracking on different sides of the traces.
% Also change 'b' to see how this influences the solution for stbar. You
% want to make sure you're doing things in the proper directions in each
% step!!!

%% plot the results

tvec  = ( 0 : npts-1 ) .* dt; % make the time axis
tvec2 = tvec + stbar .* dt; % make the warped time axis

figure;
% plot the distance function
subplot(1,2,1);
imagesc(tvec,lvec,dist'); 
title('Distance function'); c=colorbar; axis('tight'); axis xy;
xlabel('Time (s)'); ylabel('Lag (nsamp)');
% plot real shifts against estimated shifts
subplot(1,2,2);
plot(tvec,st(idx),'ko'); hold on;
plot(tvec,stbar,'r*'); legend('Actual','Estimated'); legend boxoff;
title('Estimated shifts'); 
%axis([tvec(1), tvec(end), -shift1-5, 0]);
xlabel('Time (s)'); ylabel('Lag (nsamp)');

figure;
% plot input traces
subplot(1,2,1);
plot(tvec,u0,'b',tvec,u1,'r--'); legend('Raw','Shifted'); legend boxoff;
title('Input traces for dynamic time warping')
xlabel('Time (s)'); ylabel('Amplitude (a.u.)');
% plot warped trace to compare
subplot(1,2,2);
plot(tvec,u0,'b'); hold on;
plot(tvec2,u1,'r--'); axis('tight');
legend('Raw','Warped'); legend boxoff;
title('Output traces for dynamic time warping')
xlabel('Time (s)'); ylabel('Amplitude (a.u.)');


%% Apply dynamic time warping in the backward direction

dist  = accumulateErrorFunction( 1, err, npts, lag, b ); % backward accumulation to make distance function
stbar = backtrackDistanceFunction( -1, dist, err, -lag, b ); % find shifts
% it is instructive to flip the sign of +/-1 here to see how the function
% changes as we start the backtracking on different sides of the traces.
% Also change 'b' to see how this influences the solution for stbar

%% plot the results

tvec  = ( 0 : npts-1 ) .* dt; % make the time axis
tvec2 = tvec + stbar .* dt; % make the warped time axis

figure;
% plot the distance function
subplot(1,2,1);
imagesc(tvec,lvec,dist'); 
title('Distance function'); c=colorbar; axis('tight'); axis xy;
xlabel('Time (s)'); ylabel('Lag (nsamp)');
% plot real shifts against estimated shifts
subplot(1,2,2);
plot(tvec,st,'ko'); hold on;
plot(tvec,stbar,'r*'); legend('Actual','Estimated'); legend boxoff;
title('Estimated shifts'); 
% axis([tvec(1), tvec(end), -shift1-5, 0]);
xlabel('Time (s)'); ylabel('Lag (nsamp)');

figure;
% plot input traces
subplot(1,2,1);
plot(tvec,u0,'b',tvec,u1,'r--'); legend('Raw','Shifted'); legend boxoff;
title('Input traces for dynamic time warping')
xlabel('Time (s)'); ylabel('Amplitude (a.u.)');
% plot warped trace to compare
subplot(1,2,2);
plot(tvec,u0,'b'); hold on;
plot(tvec2,u1,'r--'); axis('tight');
legend('Raw','Warped'); legend boxoff;
title('Output traces for dynamic time warping')
xlabel('Time (s)'); ylabel('Amplitude (a.u.)');

%% Apply dynamic time warping in both directions to smooth. (Example in Hale 2013)

dist1 = accumulateErrorFunction( -1, err, npts, lag, b ); % forward accumulation to make distance function
dist2 = accumulateErrorFunction( 1, err, npts, lag, b ); % backwward accumulation to make distance function
dist = dist1 + dist2 - err; % add them and remove 'err' to not counted twice
stbar = backtrackDistanceFunction( -1, dist, err, -lag, b ); % find shifts

% !! Notice now that you can backtrack in either direction and get the same
% result.

%% plot the results

tvec  = ( 0 : npts-1 ) .* dt; % make the time axis
tvec2 = tvec + stbar .* dt; % make the warped time axis

figure;
% plot the distance function
subplot(1,2,1);
imagesc(tvec,lvec,dist'); 
title('Distance function'); c=colorbar; axis('tight'); axis xy;
xlabel('Time (s)'); ylabel('Lag (nsamp)');
% plot real shifts against estimated shifts
subplot(1,2,2);
plot(tvec,st,'ko'); hold on;
plot(tvec,stbar,'r*'); legend('Actual','Estimated'); legend boxoff;
title('Estimated shifts'); 
% axis([tvec(1), tvec(end), -shift1-5, 0]);
xlabel('Time (s)'); ylabel('Lag (nsamp)');

figure;
% plot input traces
subplot(1,2,1);
plot(tvec,u0,'b',tvec,u1,'r--'); legend('Raw','Shifted'); legend boxoff;
title('Input traces for dynamic time warping')
xlabel('Time (s)'); ylabel('Amplitude (a.u.)');
% plot warped trace to compare
subplot(1,2,2);
plot(tvec,u0,'b'); hold on;
plot(tvec2,u1,'r--'); axis('tight');
legend('Raw','Warped'); legend boxoff;
title('Output traces for dynamic time warping')
xlabel('Time (s)'); ylabel('Amplitude (a.u.)');




