# DynamicWarping
Code for Dynamic Time Warping

You can run the run_example.m with two example shifts.
1) A step function of shifts.
2) A sine function of shitts.

In (1) we see that we do not well match the shift in the area the shift 
occurs. This is because the single shift in the step function is larger 
than dt. Therefore is has to spread this shift out over a number of samples.

In (2) we see that we recover the fits well with b=1. If we go to larger b 
values we see that we no longer recover the shifts, because we do not allow
large enough steps.

Read the comments in the script run_example.m to see other things like
forward and backward accumulation operations. 
