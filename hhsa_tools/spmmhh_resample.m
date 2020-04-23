function [Y,alpha] = spmmhh_resample(X,alpha)
% [Jean:] Basic resample function (when no Signal Proc. Toolbox)
% FORMAT Y = spm_resample(X,alpha)
% IN:
%   - X: a nXm matrix of n time series
%   - alpha: the ration of input versus output sampling frequencies. If
%   alpha>1, rs(X,alpha) performs upsampling of the time series.
% OUT:
%   - Y: nX[alpha*m] matrix of resampled time series
%   - alpha: true alpha used (due to rational rounding)
% This function operates on rows of a signal matrix. This means it can be
% used on a block of channels.
N0     = size(X,2);
N=(N0-1)*alpha+1;
tr=1:N0;
ts=linspace(1,N0,N);
Y = spline(tr,X,ts);