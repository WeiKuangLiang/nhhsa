%function  [f,a,ph]=FAquadrature(data,dt)
% This new FAquadrature uses a mask emd filter to avoid the need of median filter in the following
% The usage of mask emd filter to get non-negative IF is proposed by Wei-Kuang Liang on 2017
% The main purpose of Liang's update is to obtain stable phases for further analysis(e.g. FM).

% The function FAquadrature generates a frequency and amplitude using Quadrature method 
% applied to data(n,k), where n specifies the length of time series, 
% and k is the number of IMFs.
% Non MATLAB Library routine used in the function is: FINDCRITICALPOINTS.
%
% Calling sequence-
% [f,a]=faquadrature(data,dt)
%
% Input-
%	  data	- 2-D matrix of IMF components 
%	  dt	    - time increment per point
% Output-
%	  f	    - 2-D matrix f(n,k) that specifies frequency
%	  a	    - 2-D matrix a(n,k) that specifies amplitude
%
%Note:
% Used by- 
%	FA
%History of the code:
% FAquadrature code modified for obtaining stable IF by Wei-Kuang Liang (2017 10)
% fa.m code modification by Xianyao Chen (RCADA, FIO)
% Discussion result about fa.m in 2009 HHT class with Norden E. Huang 
%written by 
% S.C.Su    Sep. 2009(NCU Rcada)combine noh.m and calcqf.m together
%           as a new mfile FAquadrature.m  
%  in order to integrate all m file in a group type         
%footnote: S.C.Su (2009/09/07)
%
% 1. Two extra normalization before quad
% 2.use function noh to form the analytical signal with quadratrue formula
% 3.use function calcqf to calculate the Instantaneous frequency
%

function  [f,a,ph]=FAquadrature(data,dt, ratio_boundary, ns)
% This new quadrature offer the freq expansion
%1. Two extra normalization before quad 
%before calculate I.F with quadrature mathod 
%    Normalize twice to make sure the values lie between +1~-1
    [data,a1]=pchipnormalize(data,4);
    %[data,a2]=pchipnormalize(data,4);
    a=a1;
nIMF=data;
%2.use function noh to form the analytical signal with quadratrue formula
quadrature = noh(nIMF);

%3.use function calcqf to calculate the Instantaneous frequency
%smooth=5;
if nargin < 3
    ratio_boundary=[];
end
if nargin < 4
    ns=[];
end
if isempty(ratio_boundary)
    ratio_boundary=0.05;
end
if isempty(ns)
    ns=2;
end
if nargout>2
    [f,ph]=calcqf(quadrature, dt, ratio_boundary, ns);
else
    f = calcqf(quadrature, dt, ratio_boundary, ns);
end


% 
function quadrature = noh(nIMF)

%1.check input,flip signal if necessary
%----- Get the dimension
[npt,ncol] = size(nIMF);

%----- Initialize and flip data if needed 
flipped=0;
if (ncol > npt)
    flipped=1;
    nIMF=nIMF';
    [npt,ncol] = size(nIMF);
end

%2.Calculate the quadrature value--loop start
 
quadrature=zeros(npt,ncol); 
unIMF = spmmhh_resamplev(nIMF,2); % upsample IMFs to overcome the possible jitter when data~=1/-1, for 2020
df_unIMF=diff(unIMF);
for i=1:ncol
        data = nIMF(:,i);
        mask=zeros(size(data));
        df_udata=df_unIMF(:,i); % for 2020
%3.add mask-for positive/negative signal    
        %creates a mask for the signal 
        %the 1st & 2nd Q = 1 
        %and the 3rd & 4th Q = -1
%%% before 2018        
%         mask(1:end-1) = ((diff(data)>0) .* -2) + 1;
%         mask(end) = mask(end-1);
%%% 2019
%         diffx=diff(data);
%         diffxp=zeros(size(data));
%         diffxp(1)=diffx(1);
%         diffxp(npt)=diffx(npt-1);
%         diffxp(2:(npt-1))=0.5*(diffx(1:(npt-2))+diffx(2:(npt-1)));
%         mask(1:end) = ((diffxp>0) .* -2) + 1;
%%% 2020
        diffxp=zeros(size(data));
        diffxp(1)=df_udata(1)*2;
        diffxp(npt)=df_udata(end)*2;
        tmpdiff=sum(reshape(df_udata(2:(end-1)),2,[]));
        diffxp(2:(npt-1))=tmpdiff';
        mask(1:end) = ((diffxp>0) .* -2) + 1;
%4.the calculation of quadrature value
        y = real(sqrt(1-data.^2));
        %y = abs(sqrt(1-data.^2));
%5.multiplying by mask flips data in 3rd & 4th quadrature
        q = y .* mask;
        quadrature(:,i) = complex(data, q);
        clear mask y q
end
 if flipped==1
    quadrature = quadrature';
 end
%2.Calculate the quadrature value--loop end
%end of noh.m   

%#########################################################################

%function instfreq = calcqf(quadrature, dt, smooth)
%
%
% INPUT:
%        quadrature  - is the complex quadrature signal
%        dt          - is the time increment of the data
%        smooth      - the number of digital values to do median smooth
% OUTPUT: 
%        instfreq    - is the instantaneous frequency of the quadrature signal
%
% NOTE: 
%     this code calculates the instantaneous frequency from the analytic signal
%       Quadrature- a new method to calculate instantaneous frequency propose by Norden E. Huang
%                   There are 3 steps .
%                    1.normalize the IMF, this helps IMF become AM/FM disjointed 
%                    2.form the analytical signal,the normalized-IMF(x), be the real part
%                                                 the quadrature of normalized-IMF(x)- sqrt(1-xx) , be the imaginary part
%                    3.calculate the nstantaneous frequency,use real/imaginary part to find phase angle,and the time derivative is I.F
%       this code is some part of the quadrature calculations
%
% References:   
%  N. E Huang (2008),NCU hht class lecture 
%  6 Instantaneous Frequency.ppt
 
%
function [instfreq, instph, ph_exp] = calcqf(quadrature, dt, ratio_boundary, ns)
%

%1.check input,flip signal if necessary
%----- Get the dimension
if nargin < 3
    ratio_boundary=[];
end
if nargin < 4
    ns=[];
end
if isempty(ratio_boundary)
    ratio_boundary=0.05;
end
if isempty(ns)
    ns=2;
end
[npt,ncol] = size(quadrature);

%----- Initialize and flip data if needed 
flip=0;
if (ncol > npt)
    flip=1;
    quadrature=quadrature';
    [npt,ncol] = size(quadrature);
end
TNM=fix(log2(npt));
instfreq=zeros(npt,ncol);
if nargout>1
   instph=zeros(npt,ncol);
end
if nargout>2
   ph_exp=zeros(npt,TNM,ncol);
end
expt=floor(ratio_boundary*npt);
%2.Calculate the phase and it's derivative  value--loop start
for i=1:ncol
    uwph=unwrap(angle(quadrature(:,i)));
    f=zeros(size(uwph));
    diffuwph=diff(uwph);
    if nargout<3 && all(diffuwph(expt+1:end-expt)>=0) % exclude boundary
        diffph=diffuwph;
        suwph=uwph;
    else
        monoph=linspace(uwph(1),uwph(end),length(uwph))';
        
        imn=cmask_emdn((uwph-monoph),-1,ns);
        rimn=imn; % rimn is for the following manipulation
        scmn = cum_mx( imn); % imf residuals
        cmn=scmn+repmat(monoph,1,size(scmn,2));
        diffcmn=diff(cmn(expt+1:end-expt,:));% exclude boundary
        nonegids=all(diffcmn>=0);
        ind=find(nonegids,1); % get first non-negative residual
        if ~isempty(ind)
            suwph=cmn(:,ind);
            if ind>1
                rimn(:,1:(ind-1))=0.001*imn(:,1:(ind-1)); % set phase imf above ind to ~0
            end
        else
            suwph=monoph;
        end
        diffph=diff(suwph);
    end
    f(1)=diffph(1)/(2*pi*dt);
    f(end)=diffph(end)/(2*pi*dt);
    f(2:end-1) = 0.5*(diffph(1:end-1)+diffph(2:end))./(2*pi*dt);
%     f(1)=diffph(1)/(2*pi*dt);
%     f(end)=diffph(end)/(2*pi*dt);
%     f(1:end-1) = diffph(1:end)/(2*pi*dt);
    instfreq(:,i) = f;
    if nargout>1
       instph(:,i)=suwph;
    end
    if nargout>2
        %rimn=imn;
        rimn(:,end)=imn(:,end)+monoph;
        ph_exp(:,1:size(rimn,2),i)=rimn;
    end
end
%2.Calculate the phase and it's derivative  value--loop end
if flip ==1
   instfreq = instfreq'; 
   if nargout>1
       instph= instph';
   end
end    
%end of calcqf.m
function cmn = cum_mx( imn )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cmn=imn;
for i=1:size(imn,2)
    cmn(:,i)=sum(imn(:,i:end),2);
end

function [Y,alpha] = spmmhh_resamplev(X,alpha)
% vertical version resample function 
% FORMAT Y = spm_resample(X,alpha)
% IN:
%   - X: a nXm matrix of m time series
%   - alpha: the ration of input versus output sampling frequencies. If
%   alpha>1, rs(X,alpha) performs upsampling of the time series.
% OUT:
%   - Y: nX[alpha*m] matrix of resampled time series
%   - alpha: true alpha used (due to rational rounding)
% This function operates on rows of a signal matrix. This means it can be
% used on a block of channels.
N0     = size(X,1);
N=(N0-1)*alpha+1;
tr=1:N0;
ts=linspace(1,N0,N);
X_2lastdim=shiftdim(X,1);
Y_p = spline(tr,X_2lastdim,ts);
Y=permute(Y_p,[ndims(Y_p),1:(ndims(Y_p)-1)]);