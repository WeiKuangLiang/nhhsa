%function [f, a, ph] = FAimphilbert(data,dt)
% This new FAimphilbert uses a mask emd filter to avoid the need of median filter in the following
% The usage of mask emd filter to get non-negative IF is proposed by Wei-Kuang Liang on 2017
% The function FAimphilbert calculates frequency and amplitude
% of data(n,k) using an improved Hilbert method, where n specifies the length
%  of time series, and k is the number of IMF components.
%    hilbt.m is used to perform Hilbert transform instead of hilbert.m,
%     which mainly reduce the Gipps impacts from the end points.
% 
% hilbtm.m is used to perform an improved Hilbert transform.(Gibbs phenomena is fixed)
% Note: FAH allows the instantaneous frequency to be negative. 
%
% Calling sequence-
% [f,a] = FAimphilbert(data,dt)
%
% Input-
%   data	- 2-D matrix data(n,k) of IMF components 
%	  dt	  - sampling period in seconds
% Output-
%	  f	    - 2-D matrix f(n,k) that specifies the Hilbert frequency in Hz
%	  a	    - 2-D matrix a(n,k) that specifies the Hilbert amplitude
%
% Used by-
% 	FA, NSPABMUN.
%
%History:  
% Norden Huang (NASA GSFC) Junw 2,2002 :Initial
% Xianyao Chen, september, 2008 : modified
% S.C.Su ,Sep ,2009, rename all the functions
% Code modified for obtaining stable IF by Wei-Kuang Liang (2017 10)
function [fr, ar, phr] = FAimphilbert(datar,dt,ratio_boundary, ns)

%----- Get the dimension
%[nPoints, nIMF] = size(data);
[nIMF,npt] = size(datar); 
flip=0;  
if nIMF>npt
    data=datar';
    [nIMF,npt] = size(data);
    flip=1;
else
    data=datar;
end   
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
%----- Apply improved Hilbert transform
h=hilbtm(data);

%----- Get the instantaneous frequency
%f(1:nIMF,1:npt)=0;
TNM=fix(log2(npt));
f=zeros(nIMF,npt);
%if nargout>2
   ph=zeros(nIMF,npt);
%end
expt=floor(ratio_boundary*npt);
for i=1:nIMF
    %temp=diff(unwrap(angle(h(j1,:))))./(2*pi*dt);
    uwph=unwrap(angle(h(i,:)));
    monoph=linspace(uwph(1),uwph(end),length(uwph));
    imn=cmask_emdn((uwph-monoph),-1,ns);
    rimn=imn;
    scmn = cum_mx( imn); % imf residuals
    cmn=scmn+repmat(monoph',1,size(scmn,2));
    diffcmn=diff(cmn(expt+1:end-expt,:));% exclude boundary
    nonegids=all(diffcmn>=0);
    indp=find(nonegids,1); % get first non-negative residual
    if ~isempty(indp)
       ind=indp(1);
    else
        ind=-1;
    end
    if ind>0
        suwph=cmn(:,ind)';
        if ind>1
            rimn(:,1:(ind-1))=imn(:,1:(ind-1))*0.001; % set phase imf above ind to ~0
        end
    else
        suwph=monoph;
    end
    inst_f=zeros(size(uwph));
    diffph=diff(suwph);
    inst_f(1)=diffph(1)/(2*pi*dt);
    inst_f(end)=diffph(end)/(2*pi*dt);
    inst_f(2:end-1) = 0.5*(diffph(1:end-1)+diffph(2:end))./(2*pi*dt);
    f(i,:) = inst_f;
%    if nargout>2
       ph(i,:)=suwph(:)';
    %end
end

%----- Get the amplitude
a=abs(h);
if flip ==1
    fr=f';
    ar=a';
%     if nargout>2
       phr=ph';
%     end
else
    fr=f;
    ar=a;
    phr=ph;
end
function cmn = cum_mx( imn )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cmn=imn;
for i=1:size(imn,2)
    cmn(:,i)=sum(imn(:,i:end),2);
end
