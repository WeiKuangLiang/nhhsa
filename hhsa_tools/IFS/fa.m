%function [f,a] = fa(data,dt,ifmethod,normmethod,nfilter)
% This version gets rid of the median filter for quadrature, Wei-Kuang Liang
% 2020
% Also, this version provides 'phase' information from hilbertm, quad.
% The function FA computes a frequency and amplitude of data(n,k), where 
% n specifies the length of time series, and k is the number of IMFs.
% The user has a choice to choose the instantaneous frequency and
% normalization methods. Nature of the arcosine method suggests not 
% to use the latter to process the residue component.
% First 2 arguments are required. If not passed the rest is set to have 
% default values: the frequency and amplitude is calculated using Hilbert
% method with spline normalization and no smoothing of data.
%
% Calling sequence-
% [f,a] = fa(data,dt[,ifmethod][,normmethod][,nfilter])
%
% Input-
%	  data		  - 2-D matrix of IMF components
%	  dt		    - time increment per point
%	  ifmethod	- method of determining an instantaneous frequency (and amplitude)
%	  normmethod	- normalization method
%	  nfilter		- number of points to use for filter
% Output-
%	  f		    - 2-D matrix f(n,k) that specifies frequency
%	  a		    - 2-D matrix a(n,k) that specifies amplitude
%
% Ifmethod options:
%    'hilbtm'  :use Hilbert transform, but hilbtm.m is used instead of the standard hilbert.m(function FAimpHilbert )
%                  normalization of input data recommended but not required
%    'quad'    :use quadrature method (function FAquadrature),
%                  normalization of input data required
%
%
% Normmethod options:
%    'none'	   :no normalization, recommended for 'zc' option 
%    'spline'	 :spline normalization, recommended for 'hilbert' and 'acos' options (function splinenormalize )
%               not recommended for Ensemble EMD method due to the possible overshot
%    'pchip'   :cubic hermite spline normalization , recommended for normalization when (function pchipnormalize )
%               using Ensemble EMD
%
%
% Non MATLAB Library routines used in the function are:
%	BLOCKNORMALIZE, SPLINENORMALIZE, HILBERTNORMALIZE, MEDIANFILTER.
%
% Earlier version of this code was written by  
% Kenneth Arnold (NASA GSFC)	Summer 2003, Initial
% Karin Blank (NASA GSFC)       5/17/05, edited to add Quad
% Xianyao Chen (RCADA, FIO)   Sep. 20, 2008
% Sheng-Chung.Su 2009/09/10

% The current version is written by 
% Wei-Kuang Liang (ICN, NCU, Taiwan) Feb. 1, 2020
function [f,a,ph] = fa(data,dt,ifmethod,normmethod,nfilter,ns)
% f_exp output the FM expansion by maskemd
%0.Initial the parameters and default settings
  %----- Define default parameters
  
  if nargin<3
      ifmethod = [];
  end
  if nargin<4
      normmethod = [];
  end
  if nargin<5
      nfilter=[];
  end
  if nargin<6
      ns=[];
  end
  if isempty(ifmethod)
      ifmethod = 'quad';
  end
  if isempty(normmethod)
      if (isequal(ifmethod, 'hilbtm'))
          normmethod = 'spline';
      elseif (isequal(ifmethod, 'quad'))
          normmethod = 'pchip';
      else
          normmethod = 'none';
      end
  end
  if isempty(nfilter)
      nfilter=5;
      medf=0;
  else
      nfilter=round(nfilter);
      medf=1;
  end
  if isempty(ns)
      ns=2;
  end
  
  %----- Get the dimensions
  [npt,nIMF]=size(data);
  flipped=0;
  
  if npt<nIMF
      %----- Flip the data
      data=data';
      flipped=1;
      [npt, nIMF]=size(data);
  end

%1.Normalize the data by specific method
%2009/10/14 Norden E Huang decided to operate normalize procedure 3 times in fa.m
  %----- Normalize data if requested    
 if isequal(normmethod, 'spline')
      [data1,na1]=splinenormalize(data);
      [data2,na2]=splinenormalize(data1);
      [data3,na3]=splinenormalize(data2);
      data=data3;
      na=na1.*na2.*na3;
 elseif isequal(normmethod, 'pchip')
      [data1,na1]=pchipnormalize(data,4);  % there is an additional pchip in FAquadrature
      [data2,na2]=pchipnormalize(data1,4); 
      data=data2;
      na=na1.*na2;
  else
      [data1,na1]=splinenormalize(data);
      [data2,na2]=splinenormalize(data1);
      [data3,na3]=splinenormalize(data2);
      data=data3;
      na=na1.*na2.*na3;
      disp ('fa: unknown normalization method,use default ');
  end % ... use spline normalize as default
  
%2.Calculate the [f,a] by specific method
  %----- Calculate the frequency and amplitude
  if (isequal(ifmethod, 'hilbtm'))
      if nargout>2
          [f, ~,ph] = FAimphilbert(data,dt,0.05,ns);
      else
          [f, ~] = FAimphilbert(data,dt,0.05, ns);
      end
  elseif (isequal(ifmethod, 'quad')) % quad, supporting FM expansion now 
      if nargout>2
          [f,~,ph]  = FAquadrature(data,dt,[],ns);
      else
          [f,~]  = FAquadrature(data,dt,[],ns);
      end
  end
  
  if ~isempty(na)
     %----- Throw away Hilbert etc. amplitude, use normalized amplitude
     a=na;
  end
  
%3.use median-filter to process the I.F values 
  %----- Filter the frequency if requested
%  if (~isequal(ifmethod, 'quad') && ~isequal(ifmethod, 'hilbtm') && ~isequal(ifmethod, 'qzc') && ~isequal(ifmethod, 'zc')) && nfilter>0 % get rid of median filter for quadrature, Wei-Kuang Liang 2017
  if medf>0 && nfilter>0
      for i=1:nIMF
          f(:,i)=medianfilter(f(:,i),nfilter);
          %f(:,i)=medianfilter2(f(:,i),nfilter); % using Kevin's customized median, instead of matlab build-in median
      end
  end

  %----- Flip again if data was flipped at the beginning
  if flipped
      f=f';
      a=a';
      if nargout>2
          ph=ph';
      end
  end