%  function [spmax, spmin, flag]= extrema(in_data)
%
% This is a utility program for cubic spline envelope,
%   the code is to  find out max values and max positions
%                            min values and min positions
%    (then use matlab function spline to form the spline)
%
%   function [spmax, spmin, flag]= extrema(in_data)
%
% INPUT:
%       in_data: Inputted data, a time series to be sifted;
% OUTPUT:
%       spmax: The locations (col 1) of the maxima and its corresponding
%              values (col 2)
%       spmin: The locations (col 1) of the minima and its corresponding
%              values (col 2)
%
% NOTE:
%      EMD uses Cubic Spline to be the Maximun and Minimum Envelope for
%        the data.Besides finding spline,end points should be noticed. 
%
%References:  ? which paper?
% 
%
%
% code writer: Zhaohua Wu. 
% footnote:S.C.Su
%
% There are two seperste loops in this code .
% part1.-- find out max values and max positions 
%          process the start point and end point  
% part2.-- find out min values and max positions 
%          process the start point and end point  
% Those parts are similar.
%
% Association:eemd.m
% this function ususally used for finding spline envelope
%
% Concerned function: no
%                     (all matlab internal function)

function [spmax, spmin]= extrema_x(in_data)

%flag=1;
m=length(in_data);
if size(in_data,1)>size(in_data,2)
    x=in_data';
else
    x=in_data;
end
%part1.--find local max value and do end process

%start point 
%spmax(1,1)-the first 1 means first point max value,the second 1 means first index
%spmax(1,2)-the first 1 means first point max value,the second 2 means first index
%spmax(1,1)-for position of max 
%spmax(1,2)-for value    of max
d = diff(x);

n = length(d);
d1 = d(1:n-1);
d2 = d(2:n);
% block to find extrema
extrbool=(d1.*d2)<0;
extrbool_max=extrbool & (d1>0);
extrbool_min=extrbool & (d1<0);
if any(extrbool_max)
   indmax = find(extrbool_max)+1;
else
    indmax=zeros(1,0);
end
if any(extrbool_min)
   indmin = find(extrbool_min)+1;
else
    indmin=zeros(1,0);
end
if any(d==0)
    bad = (d==0);
    dd = diff([0 bad 0]);
    debs = find(dd == 1);
    fins = find(dd == -1);
    if debs(1) == 1
        if length(debs) > 1
            debs = debs(2:end);
            fins = fins(2:end);
        else
            debs = zeros(1,0);
            fins = zeros(1,0);
        end
    end
    if ~isempty(debs)
        if fins(end) == m
            if length(debs) > 1
                debs = debs(1:(end-1));
                fins = fins(1:(end-1));
            else
                debs = zeros(1,0);
                fins = zeros(1,0);
            end
        end
    end
    lc = length(debs);
    if lc > 0
        imax0=zeros(1,lc);
        imin0=zeros(1,lc);
        for k = 1:lc
            if d(debs(k)-1) > 0 && d(fins(k)) < 0
               %imax = [imax round((fins(k)+debs(k))/2)];
               imax0(k)= round((fins(k)+debs(k))/2);    
            end 
            if d(debs(k)-1) < 0 && d(fins(k)) > 0 
               %imin = [imin round((fins(k)+debs(k))/2)];
               imin0(k)= round((fins(k)+debs(k))/2); 
            end
        end
    
        imax=imax0(imax0>0);
        imin=imin0(imin0>0);
        if ~isempty(imax)
            indmax = sort([indmax imax]);
        end
        if ~isempty(imin)
            indmin = sort([indmin imin]);
        end
    end
end
% end of block "find extrema"
spmax=zeros(length(indmax)+2, 2);
spmin=zeros(length(indmin)+2, 2);
spmax(1,1) = 1;
spmax(1,2) = x(1);
spmax(end,1) = m;
spmax(end,2) = x(m);
spmax(2:end-1,1)=indmax';
spmax(2:end-1,2)=x(indmax)';
spmin(1,1) = 1;
spmin(1,2) = x(1);
spmin(end,1) = m;
spmin(end,2) = x(m);
spmin(2:end-1,1)=indmin';
spmin(2:end-1,2)=x(indmin)';

%Loop --start find max by compare the values 
%when [ (the jj th value > than the jj-1 th value ) AND (the jj th value > than the jj+1 th value )
%the value jj is the position of the max
%the value in_data (jj) is the value of the max
%do the loop by index-jj
%after the max value is found,use index -kk to store in the matrix
%kk=1,the start point
%the last value of kk ,the end point 



%End point process-please see reference about spline end effect
%extend the slpoe of neighbor 2 max value ---as extend value
%original value of end point -----as original value
%compare extend and original value 
kk=size(spmax,1);
if kk>=4
    slope1=(spmax(2,2)-spmax(3,2))/(spmax(2,1)-spmax(3,1));
    tmp1=slope1*(spmax(1,1)-spmax(2,1))+spmax(2,2);
    if tmp1>spmax(1,2)
        spmax(1,2)=tmp1;
    end

    slope2=(spmax(kk-1,2)-spmax(kk-2,2))/(spmax(kk-1,1)-spmax(kk-2,1));
    tmp2=slope2*(spmax(kk,1)-spmax(kk-1,1))+spmax(kk-1,2);
    if tmp2>spmax(kk,2)
        spmax(kk,2)=tmp2;
    end
% else
%     flag=-1;
end



%part2.--find local min value and do end process
%the syntax are all similar with part1.
%here-explan with beginning local max-find upper starting envelope
%the end process procedure-find out the neighbor 2 local extrema value
%connect those 2 local extrema and extend the line to the end
%make judgement with 1).line extend value  2).original data value
%the bigger value is chosen for upper envelope end control point

%local min 

kk=size(spmin,1);
if kk>=4
    slope1=(spmin(2,2)-spmin(3,2))/(spmin(2,1)-spmin(3,1));
    tmp1=slope1*(spmin(1,1)-spmin(2,1))+spmin(2,2);
    if tmp1<spmin(1,2)
        spmin(1,2)=tmp1;
    end

    slope2=(spmin(kk-1,2)-spmin(kk-2,2))/(spmin(kk-1,1)-spmin(kk-2,1));
    tmp2=slope2*(spmin(kk,1)-spmin(kk-1,1))+spmin(kk-1,2);
    if tmp2<spmin(kk,2)
        spmin(kk,2)=tmp2;
    end
% else
%     flag=-1;
end

% flag=1;