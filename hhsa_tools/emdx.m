function allmodep=emdx(Yp,TNMp)
% This is an EMD program
%
% INPUT:
%       Yp: Inputted data;1-d (1*N) data only
%       TNMp: maximum number of IMFs, including residual,if it is less than zero, then automatically determine the number
% OUTPUT:
%       A matrix of N*(m) matrix, where N is the length of the input
%       data Y, and m=fix(log2(N)). 
%part1.read data, find out standard deviation ,devide all data by std
Y=Yp';
xsize=length(Y);
dd=(1:xsize);
%Ystd=zeros(1,1);
Ystdp=std(Y);
Ystd=Ystdp(1);
Y=Y/Ystd;

%part2.evaluate TNM as total IMF number

maxNM=fix(log2(xsize));
%fprintf('maxNM:%f\n',maxNM);
if (TNMp<1) || (TNMp>maxNM)
    TNM=maxNM;
else
    TNM=fix(TNMp);
end
%fprintf('TNM:%f\n',TNM);
allmodep=zeros(xsize,TNM);
%allmode=zeros(xsize,TNM);
xend=Y;
xstart=zeros(size(Y));
%start to find an IMF-----IMF loop start
for nmode =1: (TNM-1)
    xstart(:) = xend(:); 
    %--sift 10 times to get IMF---sift loop  start 
    for iter=1:10
        [spmax, spmin]=extrema_x(xstart');  %call function extrema 
        %the usage of  spline ,please see part11.  
        upper= spline(spmax(:,1),spmax(:,2),dd(:)); %upper spline bound of this sift 
        lower= spline(spmin(:,1),spmin(:,2),dd(:)); %lower spline bound of this sift 
        mean_ul = (upper + lower)/2;%spline mean of upper and lower 
        xstart(:) = xstart(:) - mean_ul(:);%extract spline mean from Xstart
    end
    %--after sift 10 times,that xstart is this time IMF
    allmodep(:,nmode)=xstart(:);
    %part8--subtract IMF from data ,then let the residual xend to start to find next IMF 
    xend(:) = xend(:) - xstart(:);
end
allmodep(:,TNM)=xend(:); % residual
allmodep=allmodep*Ystd;
end

%part11--the syntax of the matlab function spline
%yy= spline(x,y,xx); this means
%x and y are matrixs of n1 points ,use n1 set (x,y) to form the cubic spline
%xx and yy are matrixs of n2 points,we want know the spline value yy(y-axis) in the xx (x-axis)position
%after the spline is formed by n1 points ,find coordinate value on the spline for [xx,yy] --n2 position. 
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
dp = diff(x,1,2);
d=dp(:)';
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
end
