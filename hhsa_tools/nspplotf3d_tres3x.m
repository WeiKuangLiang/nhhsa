function [All_nt,fscale,Fscale]=nspplotf3d_tres3x(fm,FM,AM,ntp0,ntp1,fres,Fres,fw0,fw1,Fw0,Fw1,tres,S)

% Input-
%	fm	    - 2-D matrix that specifies the frequency values
%	FM      - 3-D matrix that specifies the frequency values of envelope 
%	AM      - 3-D matrix that specifies the amplitude values of envelope 
%	ntp0	- the start time (point)
%	ntp1	- the end time  (point)
%	fres	- the frequency resolution
%	Fres	- the frequency resolution of AM (envelope) 
%	fw0	    - the minimum frequency
%	fw1	    - the maximum frequency
%	Fw0	    - the minimum frequency of AM 
%	Fw1	    - the maximum frequency of AM 
%   tres    - the time resolution
% Output-

%	All_nt	- 4-D matrix of the HHS spectrum, where
%		    1st dimension specifies the number of frequencies,
%		    2nd dimension specifies the number of frequencies of envelope
%           3rd dimension specifies the number of time points
%           4th dimension what number IMF
%	fscale	- vector that specifies the frequency-axis values
%	Fscale	- vector that specifies the AM frequency-axis values


%----- Check the input arguments
if nargin<3
    error('nspplot: both frequency and amplitude matrices required');
end

%----- Specify default variables
if nargin < 5
    ntp0 = [];
    ntp1 = [];
end
if nargin<7
    fres=[];
    Fres=[];
end
if nargin<9
    fw0=[];
    fw1=[];
end
if nargin<11
    Fw0=[];
    Fw1=[];
end
if nargin<12
    tres=[];
end
if nargin<13
    S=[];
end
[npt,nimf]=size(fm);
%----- Initialize default variables
if isempty(ntp0)
    ntp0=1;
    ntp1=npt;
end
if isempty(S)
    S.dyadic=0;
    S.collapse=0;
end
% if ~isfield(S,'CW')
%     S.CW=[];
% end
% if isempty(S.CW)
%     S.CW=0;
% end
dyadic=S.dyadic;
collapse=S.collapse;
if dyadic==0
    if isempty(fres)
        fres=129;
    end
    if isempty(Fres)
        Fres=65;
    end
    if isempty(fw0)
        fw0=0;
        fw1=64;
    end
    if isempty(Fw0)
        Fw0=0;
        Fw1=32;
    end
    if isempty(tres)
        tres=npt;
    end
    if fw0<0
        fw0=0;
        disp('Without log scaling, freq lower bound cannot be less than 0! Set to 0 for you!');
    end
    if Fw0<0
        Fw0=0;
        disp('Without log scaling, AM Freq lower bound cannot be less than 0! Set to 0 for you!');
    end
elseif dyadic==1
    btwFreqs=S.dyad_btw;
    if isempty(fw0)
        fw0=-1;
        fw1=6;
    end
    if isempty(Fw0)
        Fw0=-1;
        Fw1=5;
    end
    if isempty(fres)
        fres=(fw1-fw0)*btwFreqs+1;
    end
    if isempty(Fres)
        Fres=(Fw1-Fw0)*btwFreqs+1;
    end
    
    if isempty(tres)
        tres=npt;
    end
end


fw=fw1-fw0;
Fw=Fw1-Fw0;
t=ceil((1:npt)*tres/npt); %t is the mapping position of time values into the time axis grid
%snt = zeros(fres,Freq);
%TNM = nimf;
%Tlen = npt;

if dyadic==0 && Fw0==0
    collapse=0;
    disp('The lower bound of AM has been set to zero, no need to collapsing data!');
end
if collapse
       All_nt=zeros(fres,Fres+1,tres,nimf); % the addtional 1 added to Fres is for DC_AM
else
       All_nt=zeros(fres,Fres,tres,nimf);
end

if collapse
    for i_imf = 1:nimf
        %nt=zeros(fres,Freq);
        f = fm(:,i_imf);
        A = AM(:,:,i_imf);
        F = FM(:,:,i_imf);
        if dyadic==0
            pf = round((fres-1)*(f-fw0)/fw)+1;
            pF=round((Fres-1)*(F-Fw0)/Fw)+1;
            pF0=round((Fres-1)*(0-Fw0)/Fw)+1;
        elseif dyadic==1
            pf=NaN(size(f));
            pF=NaN(size(F));
            f_ex0=(f>=0);
            F_ex0=(F>=0);
            pf(f_ex0) = round((fres-1)*(log2(f(f_ex0))-fw0)/fw)+1;
            pF(F_ex0)=round((Fres-1)*(log2(F(F_ex0))-Fw0)/Fw)+1;
            pF0=-Inf;
        end

        [ntp2,TNM2]=size(F);
        X=ntp0:ntp1;
        freqidx = pf(X);
        for ii_imf=1:TNM2
            Freqidx=pF(X,ii_imf);
            crit_freq=(freqidx >= 1 & freqidx <= fres);
            crit_Freq=(Freqidx >= 1 & Freqidx <= Fres);
            crit_cross=f(X) > F(X,ii_imf);
            crit_total=crit_Freq & crit_freq & crit_cross;
            if dyadic==0
                crit_totalc= (Freqidx < 1) & (Freqidx > pF0) & crit_freq & crit_cross;
            else
                crit_totalc= (Freqidx < 1) & (~isinf(Freqidx)) & crit_freq & crit_cross;
            end
            if dyadic==0
                crit_totaldc= (Freqidx < 1) & (Freqidx == pF0) & crit_freq & crit_cross;
            else
                crit_totaldc= (Freqidx < 1) & isinf(Freqidx) & crit_freq & crit_cross;
            end
            tX=t(X-ntp0+1);
            for tid=1:length(X)
                if crit_total(tid)
                   All_nt(freqidx(tid),Freqidx(tid)+1,tX(tid),i_imf)=All_nt(freqidx(tid),Freqidx(tid)+1,tX(tid),i_imf)+ 0.5*A(X(tid),ii_imf)^2; % 0.5 is given after 20190408
                elseif crit_totalc(tid)
                   All_nt(freqidx(tid),2,tX(tid),i_imf)=All_nt(freqidx(tid),2,tX(tid),i_imf)+ 0.5*A(X(tid),ii_imf)^2;
                elseif crit_totaldc(tid)
                   All_nt(freqidx(tid),1,tX(tid),i_imf)=All_nt(freqidx(tid),1,tX(tid),i_imf)+ A(X(tid),ii_imf)^2;   
                end
            end
        end
    end
else
    for i_imf = 1:nimf
        %nt=zeros(fres,Freq);
        f = fm(:,i_imf);
        A = AM(:,:,i_imf);
        F = FM(:,:,i_imf);
        if dyadic==0
            pf = round((fres-1)*(f-fw0)/fw)+1;
            pF=round((Fres-1)*(F-Fw0)/Fw)+1;
        elseif dyadic==1
            pf=NaN(size(f));
            pF=NaN(size(F));
            f_ex0=(f>=0);
            F_ex0=(F>=0);
            pf(f_ex0) = round((fres-1)*(log2(f(f_ex0))-fw0)/fw)+1;
            pF(F_ex0)=round((Fres-1)*(log2(F(F_ex0))-Fw0)/Fw)+1;
        end

        [ntp2,TNM2]=size(F);
        X=ntp0:ntp1;
        freqidx = pf(X);
        for ii_imf=1:TNM2
            Freqidx=pF(X,ii_imf);
            crit_freq=(freqidx >= 1 & freqidx <= fres);
            crit_Freq=(Freqidx >= 1 & Freqidx <= Fres);
            crit_cross=f(X) > F(X,ii_imf);
            crit_total=crit_Freq & crit_freq & crit_cross;
            %crit_totalc=(Freqidx < 1 & collapse==1) & crit_freq & crit_cross;
            tX=t(X-ntp0+1);
            for tid=1:length(X)
                if crit_total(tid)
                   All_nt(freqidx(tid),Freqidx(tid),tX(tid),i_imf)=All_nt(freqidx(tid),Freqidx(tid),tX(tid),i_imf)+ 0.5*A(X(tid),ii_imf)^2; % 0.5 is given after 20190408
                end
            end
        end
    end

end
Fscale=linspace(Fw0,Fw1,Fres)';
fscale=linspace(fw0,fw1,fres)';
%snt=flipud(rot90((snt)));