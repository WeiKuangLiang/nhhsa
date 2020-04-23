function [C,Cs,P]=cmask_emd3GU(S, TNM, varargin)
% This code implement the enhanced algorithm of EMD(masking EMD)
% History: This code was originally written by Jia-Rong Yeh, followed by an
% improvement by Wei-Kuang Liang
% INPUT : 
%        Y : input signal
%        TNM : number of required imf; if it is less than zero, then automatically determine the number
%
%-------------Additional Input Properties-----------------------------------------
%        
%        odrMask: 0=> 4 phases masking, 1=> 8 phases masking, 2=> 16 phases
%        masking and so on.
%        toFlip : 0=> Original EEMD, References[2] ; 1=> Add anti-phase noise into signal, the way is the same as CEEMD, References[3]
%        numIteration : number of sifting iteration
%        typeSpline : 1=> clamped spline; 2=> not a knot spline;
%        toModify : 0=> None ; 1=> Apply modified linear extrapolation to boundary ; 2 => Mirror Boundary
%        randType : 1=> uniformly distributed white noise; 2=> gaussian white noise
%        seedNo : random seed used for white noise; The value of seed must be an integer between 0 and 2^32 - 1
%        checkSignal : 1=> verify if input signal had NaN or Infinity elements; 0=> Not verify input signal 
%
% OUTPUT :
%        allmode : returned imf
%
%C = rcada_eemd(S,0,1,1);
[S,odrMask,shiftLevel,amp_ratio,ws_r,NoiseLevel,NE,TNM,toFlip,numIteration,typeSpline,toModifyBC,randType,seedNo,IsInputOkay] = parse_checkProperty(S, 0, 1, TNM, varargin);
if shiftLevel>0
    S2 = spmmhh_resample(S,2^shiftLevel);
else
    S2=S;
end
%allmode = emdx(S2, 2);
allmode = rcada_emd(S2, toModifyBC, typeSpline, 2, numIteration);
thd=1; % thd*SD is used for determine the lower threshold

n=size(S2);
t=(1:length(S2))';
if n(1) < n(2)
    S2=S2';
end
cnt=zero_cross_cnt(allmode(:,1)); 
if ws_r<0
   %ws=1.189*pi*cnt/(length(S2)-1);
   ws=1*pi*cnt/(length(S2)-1);
else
   ws=ws_r*pi*cnt/(length(S2)-1); 
end
w=(2^shiftLevel)*ws; 
P.w1=w;
sd=std(allmode(:,1));
th_Amp=thd*sd;
Amp=amp_ratio*(ones(size(allmode,1),1)*th_Amp);
    TNMc=floor(log2(cnt))+shiftLevel; % 20160921
TNMs=TNM+shiftLevel;
Cs=zeros(length(S2),TNMs);
C=zeros(length(S),TNM);
if cnt<=2
    C(:,1)=S;
    return
end
if TNMc>(TNMs-1)
    TNMc=(TNMs-1);
end
for i=1:TNMc
    denominator=2^odrMask;

    for j = 0:(denominator-1)
        %fprintf('Initial phase is %f \n',j/(2*denominator));
        MS1=Amp.*cos(w*t+j*pi/(denominator*2));
        MS2=Amp.*sin(w*t+j*pi/(denominator*2));
           nc=2; 
        %Ctmp=emdx((S2+MS1)', nc);
        Ctmp=rcada_emd((S2+MS1)', toModifyBC, typeSpline, nc, numIteration);
        Cs(:,i)=Cs(:,i)+Ctmp(:,1);
        %Ctmp=emdx((S2-MS1)', nc);
        Ctmp=rcada_emd((S2-MS1)', toModifyBC, typeSpline, nc, numIteration);
        Cs(:,i)=Cs(:,i)+Ctmp(:,1);
        %Ctmp=emdx((S2+MS2)', nc);
        Ctmp=rcada_emd((S2+MS2)', toModifyBC, typeSpline, nc, numIteration);
        Cs(:,i)=Cs(:,i)+Ctmp(:,1);
        %Ctmp=emdx((S2-MS2)', nc);
        Ctmp=rcada_emd((S2-MS2)', toModifyBC, typeSpline, nc, numIteration);
        Cs(:,i)=Cs(:,i)+Ctmp(:,1);
    end
    clear Ctmp;
    w=w/2;
    %Amp=Amp/1.414;
    Cs(:,i)=Cs(:,i)/(4*denominator);
    S2=S2-Cs(:,i);
    allmode = emdx(S2', 2);
    
    sd=std(allmode(:,1));
    th_Amp=thd*sd;
    Amp=amp_ratio*(ones(size(allmode,1),1)*th_Amp);
    
end
   Cs(:,TNMc+1)=S2;
   pC=spmm_downsample(Cs(:,shiftLevel+1:end),2^shiftLevel);
   C(1:size(pC,1),1:size(pC,2))=pC;
   if size(pC,1)<size(C,1)
       redun=size(C,1)-size(pC,1)
       C(size(pC,1)+1:end,:)=repmat(pC(size(pC,1),:),redun,1);
   end
end

function [Y, odrMask, shiftLevel,amp_ratio,ws_r, NoiseLevel, NE, TNM, toFlip, numIteration, typeSpline,toModifyBC,randType,seedNo, IsInputOkay] = parse_checkProperty(Y, NoiseLevel, NE, TNM, varargin)
% Default Parameters
odrMask = 0;
toFlip = 0; % Original EEMD
numIteration = 10; % numIteration = 10
typeSpline = 2;
toModifyBC = 1;
randType = 2;
seedNo = now;
checkSignal = 0;
IsInputOkay = true;
shiftLevel=0;
amp_ratio=0.20;
ws_r=-1;
if(~isempty(varargin{1}))

for iArg = 1 : length(varargin{1});
    
if(iArg == 1)
   odrMask = varargin{1}{iArg};
   if odrMask < 0 && ~isinteger(odrMask)
    fprintf('ERROR : order of Masking signals must be 0 (4) or positive integer n (2^(n+2)).\n');
    IsInputOkay = false;
    return;
   end
end
if(iArg == 2)
   shiftLevel = varargin{1}{iArg};
   
end
if(iArg == 3)
   amp_ratio = varargin{1}{iArg};
   
end
if(iArg == 4)
   ws_r = varargin{1}{iArg};
   
end
if(iArg == 5)
   toFlip = varargin{1}{iArg};
   if(toFlip ~= 0 && toFlip ~= 1)
    fprintf('ERROR : toFlip must be 0 (Off) or 1 (On).\n');
    IsInputOkay = false;
    return;
   end
end
if(iArg == 6)
   numIteration = varargin{1}{iArg};
   if(numIteration < 1 || (mod(numIteration, 1) ~= 0))
    fprintf('ERROR : Number of Iteration must be an integer more than 0.\n');
    IsInputOkay = false;
    return;
   end
end
if(iArg == 7)
    typeSpline = varargin{1}{iArg};
    if(typeSpline ~= 1 && typeSpline ~= 2 && typeSpline ~= 3)
    fprintf('ERROR : typeSpline must be 1 (clamped spline); 2 (not a knot spline).\n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 8)
    toModifyBC = varargin{1}{iArg};
    if(toModifyBC ~= 0 && toModifyBC ~= 1 && toModifyBC ~= 2)
    fprintf('ERROR : toModifyBC must be 0 (None) ; 1 (modified linear extrapolation); 2 (Mirror Boundary)\n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 9)
    randType = varargin{1}{iArg};
    if(randType ~= 1 && randType ~= 2)
    fprintf('ERROR : randType must be 1 (uniformly distributed white noise) ; 2 (gaussian white noise).\n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 10)
    seedNo = varargin{1}{iArg};
    if(seedNo < 0 || seedNo >  2^32-1 || (mod(seedNo, 1) ~= 0))
    fprintf('ERROR : The value of seed must be an integer between 0 and 2^32 - 1. \n');
    IsInputOkay = false;
    return;
    end
end
if(iArg == 11)
    checkSignal = varargin{1}{iArg};
    if(checkSignal ~= 0 && checkSignal ~= 1)
    fprintf('ERROR : Number of checksignal must be 1 (Yes) or 0 (No).\n');
    IsInputOkay = false;
    return;
    end
end

end

end

if(NoiseLevel == 0)
    %fprintf('If NoiseLevel is ZERO, EEMD algorithm will be changed to EMD algorithm.\n');
end
if ((NE < 1) || (mod(NE, 1) ~= 0))
    fprintf('ERROR : Number of Ensemble must be integer more than 0.\n');
    IsInputOkay = false;
    return;
end


[m,n] = size(Y);
if(m ~= 1)
    if((n ~= 1))
       fprintf('ERROR : EMD could not input matrix array !\n');
       IsInputOkay = false;
       return;
    else
        Y =Y';
        xsize = m;
    end
else
    xsize = n;
end

if (checkSignal == 1)
    if((any(isinf(Y(:)) == 1)) || (any(isnan(Y(:)) == 1)))
        fprintf('ERROR : The input signal has NaN or Infinity elements.\n');
        IsInputOkay = false;
        return;
    end
end

if(mod(TNM, 1) ~= 0)
    fprintf('ERROR : TNM must be an integer more than 0. \n');
    IsInputOkay = false;
    return;
end
shiftLevel=fix(shiftLevel);

if (TNM <= 0) % automatic estimating number of imf 
    TNM=fix(log2(xsize));
end

end

function cnt=zero_cross_cnt(s)
n=length(s);
d1=zeros(size(s));
d2=zeros(size(s));
d1(1:n-1)=diff(s);
d2(2:n)=diff(s);
d1(n)=0;
d2(1)=0;
a=d1.*d2 < 0;
cnt=sum(a);
end
