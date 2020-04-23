function [C]=cmask_emdn(S, TNM, noise ,varargin)
%---------------------------------------------------------------
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
[S,odrMask,TNM,numIteration,typeSpline,toModifyBC] = parse_checkProperty(S,  TNM, varargin);
%allmode = rcada_emd(S, toModifyBC, typeSpline, 2, numIteration);
%cnt=zero_cross_cnt(allmode(:,1));

n=size(S);
 
t=(1:length(S))';
if n(1) < n(2)
    S=S';
end
cnt=length(S)-2; % this n version start from Niquist
%cnt=0.64*(length(S)-2); % this n version start from 0.32*sampling rate
% if cnt<2
%     C=S;
%     return
% end
%SD0=std(S)*0.03;
sd_S=std(S);
w=pi*cnt/(length(S)-1); 
%w=pi;
if sd_S~=0
   Amp=noise*sd_S;
else
    Amp=noise;
    %disp('There is null signals');
end
%Amp=2*Amp;
% TNMc=-1;
% if TNMc < 0
    %TNMc=floor(log2(cnt)-1);
    TNMc=floor(log2(cnt)); % 20160429
% end

C=zeros(length(S),TNM);
if cnt<2
    C(:,1)=S;
    return
end
if TNMc>(TNM-1)
    TNMc=(TNM-1);
end
for i=1:TNMc
    denominator=2^odrMask;

    for j = 0:(denominator-1)
        %fprintf('Initial phase is %f \n',j/(2*denominator));
        MS1=Amp*cos(w*t+j*pi/(denominator*2));
        MS2=Amp*sin(w*t+j*pi/(denominator*2));
%         nc=TNMc-i+1;
%         if nc>2
           nc=2; 
%         end
        %Ctmp=emdx((S+MS1)',nc);
        Ctmp=rcada_emd((S+MS1)', toModifyBC, typeSpline, nc, numIteration);
        C(:,i)=C(:,i)+Ctmp(:,1);
        %Ctmp=emdx((S-MS1)',nc);
        Ctmp=rcada_emd((S-MS1)', toModifyBC, typeSpline, nc, numIteration);
        C(:,i)=C(:,i)+Ctmp(:,1);
        %Ctmp=emdx((S+MS2)',nc);
        Ctmp=rcada_emd((S+MS2)', toModifyBC, typeSpline, nc, numIteration);
        C(:,i)=C(:,i)+Ctmp(:,1);
        %Ctmp=emdx((S-MS2)',nc);
        Ctmp=rcada_emd((S-MS2)', toModifyBC, typeSpline, nc, numIteration);
        C(:,i)=C(:,i)+Ctmp(:,1);
    end
    %clear Ctmp;
    w=w/2;
    Amp=Amp/1.414;
    C(:,i)=C(:,i)/(4*denominator);
    S=S-C(:,i);
end
% if (TNMc+1)>=1
   C(:,TNMc+1)=S;
% else
%     C(:,1)=S;
% end

end

function [Y, odrMask, TNM,  numIteration, typeSpline,toModifyBC] = parse_checkProperty(Y, TNM, varargin)
% Default Parameters
odrMask = 0;
% shiftLevel=0;
% amp_ratio=2;
numIteration = 10; % numIteration = 10
typeSpline = 2;
toModifyBC = 1;


if(~isempty(varargin{1}))
    for iArg = 1 : length(varargin{1})
        if(iArg == 1)
           odrMask = varargin{1}{iArg};
           if odrMask < 0 && ~isinteger(odrMask)
            fprintf('ERROR : order of Masking signals must be 0 (4) or positive integer n (2^(n+2)).\n');
            return;
           end
        end
        if(iArg == 2)
           numIteration = varargin{1}{iArg};
           if(numIteration < 1 || (mod(numIteration, 1) ~= 0))
            fprintf('ERROR : Number of Iteration must be an integer more than 0.\n');
            return;
           end
        end
        if(iArg == 3)
            typeSpline = varargin{1}{iArg};
            if(typeSpline ~= 1 && typeSpline ~= 2 && typeSpline ~= 3)
            fprintf('ERROR : typeSpline must be 1 (clamped spline); 2 (not a knot spline).\n');
            return;
            end
        end
        if(iArg == 4)
            toModifyBC = varargin{1}{iArg};
            if(toModifyBC ~= 0 && toModifyBC ~= 1 && toModifyBC ~= 2)
            fprintf('ERROR : toModifyBC must be 0 (None) ; 1 (modified linear extrapolation); 2 (Mirror Boundary)\n');
            return;
            end
        end
    end
end
[m,n] = size(Y);
if(m ~= 1)
    if((n ~= 1))
       fprintf('ERROR : EMD could not input matrix array !\n');
       return;
    else
        Y =Y';
        xsize = m;
    end
else
    xsize = n;
end


if(mod(TNM, 1) ~= 0)
    fprintf('ERROR : TNM must be an integer more than 0. \n');
    return;
end
%shiftLevel=fix(shiftLevel);

if (TNM <= 0) % automatic estimating number of imf 
    TNM=fix(log2(xsize));
end

end


