load hhsa_ex1_data.mat
fs=1000;  % specify the sampling rate
TNM=-1; TNM2=-1; % determine the maximum number of IMFs automatically 
S.ifmethod='quad';  % use ¡§Direct Quadrature¡¨ to estimate 1st layer IF.
S.ifmethod2='quad';  % use ¡§Direct Quadrature¡¨ to estimate 2nd layer IF.
S.NEnsemble=0;  % set the order of Masking EMD for 1st layer EMD
S.NEnsemble2=0;  % set the order of Masking EMD for 2nd layer EMD
S.ENoise=2;  % set the masking amplitude factor for 1st layer EMD
S.ENoise2=2;  % set the masking amplitude factor for 2nd layer EMD
S.shiftLevel=1; % set the upsampling level, the new sampling rate is: 
% (2^ S.shiftLevel)*(original sampling rate)
%% Here, perform the 2-layer EMD
[fm, am, FM, AM, IMF, IMF2] = multi_EMD_DCM_SV(data,fs,TNM,TNM2,S);
%% Here, project the result of the 2-layer EMD to a 3D spectral space: (carrier frequency, AM frequency, time)
S.dyadic=1; % frequency represented in log2 scale
S.dyad_btw=8; % 8 bins within a log2 frequency scale. e.g. between 2^3 and 2^4 hz, there are 8 frequency bins 
S.collapse=1; 
ntp0=501; ntp1=2500; tres=500; % between 501ms and 2500ms, the time resolution is 500
fw0=0; fw1=7; fres=[];  % carrier frequency is from 2^0 to 2^7. fres is empty, because we already have S.dyad_btw=8
Fw0=-1; Fw1=6; Fres=[]; % AM frequency is from 2^-1 to 2^6
%[All_nt, fscale, Fscale]= nspplotf3d_tres2x( fm(:,1:end-1), FM(:,1:end,1:end-1), AM(:,1:end,1:end-1), ntp0, ntp1, fres, Fres, fw0, fw1, Fw0, Fw1, tres, S);
[All_nt, fscale, Fscale]= nspplotf3d_tres3x( fm(:,1:end-1), FM(:,1:end,1:end-1), AM(:,1:end,1:end-1), ntp0, ntp1, fres, Fres, fw0, fw1, Fw0, Fw1, tres, S);
All_nt3=sum( All_nt, 4); % obtain the 3D full HHS by summing over all IMFs
All_nt2=sum( All_nt3, 3); % get the 2D (AM frequency ¡Ñ frequency) HHS by summing over all time samples
% Here we have finished the HHSA procedure.
% In the following, we simply use a contour plot to represent the HHS.
% However, there should be a lot of alternatives for ploting the result.
% q1=fspecial('gaussian', 3, 0.6); % if MATLAB: image processing toolbox is available
% q2=fspecial('gaussian', 3, 0.6);
% if MATLAB: image processing toolbox is not available
q1=[0.0276818087794658,0.111014893010991,0.0276818087794658;0.111014893010991,0.445213192838173,0.111014893010991;0.0276818087794658,0.111014893010991,0.0276818087794658];
q2=q1;

ps_q1=filter2(q1, All_nt2'); % a gentle smoothing procedure for the spectrum
ps_q2=filter2(q2, ps_q1);
ymin=-1;
ymax=6;
xmin=0;
xmax=7;
cfg.xmin=xmin; cfg.xmax=xmax; cfg.ymin=ymin; cfg.ymax=ymax;
cfg.clim=[-8,4];
cfg.fF=1; cfg.dyadic=S.dyadic; cfg.collapse=S.collapse;
spect.T=log(ps_q2);
spmm_tfplot_r(cfg, spect);

% figure; contour( fscale, Fscale, log(ps_q2(2:end,:)), 30); 
% axis xy;
% xlabel('carrier frequency(Hz)'); ylabel('AM frequency(Hz)');% zlabel('Energy of envelope');
% ytickmin=ceil(ymin/1)*1;
% ytickmax=floor(ymax/1)*1;
% ytick_value=ytickmin:1:ytickmax;
% set(gca,'YTick',ytick_value);
% set(gca,'YTickLabel',2.^ytick_value);
% set(gca,'fontsize',14);
% tickmin=ceil(xmin/1)*1;
% tickmax=floor(xmax/1)*1;
% xtick_value=tickmin:1:tickmax;
% set(gca,'XTick',xtick_value);
% set(gca,'XTickLabel',2.^xtick_value);
% set(gca,'clim',[-8 2]);
% hold on
% ssy=linspace(xmin,ymax);
% ssx=ssy;
% plot(ssx,ssy,'k:','LineWidth',2);
% colorbar;