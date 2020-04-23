function [cfg] = spmm_tfplot_r(cfg, varargin)

% selcomp=1:size(varargin{1}.topo,2);
% Fs=varargin{1}.Fs; % am loc
% fs=varargin{1}.fs; % fm loc
xmin=cfg.xmin;
xmax=cfg.xmax;
ymin=cfg.ymin;
ymax=cfg.ymax;
clim=cfg.clim;
fF=cfg.fF;
dyadic=cfg.dyadic;
if ~isfield(cfg,'op')
    op=0;
else
    op=cfg.op;
end
% if length(varargin)>1
    % coding P_masked here 
% else
    P_masked=zeros(size(varargin{1}.T));
%end
%chs=cfg.chs;
if ~isfield(cfg, 'collapse')
    collapse=0;
else
    collapse=cfg.collapse;
end
%if ~isfield(cfg, 'style')
    cfg.style=1;
%end
sptdata=varargin{1}.T;
% if isfield(varargin{1},'chlbls')
%     lbls=varargin{1}.chlbls;
% else
%     lbls=cellstr(string(1:size(sptdata,1)));
% end
clev=linspace (clim(1), clim(2),33);
%load colormap_myhot;
%load colormap_myjet;
%myhot=myhot(1:2:63,:);
%tint_jet=tint(jet,0.9);
tint_jet=jet;
if ~collapse

            %T=shiftdim(sptdata(chid,:,:),1);
            T=sptdata(:,:);
            figure
            %colormap(myjet)
            switch cfg.style
                case 1
                      %colormap(myjet);
                      colormap(tint_jet);
                      ham=imagesc([xmin, xmax],[ymin, ymax],T,clim);
                case 2
            end
            %title(lbls{chid});
            set(gca,'YDir','normal','FontSize',20,'FontName','Arial');
            set(gca,'box','on');
            if isfield(cfg,'title')
               title(cfg.title,'FontSize',20); 
            end
            if fF
                if dyadic
                    ytickmin=ceil(ymin/1)*1;
                    ytickmax=floor(ymax/1)*1;
                    if (ymax-ytickmax)>0.1
                        halfypix=0.5*(ymax-ymin)/(size(T,1)-1);
                        ymaxp=ymax+halfypix;
                        ytick_value=[ytickmin:1:ytickmax ymaxp];
                    else
                        ytick_value=ytickmin:1:ytickmax;
                    end
                    set(gca,'YTick',ytick_value);
                    ylblval=2.^ytick_value;
                    ylblval(end)=round(ylblval(end),1);
                    set(gca,'YTickLabel',ylblval);
                    tickmin=ceil(xmin/1)*1;
                    tickmax=floor(xmax/1)*1;
                    if (xmax-tickmax)>0.1
                        halfxpix=0.5*(xmax-xmin)/(size(T,2)-1);
                        xmaxp=xmax+halfxpix;
                        xtick_value=[tickmin:1:tickmax xmaxp];
                    else
                        xtick_value=tickmin:1:tickmax;
                    end
%                     set(gca,'XTick',xtick_value);
%                     set(gca,'XTickLabel',2.^xtick_value);
                    set(gca,'XTick',xtick_value);
                    xlblval=2.^xtick_value;
                    %xlblval(end)=round(xlblval(end),1);
                    xlblval(end)=round(xlblval(end)); % reset it to previous line after drafting
                    set(gca,'XTickLabel',xlblval);
                else
                    tickmin=ceil(xmin/10)*10;
                    tickmax=floor(xmax/10)*10;
                    set(gca,'XTick',tickmin:10:tickmax);
                    set(gca,'XTickLabel',tickmin:10:tickmax);
                end
                set(get(gca, 'XLabel'),'String','carrier frequency (Hz)','FontSize',20); 
                set(get(gca, 'YLabel'),'String','AM frequency (Hz)','FontSize',20);
                hold on
                ssy=linspace(xmin,ymax);
                %ssx=zeros(size(ssy));
                ssx=ssy;
                plot(ssx,ssy,'k:','LineWidth',2);

            else
                if dyadic
                    ytickmin=ceil(ymin/1)*1;
                    ytickmax=floor(ymax/1)*1;
                    if (ymax-ytickmax)>0.1
                        halfypix=0.5*(ymax-ymin)/(size(T,1)-1);
                        ymaxp=ymax+halfypix;
                        ytick_value=[ytickmin:1:ytickmax ymaxp];
                    else
                        ytick_value=ytickmin:1:ytickmax;
                    end
                   set(gca,'YTick',ytick_value);
                    ylblval=2.^ytick_value;
                    ylblval(end)=round(ylblval(end),1);
                    set(gca,'YTickLabel',ylblval);
                end
                tickmin=ceil(xmin/500)*500;
                tickmax=floor(xmax/500)*500;
                set(gca,'XTick',tickmin:500:tickmax);
                set(gca,'XTickLabel',tickmin:500:tickmax);
                set(get(gca, 'XLabel'),'String','time (ms)','FontSize',20); 
                set(get(gca, 'YLabel'),'String','carrier frequency (Hz)','FontSize',20);
                curr_xlim=get(gca,'xlim');
                hold on
                ssy=linspace(ymin,ymax);
                ssx=zeros(size(ssy));
                %ssx=ssy;
                plot(ssx,ssy,'k:','LineWidth',2);
                % delete after the 2012 tDCS WM project
                    ssx_m1100=-1100*ones(size(ssy));
                    plot(ssx_m1100,ssy,'k:','LineWidth',2);
                    ssx_m900=-900*ones(size(ssy));
                    plot(ssx_m900,ssy,'k:','LineWidth',2);
                % delete after the binding WM project
                    ssx_1000=1000*ones(size(ssy));
                    plot(ssx_1000,ssy,'k:','LineWidth',2);
                    ssx_3000=3000*ones(size(ssy));
                    plot(ssx_3000,ssy,'k:','LineWidth',2);
                set(gca,'xlim',curr_xlim);    
            end
            if length(varargin)>1
                xs=linspace(xmin,xmax,size(T,2));
                ys=linspace(ymin,ymax,size(T,1));
                cx=repmat(xs,size(T,1),1);
                cy=repmat(ys',1,size(T,2));

                maskch=P_masked(:,:);
                % separate + & -
                maskch_p=maskch & (T>0);
                maskch_n=maskch & (T<0);
                switch cfg.style
                case 1
                      if ~op
                          contour(cx,cy,maskch_p,1,'w','LineWidth',3);
                          contour(cx,cy,maskch_n,1,'w','LineWidth',3);
                      else
                          maskch=maskch+0;
                          maskch(maskch==0)=0.3;
                          set(ham,'alphadata',maskch);
                      end
                case 2
                     %contour(cx,cy,maskch,1,'r','LineWidth',2);
                end
            end
       old_pos=get(gcf,'position');
       new_pos=old_pos;
       new_pos(4)= old_pos(4)*1.1;
       new_pos(2)= old_pos(2)*0.8;
       set(gcf,'position',new_pos);
       set(gcf,'position',[1,1,1.4,1.15].*get(gcf,'position')); 
       cax_pos=get(gca,'position');
            
        amax_opos=get(ham.Parent,'position');
        amax_npos=amax_opos;
        amax_npos(3)=cax_pos(3);
        amax_npos(1)=cax_pos(1);
        set(ham.Parent,'position',amax_npos);
else

            T=sptdata(2:end,:);

            figure
            subplot(9,1,1:8);
            switch cfg.style
                case 1
                      %colormap(myjet);
                      colormap(tint_jet);
                      ham=imagesc([xmin, xmax],[ymin, ymax],T,clim);
                case 2
                      
            end
            %title(lbls{chid});
            set(gca,'YDir','normal','FontSize',20,'FontName','Arial');
            set(gca,'box','on');
            if isfield(cfg,'title')
               title(cfg.title,'FontSize',20); 
            end

            %colorbar
            if fF
                if dyadic
                    ytickmin=ceil(ymin/1)*1;
                    ytickmax=floor(ymax/1)*1;
                    if (ymax-ytickmax)>0.1
                        halfypix=0.5*(ymax-ymin)/(size(T,1)-1);
                        ymaxp=ymax+halfypix;
                        ytick_value=[ytickmin:1:ytickmax ymaxp];
                    else
                        ytick_value=ytickmin:1:ytickmax;
                    end
                    set(gca,'YTick',ytick_value);
                    ylblval=2.^ytick_value;
                    ylblval(end)=round(ylblval(end),1);
                    set(gca,'YTickLabel',ylblval);
                    tickmin=ceil(xmin/1)*1;
                    tickmax=floor(xmax/1)*1;
                    if (xmax-tickmax)>0.1
                        %xtick_value=[tickmin:1:tickmax xmax];
                        halfxpix=0.5*(xmax-xmin)/(size(T,2)-1);
                        xmaxp=xmax+halfxpix;
                        xtick_value=[tickmin:1:tickmax xmaxp];
                    else
                        xtick_value=tickmin:1:tickmax;
                    end
                    set(gca,'XTick',xtick_value);
                    %set(gca,'XTickLabel',2.^xtick_value);
                    set(gca,'XTickLabel',{});
                    %set(gca,'xlim',[xmin xmax]);
                else
                    tickmin=ceil(xmin/10)*10;
                    tickmax=floor(xmax/10)*10;
                    set(gca,'XTick',tickmin:10:tickmax);
                    %set(gca,'XTickLabel',tickmin:10:tickmax);
                    set(gca,'XTickLabel',{});
                    %set(gca,'xlim',[xmin xmax]);
                end
                %set(get(gca, 'XLabel'),'String','frequency (Hz)','FontSize',14); 
                set(get(gca, 'YLabel'),'String','AM frequency (Hz)','FontSize',20);
                hold on
                ssy=linspace(xmin,ymax);
                %ssx=zeros(size(ssy));
                ssx=ssy;
                plot(ssx,ssy,'k:','LineWidth',2);

            else
                if dyadic
                    ytickmin=ceil(ymin/1)*1;
                    ytickmax=floor(ymax/1)*1;
%                     if (ymax-ytickmax)>0.1
%                         ytick_value=[ytickmin:1:ytickmax ymax];
%                     else
%                         ytick_value=ytickmin:1:ytickmax;
%                     end
                    if (ymax-ytickmax)>0.1
                        halfypix=0.5*(ymax-ymin)/(size(T,1)-1);
                        ymaxp=ymax+halfypix;
                        ytick_value=[ytickmin:1:ytickmax ymaxp];
                    else
                        ytick_value=ytickmin:1:ytickmax;
                    end
                    set(gca,'YTick',ytick_value);
                    ylblval=2.^ytick_value;
                    ylblval(end)=round(ylblval(end),1);
                    set(gca,'YTickLabel',ylblval);
                end
                tickmin=ceil(xmin/500)*500;
                tickmax=floor(xmax/500)*500;
                set(gca,'XTick',tickmin:500:tickmax);
                %set(gca,'XTickLabel',tickmin:500:tickmax);
                set(gca,'XTickLabel',{});
                %set(gca,'xlim',[xmin xmax]);
                %set(get(gca, 'XLabel'),'String','Time (ms)','FontSize',14); 
                set(get(gca, 'YLabel'),'String','AM frequency (Hz)','FontSize',20);
                xlimv=get(gca,'xlim');
                hold on
                ssy=linspace(ymin,ymax);
                ssx=zeros(size(ssy));
                %ssx=ssy;
                plot(ssx,ssy,'k-','LineWidth',1);
                % delete after the 2012 tDCS WM project
                    ssx_m1100=-1100*ones(size(ssy));
                    plot(ssx_m1100,ssy,'k:','LineWidth',2);
                    ssx_m900=-900*ones(size(ssy));
                    plot(ssx_m900,ssy,'k:','LineWidth',2);
                % delete after the binding WM project
                    ssx_1000=1000*ones(size(ssy));
                    plot(ssx_1000,ssy,'k:','LineWidth',2);
                    ssx_3000=3000*ones(size(ssy));
                    plot(ssx_3000,ssy,'k:','LineWidth',2);
               set (gca,'xlim',xlimv);
            end
            if length(varargin)>1
                xs=linspace(xmin,xmax,size(T,2));
                ys=linspace(ymin,ymax,size(T,1));
                cx=repmat(xs,size(T,1),1);
                cy=repmat(ys',1,size(T,2));
                maskch=P_masked(2:end,:);
                % separate + & -
                maskch_p=maskch & (T>0);
                maskch_n=maskch & (T<0);
                switch cfg.style
                case 1
                    if ~ op
                        contour(cx,cy,maskch_p,1,'w','LineWidth',3);
                        contour(cx,cy,maskch_n,1,'w','LineWidth',3);
                    else
                        if (sigbool(4)==1 || sigbool(5)==1) || nonparam
                           maskch=maskch+0;
                           maskch(maskch==0)=0.3;
                           set(ham,'alphadata',maskch);
                        else
                            maskch=maskch+0;
                            maskch(:)=0.3;
                            set(ham,'alphadata',maskch);
                            contour(cx,cy,maskch_p,1,'r:','LineWidth',3);
                            contour(cx,cy,maskch_n,1,'r:','LineWidth',3);
                        end
                    end
                case 2
                    %contour(cx,cy,maskch,1,'r','LineWidth',2);    
                end

            end
            %axis tight
            T_dc=sptdata(1,:);
            T_dc3=repmat(T_dc,3,1);
            subplot(9,1,9);

            switch cfg.style
                case 1
                      hdc=imagesc([xmin, xmax],[-1, 1],T_dc3,clim);
                case 2
            end
            set(gca,'YDir','normal','FontSize',20,'FontName','Arial');
            set(gca,'box','on');
            axis tight
            %hb2=colorbar;

            if fF
                if dyadic
                    ytickmin=-1.5;
                    ytickmax=1.5;
                    ytick_value=0;
                    set(gca,'YTick',ytick_value);
                    set(gca,'YTickLabel','trend_A_M');
                    tickmin=ceil(xmin/1)*1;
                    tickmax=floor(xmax/1)*1;
                    if (xmax-tickmax)>0.1
                        halfxpix=0.5*(xmax-xmin)/(size(T_dc3,2)-1);
                        xmaxp=xmax+halfxpix;
                        xtick_value=[tickmin:1:tickmax xmaxp];
                    else
                        xtick_value=tickmin:1:tickmax;
                    end
                    set(gca,'XTick',xtick_value);
                    xlblval=2.^xtick_value;
                    %xlblval(end)=round(xlblval(end),1);
                    xlblval(end)=round(xlblval(end)); % reset it to previous line after drafting
                    set(gca,'XTickLabel',xlblval);
                    %set(gca,'xlim',[xmin xmax]);
                else
                    ytickmin=-1.5;
                    ytickmax=1.5;
                    ytick_value=0;
                    set(gca,'YTick',ytick_value);
                    set(gca,'YTickLabel','trend_A_M');
                    tickmin=ceil(xmin/10)*10;
                    tickmax=floor(xmax/10)*10;
                    set(gca,'XTick',tickmin:10:tickmax);
                    set(gca,'XTickLabel',tickmin:10:tickmax);
                    %set(gca,'xlim',[xmin xmax]);
                end
                set(get(gca, 'XLabel'),'String','carrier frequency (Hz)','FontSize',20); 
%                 set(get(gca, 'YLabel'),'String','AM (Hz)','FontSize',14);
%                 hold on
%                 ssy=linspace(xmin,ymax);
%                 %ssx=zeros(size(ssy));
%                 ssx=ssy;
%                 plot(ssx,ssy,'k:','LineWidth',2);

             else
                if dyadic
                    ytickmin=-1.5;
                    ytickmax=1.5;
                    ytick_value=0;
                    set(gca,'YTick',ytick_value);
                    set(gca,'YTickLabel','trend_A_M');
                end
                tickmin=ceil(xmin/500)*500;
                tickmax=floor(xmax/500)*500;
                set(gca,'XTick',tickmin:500:tickmax);
                set(gca,'XTickLabel',tickmin:500:tickmax);
                set(get(gca, 'XLabel'),'String','time (ms)','FontSize',20); 
            end
            %hold on

            set(gca,'ylim',[-0.9 0.9]);
           
       old_pos=get(gcf,'position');
       new_pos=old_pos;
       new_pos(4)= old_pos(4)*1.1;
       new_pos(2)= old_pos(2)*0.8;
       set(gcf,'position',new_pos);
       set(gcf,'position',[1,1,1.4,1.15].*get(gcf,'position')); 
       colorbar(ham.Parent);
       cax_pos=get(gca,'position');
            
            amax_opos=get(ham.Parent,'position');
            amax_npos=amax_opos;
            amax_npos(3)=cax_pos(3);
            amax_npos(1)=cax_pos(1);
            set(ham.Parent,'position',amax_npos);
end
