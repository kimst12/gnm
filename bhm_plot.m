function bhm_plot(infile,plt_opt,c5,mkr5,xoffset)
%------------------------------bhm_plot.m-------------------------------
% This function plots results from the branching hexagon model. The
% resulting plots differ depending on what inputs are given. For
% 'node_index' or 'node_type', the plots should be apparent. 
% For 'pressure', the following plots are generated:
%       figure 1 - a 3D plot where the node color indicates the log of pressure
%       figure 2 - scatter plot of pressure versus pial generation (y) where the marker shape/color indicates the node type
%       figure 3 - scatter plot of depth (z) versus pressure where the marker shape/color indicates the node type
%       figure 4 - error bar plot of depth (z) versus pressure where the marker shape/color indicates the node type
%       figure 5 - error bar plot of pressure versus node type, where the marker/color is specified as an input
%
% For 'volume_flow_rate', 'speed', 'Re', or 'Pe' the following plots are generated:
%       figure 1 - a 3D plot where the edge color indicates the log of the specified plt_opt
%       figure 2 - scatter plot of the specified plt_opt versus pial generation (y) where the marker shape/color indicates the edge type
%       figure 3 - scatter plot of depth (z) versus the specified plt_opt where the marker shape/color indicates the edge type
%       figure 4 - error bar plot of depth (z) versus the specified plt_opt where the marker shape/color indicates the edge type
%       figure 5 - error bar plot of the specified plt_opt versus edge type, where the marker/color is specified as an input
%
% Inputs: 
%   infile - the input file containing the results from a simulation
%   plt_opt - a string specifying which quantity/quantities to plot:
%               'node_index' - a plot with index numbers superimposed on
%                              top of every node
%               'node_type' - a plot with pial edges/nodes indicated in 
%                             black, penetrating edges/nodes indicated in 
%                             red, capillary edges/nodes indicated in
%                             green, and parenchymal edges/nodes indicated
%                             in purple
%               'pressure' - plots of the pressure
%               'volume_flow_rate' - plots of the volume flow rate
%               'speed' - plots of the speed
%               'Re' - plots of the Reynolds number
%               'Pe' - plots of the Peclet number
%     c5 - color of marker on fig 5
%     mkr5 - marker symbol on fig 5
%     xoffset - how much to shift the x axis by (value from 0 to 1) in fig 5 so markers from
%     different sets of results don't overlap. If you're only plotting one
%     set of results this number should be 0.
% Example usage: 
% bhm_plot(model_results_name,'volume_flow_rate','b','o',0) -- plots the volume flow rate in blue with circles as markers and no offset on the x-axis. 
% bhm_plot(model_results_name,'volume_flow_rate','r','x',0.5) -- plots the volume flow rate in red with x's as markers and an offset halfway to the next marker on the x-axis. 
%--------------------------------------------------------------------------
% Load infile
load(infile)

% Save figures?
savefigs=0;

% Set different types of markers
mkr_pial='o';
mkr_pen='s';
mkr_cap='^';
mkr_paren='p';
c_pial=[0 0 0]; % black
c_pen=[1 0 0]; % red
c_cap=[0 1 0]; % green
c_paren=lines(4);c_paren=c_paren(4,:); % purple

% choose viewing orientation
depthOnVert=1; %1 for depth on vertical axis; 0 for depth on horizontal axis 

% Plotting parameters
flm=0; % fill markers and color according to value (as opposed to vessel type)?
msz=60; % marker size
msz_plt=6;
yscl='log'; % 'linear' or 'log'

% Set percentiles for color limits on plots
prc=[0 100];

% Scale l_penart_plot
l_penart_plot=l_penart_plot*1e-6;

% Set rho and nu
nu=0.697e-6; % for water at 36.8 deg C (from Mestre et al 2018)
%D=6.55e-13; % for 1 micron particles (from Mestre et al 2018)
D=1.35e-10; % monomeric amyloid-beta (from Novo et al, Sci Rep 2018)
L_paren=20e-9; % approximate width of parenchymal channels

% Set viewing angle parameters
az=45;
el=30;

%% Schematic with node index
if strcmp(plt_opt,'node_index')
    f=figure;
    set(gcf,'position',[0 0 900 700]);
    hold on
    % Draw every edge
    for i=1:size(edges,1)
        if ~isnan(edges(i,1)) && ~isnan(edges(i,2)) && ~is_paren(edges(i,2))
            plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],'-k','linewidth',2);
        elseif ~isnan(edges(i,1)) && ~isnan(edges(i,2)) && is_paren(edges(i,2))
            plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],':k','linewidth',2);
        end
    end
    % Draw every node
    scatter3(xnode,ynode,znode,msz,'r','filled');
    % Add labels to the nodes
    for i=1:size(xnode,2)
        if znode(i)==0
            text(xnode(i)-0.05,ynode(i),znode(i),num2str(i),'fontsize',8,'fontweight','bold');
        end
    end
    % Format
    set(gca,'fontsize',18,'xtick',[],'ytick',[],'ztick',[]);
    daspect([1 1 1])
    box on
    view(az,el);
    
    % Save figure
    if savefigs
        export_fig([plt_opt '_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
    end

%% Schematic of network
elseif strcmp(plt_opt,'node_type')
    f=figure;
    set(gcf,'position',[0 0 900 700]);
    hold on
    
    % Draw every edge
    for i=1:size(edges,1)
        if ~isnan(edges(i,1)) && ~isnan(edges(i,2))
            if ~is_penart(edges(i,1)) & ~is_penart(edges(i,2)) & ~is_cap(edges(i,1)) & ~is_cap(edges(i,2)) & ~is_paren(edges(i,1)) & ~is_paren(edges(i,2))
                plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],'-','linewidth',4,'color',c_pial);
            elseif is_penart(edges(i,2))
                plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],'-','linewidth',3,'color',c_pen);
            elseif is_cap(edges(i,2))
                plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],'-','linewidth',2,'color',c_cap);
            elseif is_paren(edges(i,2))
                plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],':','linewidth',2,'color',c_paren);
            end
        end
    end    
    % Draw every node
    ind=find(~is_penart & ~is_cap & ~is_paren);
    scatter3(xnode(ind),ynode(ind),znode(ind),msz,c_pial,'filled','marker',mkr_pial);
    ind=find(is_penart);
    scatter3(xnode(ind),ynode(ind),znode(ind),msz,c_pen,'filled','marker',mkr_pen);
    ind=find(is_cap);
    scatter3(xnode(ind),ynode(ind),znode(ind),msz,c_cap,'filled','marker',mkr_cap);
    ind=find(is_paren);
    scatter3(xnode(ind),ynode(ind),znode(ind),msz,c_paren,'filled','marker',mkr_paren,'linewidth',2);
    
    % Format
    set(gca,'fontsize',18,'xtick',[],'ytick',[],'ztick',[]);
    daspect([1 1 1])
    box on
    view(az,el);
    set(gca,'clim',[-0.5 3.5]);
    
    % Save figure
    if savefigs
        export_fig([plt_opt '_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
    end

%% Pressure plots
elseif strcmp(plt_opt,'pressure')
    
    % Schematic plot
    f1=figure(1);
    set(gcf,'position',[0 0 900 700]);
    hold on
    
    % Loop and plot edges
    for i=1:size(edges,1)
        if ~isnan(edges(i,1)) && ~isnan(edges(i,2)) && ~is_paren(edges(i,2))
            plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],'-k','linewidth',2);
        elseif ~isnan(edges(i,1)) && ~isnan(edges(i,2)) && is_paren(edges(i,2))
            plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],':k','linewidth',2);
        end
    end
        
    % Draw every node
    ind=find(~is_penart & ~is_cap & ~is_paren);
    scatter3(xnode(ind),ynode(ind),znode(ind),msz,log10(pressures(ind)),'filled','marker',mkr_pial);
    ind=find(is_penart);
    scatter3(xnode(ind),ynode(ind),znode(ind),msz,log10(pressures(ind)),'filled','marker',mkr_pen);
    ind=find(is_cap);
    scatter3(xnode(ind),ynode(ind),znode(ind),msz,log10(pressures(ind)),'filled','marker',mkr_cap);
    ind=find(is_paren);
    scatter3(xnode(ind),ynode(ind),znode(ind),msz,log10(pressures(ind)),'filled','marker',mkr_paren);
    
    % Format
    set(gca,'fontsize',18,'xtick',[],'ytick',[],'ztick',[],'clim',[prctile(real(log10(pressures)),prc(1)) prctile(real(log10(pressures)),prc(2))]);
    cb=colorbar;
    colormap jet
    ylabel(cb,'Log_{10}[Pressure (mmHg)]','fontsize',20);
    daspect([1 1 1])
    box on
    view(az,el);

    % 2D plot (pial generation)
    f2=figure(2);
    set(gcf,'position',[0 0 900 700]);
    hold on
    % Get indices of penetrating arterioles and plot (either filled or open)
    if flm
        ind=find(~is_penart & ~is_cap & ~is_paren);
        scatter(ynode(ind),pressures(ind),msz,log10(pressures(ind)),'filled','marker',mkr_pial);
        ind=find(is_penart);
        scatter(ynode(ind),pressures(ind),msz,log10(pressures(ind)),'filled','marker',mkr_pen);
        ind=find(is_cap);
        scatter(ynode(ind),pressures(ind),msz,log10(pressures(ind)),'filled','marker',mkr_cap);
        ind=find(is_paren);
        scatter(ynode(ind),pressures(ind),msz,log10(pressures(ind)),'filled','marker',mkr_paren);
    else
        ind=find(~is_penart & ~is_cap & ~is_paren);
        scatter(ynode(ind),pressures(ind),msz,c_pial,'filled','marker',mkr_pial);
        ind=find(is_penart);
        scatter(ynode(ind),pressures(ind),msz,c_pen,'filled','marker',mkr_pen);
        ind=find(is_cap);
        scatter(ynode(ind),pressures(ind),msz,c_cap,'filled','marker',mkr_cap);
        ind=find(is_paren);
        scatter(ynode(ind),pressures(ind),msz,c_paren,'filled','marker',mkr_paren);
    end
    % Set colormap and format
    colormap jet
    set(gca,'yscale',yscl,'fontsize',18,'clim',[prctile(real(log10(pressures)),prc(1)) prctile(real(log10(pressures)),prc(2))]);%,'xtick',[],'ytick',[]);
    ylim([prctile(pressures,prc(1)) prctile(pressures,prc(2))])
    ylabel('Pressure (mmHg)','fontsize',20)
    xlabel('Pial Generation','fontsize',20)
    box on

    % 2D plot (depth below surface)
    f3=figure(3);
    set(gcf,'position',[0 0 900 700]);
    hold on
    % Get indices of penetrating arterioles and plot (either filled or open)
    if flm
        ind=find(~is_penart & ~is_cap & ~is_paren);
        scatter(abs(znode(ind))*l_penart/l_penart_plot,pressures(ind),msz,log10(pressures(ind)),'filled','marker',mkr_pial);
        ind=find(is_penart);
        scatter(abs(znode(ind))*l_penart/l_penart_plot,pressures(ind),msz,log10(pressures(ind)),'filled','marker',mkr_pen);
        ind=find(is_cap);
        scatter(abs(znode(ind))*l_penart/l_penart_plot,pressures(ind),msz,log10(pressures(ind)),'filled','marker',mkr_cap);
        ind=find(is_paren);
        scatter(abs(znode(ind))*l_penart/l_penart_plot,pressures(ind),msz,log10(pressures(ind)),'filled','marker',mkr_paren);
    else
        ind=find(~is_penart & ~is_cap & ~is_paren);
        scatter(abs(znode(ind))*l_penart/l_penart_plot,pressures(ind),msz,c_pial,'filled','marker',mkr_pial);
        ind=find(is_penart);
        scatter(abs(znode(ind))*l_penart/l_penart_plot,pressures(ind),msz,c_pen,'filled','marker',mkr_pen);
        ind=find(is_cap);
        scatter(abs(znode(ind))*l_penart/l_penart_plot,pressures(ind),msz,c_cap,'filled','marker',mkr_cap);
        ind=find(is_paren);
        scatter(abs(znode(ind))*l_penart/l_penart_plot,pressures(ind),msz,c_paren,'filled','marker',mkr_paren);
    end
    % Set colormap and format
    colormap jet
    set(gca,'yscale',yscl,'fontsize',18,'clim',[prctile(real(log10(pressures)),prc(1)) prctile(real(log10(pressures)),prc(2))]);%,'xtick',[],'ytick',[]);
    ylim([prctile(pressures,prc(1)) prctile(pressures,prc(2))])
    ylabel('Pressure (mmHg)','fontsize',20)
    xlabel('Depth below surface (\mum)','fontsize',20)
    box on
    if depthOnVert
        view([90 -90])
        set(gca, 'xdir', 'reverse')
    end
        
    % 2D plot (depth below surface w/ mean, err bars)
    f4=figure(4);
    set(gcf,'position',[0 0 900 700]);
    hold on
    % Get indices of penetrating arterioles and plot (either filled or open)
    ind=find(~is_penart & ~is_cap & ~is_paren);
    [x,y]=find_stats([abs(znode(ind))*l_penart/l_penart_plot;pressures(ind)]');
    errorbar(x,y(:,1),y(:,1)-y(:,3),y(:,2)-y(:,1),mkr_pial,'MarkerSize',msz_plt,'MarkerFaceColor',c_pial,'color',c_pial);
    ind=find(is_penart);
    [x,y]=find_stats([abs(znode(ind))*l_penart/l_penart_plot;pressures(ind)]');
    errorbar(x,y(:,1),y(:,1)-y(:,3),y(:,2)-y(:,1),mkr_pen,'MarkerSize',msz_plt,'MarkerFaceColor',c_pen,'color',c_pen);
    ind=find(is_cap);
    [x,y]=find_stats([abs(znode(ind))*l_penart/l_penart_plot;pressures(ind)]');
    errorbar(x,y(:,1),y(:,1)-y(:,3),y(:,2)-y(:,1),mkr_cap,'MarkerSize',msz_plt,'MarkerFaceColor',c_cap,'color',c_cap);
    ind=find(is_paren);
    [x,y]=find_stats([abs(znode(ind))*l_penart/l_penart_plot;pressures(ind)]');
    errorbar(x,y(:,1),y(:,1)-y(:,3),y(:,2)-y(:,1),mkr_paren,'MarkerSize',msz_plt,'MarkerFaceColor',c_paren,'color',c_paren);

    % Set colormap and format
    colormap jet
    set(gca,'yscale',yscl,'fontsize',18,'clim',[prctile(real(log10(pressures)),prc(1)) prctile(real(log10(pressures)),prc(2))]);%,'xtick',[],'ytick',[]);
    ylim([prctile(pressures,prc(1)) prctile(pressures,prc(2))])
    ylabel('Pressure (mmHg)','fontsize',20)
    xlabel('Depth below surface (\mum)','fontsize',20)
    box on
    if depthOnVert
        view([90 -90])
        set(gca, 'xdir', 'reverse')
    end
    
    % 2D plot (four points w/ mean, err bars)
    f5=figure(5);
    set(gcf,'position',[0 0 900 700]);
    hold on
    % Get indices of penetrating arterioles and plot (either filled or open)
    ind=find(~is_penart & ~is_cap & ~is_paren);
    y=[mean(pressures(ind)) max(pressures(ind)) min(pressures(ind))];
    errorbar(1+xoffset,y(1),y(1)-y(3),y(2)-y(1),mkr5,'MarkerSize',msz_plt,'MarkerFaceColor',c5,'color',c5,'linewidth',2);
    ind=find(is_penart);
    y=[mean(pressures(ind)) max(pressures(ind)) min(pressures(ind))];
    errorbar(2+xoffset,y(1),y(1)-y(3),y(2)-y(1),mkr5,'MarkerSize',msz_plt,'MarkerFaceColor',c5,'color',c5,'linewidth',2);
    ind=find(is_cap);
    y=[mean(pressures(ind)) max(pressures(ind)) min(pressures(ind))];
    errorbar(3+xoffset,y(1),y(1)-y(3),y(2)-y(1),mkr5,'MarkerSize',msz_plt,'MarkerFaceColor',c5,'color',c5,'linewidth',2);
    ind=find(is_paren);
    y=[mean(pressures(ind)) max(pressures(ind)) min(pressures(ind))];
    errorbar(4+xoffset,y(1),y(1)-y(3),y(2)-y(1),mkr5,'MarkerSize',msz_plt,'MarkerFaceColor',c5,'color',c5,'linewidth',2);

    % Set colormap and format
    colormap jet
    set(gca,'yscale',yscl,'fontsize',18,'xtick',[1:4],'xticklabel',['       Pial PVSs';'Penetrating PVSs';'  Capillary PVSs';'      Parenchyma']);
    xtickangle(45)
    xlim([0.5 4.5])
    ylim([prctile(pressures,prc(1)) prctile(pressures,prc(2))])
    ylabel('Pressure (mmHg)','fontsize',20)
    box on

    % Save figure
    if savefigs
        figure(f1)
        export_fig([plt_opt '_schematic_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
        figure(f2)
        export_fig([plt_opt '_2D_pialgen_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
        figure(f3)
        export_fig([plt_opt '_2D_depth_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
        figure(f4)
        export_fig([plt_opt '_2D_depth_avg_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
        figure(f5)
        export_fig([plt_opt '_2D_overall_avg_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);

    end

%% Set quantity to plot
else
    % Volume flow rate
    if strcmp(plt_opt,'volume_flow_rate')
        X=abs(Q')*1e3;

    % Speed (Q/A)
    elseif strcmp(plt_opt,'speed')
        cf=1.667e-8*1e6; % 1 mL/min = 1.667e-8 m^3/s, then 1e6 to convert to micron/s
        X=abs(Q')./PVS_area*cf; % This will return "nan" for parenchymal nodes
        ind=find(is_paren(edges(:,2)));
        X(ind)=abs(Q(ind))./(d_penart/2*sqrt(K(2)+1)*l_penart/ncap_per_penart)*cf; % This is volume flow rate divided by surface area of outer PVS wall

    % Reynolds number (U*L/nu)
    elseif strcmp(plt_opt,'Re')
        cf=1.667e-8; % 1 mL/min = 1.667e-8 m^3/s
        L=d_edge/2.*sqrt(K(vessel_type)+1)-d_edge/2; % L=r2-r1 for circular annulus; returns nan for parenchyma
        ind=find(is_paren(edges(:,2))); % get parenchyma node indices
        L(ind)=L_paren; % set length scale for parenchyma
        U=abs(Q')./PVS_area*cf; % velocity; returns nan for parenchyma
        U(ind)=abs(Q(ind))./(d_penart/2*sqrt(K(2)+1)*l_penart/ncap_per_penart)*cf; % parenchymal velocity; this is volume flow rate divided by surface area of outer PVS wall
        X=U.*L/nu;

    % Peclet number (U*L/D)
    elseif strcmp(plt_opt,'Pe')
        cf=1.667e-8; % 1 mL/min = 1.667e-8 m^3/s
        L=d_edge/2.*sqrt(K(vessel_type)+1)-d_edge/2; % L=r2-r1 for circular annulus; returns nan for parenchyma
        ind=find(is_paren(edges(:,2))); % get parenchyma node indices
        L(ind)=L_paren; % set length scale for parenchyma
        U=abs(Q')./PVS_area*cf; % velocity; returns nan for parenchyma
        U(ind)=abs(Q(ind))./(d_penart/2*sqrt(K(2)+1)*l_penart/ncap_per_penart)*cf; % parenchymal velocity; this is volume flow rate divided by surface area of outer PVS wall
        X=U.*L/D;
    end

    %% Generate plots

    % Compute axes range for plotting
    Xrange=[log10(prctile(X,prc(1))) log10(prctile(X,prc(2)))];
    
    % Schematic
    f1=figure(1);
    set(gcf,'position',[0 0 900 700]);
    hold on
    % Set up color map and bins to use for plotting lines
    Xbin=linspace(Xrange(1),Xrange(2),101);
    cl=jet(100);
    % Loop and plot colored lines according to value of X
    for i=1:size(edges,1)
        if ~isnan(X(i))
            logX=log10(X(i));
            ind=find(logX>=Xbin(1:end-1) & logX<=Xbin(2:end));
            if ~isempty(ind) && ~is_paren(edges(i,2)) && ~is_cap(edges(i,2)) && ~is_penart(edges(i,2))
                plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],'-','linewidth',4,'color',cl(ind,:));
            elseif ~isempty(ind) && is_penart(edges(i,2))
                plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],'-','linewidth',3,'color',cl(ind,:));
            elseif ~isempty(ind) && is_paren(edges(i,2))
                plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],':','linewidth',2,'color',cl(ind,:));
            elseif ~isempty(ind) && ~is_paren(edges(i,2))
                plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],'-','linewidth',2,'color',cl(1,:));
            elseif isempty(ind) && is_paren(edges(i,2))
                plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],':','linewidth',2,'color',cl(1,:));
            elseif isempty(ind) && ~is_paren(edges(i,2))
                plot3([xnode(edges(i,1)) xnode(edges(i,2))],[ynode(edges(i,1)) ynode(edges(i,2))],[znode(edges(i,1)) znode(edges(i,2))],'-','linewidth',2,'color',cl(1,:));
            end
        end
    end

    % Add colorbar and label it
    set(gca,'fontsize',18,'clim',Xrange);
    cb=colorbar;
    % Set colormap and format
    colormap jet
    set(gca,'fontsize',18,'xtick',[],'ytick',[],'ztick',[]);
    box on
    daspect([1 1 1])
    % Add proper label
    if strcmp(plt_opt,'volume_flow_rate');ylabel(cb,'Log_{10}[Volume flow rate (\muL/min)]','fontsize',20);
    elseif strcmp(plt_opt,'speed');ylabel(cb,'Log_{10}[Average flow speed (\mum/s)]','fontsize',20);
    elseif strcmp(plt_opt,'Re');ylabel(cb,'Log_{10}[Re]','fontsize',20);
    elseif strcmp(plt_opt,'Pe');ylabel(cb,'Log_{10}[Pe]','fontsize',20);
    end
    view(az,el);

    % 2D plot (pial generation)
    f2=figure(2);
    set(gcf,'position',[0 0 900 700]);
    hold on
    % Get indices of surface penetrating arterioles and plot volume flow rate
    if flm
        ind=find(~is_penart(edges(:,1)) & ~is_penart(edges(:,2)) & ~is_cap(edges(:,1)) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
        mh(1)=scatter(ynode(edges(ind,1)),X(ind),msz,log10(X(ind)),'filled','marker',mkr_pial);
        ind=find((is_penart(edges(:,1)) | is_penart(edges(:,2))) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
        mh(2)=scatter(ynode(edges(ind,1)),X(ind),msz,log10(X(ind)),'filled','marker',mkr_pen);
        ind=find(is_cap(edges(:,2)));
        mh(3)=scatter(ynode(edges(ind,1)),X(ind),msz,log10(X(ind)),'filled','marker',mkr_cap);
        ind=find(is_paren(edges(:,2)));
        mh(4)=scatter(ynode(edges(ind,1)),X(ind),msz,log10(X(ind)),'filled','marker',mkr_paren);
    else
        ind=find(~is_penart(edges(:,1)) & ~is_penart(edges(:,2)) & ~is_cap(edges(:,1)) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
        mh(1)=scatter(ynode(edges(ind,1)),X(ind),msz,c_pial,'filled','marker',mkr_pial);
        ind=find((is_penart(edges(:,1)) | is_penart(edges(:,2))) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
        mh(2)=scatter(ynode(edges(ind,1)),X(ind),msz,c_pen,'filled','marker',mkr_pen);
        ind=find(is_cap(edges(:,2)));
        mh(3)=scatter(ynode(edges(ind,1)),X(ind),msz,c_cap,'filled','marker',mkr_cap);
        ind=find(is_paren(edges(:,2)));
        mh(4)=scatter(ynode(edges(ind,1)),X(ind),msz,c_paren,'filled','marker',mkr_paren);
    end
    % Set colormap and format
    colormap jet
    set(gca,'yscale',yscl,'fontsize',18,'clim',Xrange);
    ylim(10.^Xrange)
    lh=legend(mh,'Pial PVSs','Penetrating PVSs','Capillary PVSs','Parenchyma');
    set(lh,'fontsize',18,'location','northeast');
    % Add proper label
    if strcmp(plt_opt,'volume_flow_rate');ylabel('Volume flow rate (\muL/min)','fontsize',20);
    elseif strcmp(plt_opt,'speed');ylabel('Average flow speed (\mum/s)','fontsize',20);
    elseif strcmp(plt_opt,'Re');ylabel('Re','fontsize',20);
    elseif strcmp(plt_opt,'Pe');ylabel('Pe','fontsize',20);
    end
    xlabel('Pial Generation','fontsize',20)
    box on

    % 2D plot (depth below surface)
    f3=figure(3);
    set(gcf,'position',[0 0 900 700]);
    hold on
    % Get indices of penetrating arterioles and plot (either filled or open)
    if flm
        ind=find(~is_penart(edges(:,1)) & ~is_penart(edges(:,2)) & ~is_cap(edges(:,1)) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
        mh(1)=scatter(abs(znode(edges(ind,1)))*l_penart/l_penart_plot,X(ind),msz,log10(X(ind)),'filled','marker',mkr_pial);
        ind=find((is_penart(edges(:,1)) | is_penart(edges(:,2))) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
        mh(2)=scatter(abs(znode(edges(ind,1)))*l_penart/l_penart_plot,X(ind),msz,log10(X(ind)),'filled','marker',mkr_pen);
        ind=find(is_cap(edges(:,2)));
        mh(3)=scatter(abs(znode(edges(ind,1)))*l_penart/l_penart_plot,X(ind),msz,log10(X(ind)),'filled','marker',mkr_cap);
        ind=find(is_paren(edges(:,2)));
        mh(4)=scatter(abs(znode(edges(ind,1)))*l_penart/l_penart_plot,X(ind),msz,log10(X(ind)),'filled','marker',mkr_paren);
    else
        ind=find(~is_penart(edges(:,1)) & ~is_penart(edges(:,2)) & ~is_cap(edges(:,1)) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
        mh(1)=scatter(abs(znode(edges(ind,1)))*l_penart/l_penart_plot,X(ind),msz,c_pial,'filled','marker',mkr_pial);
        ind=find((is_penart(edges(:,1)) | is_penart(edges(:,2))) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
        mh(2)=scatter(abs(znode(edges(ind,1)))*l_penart/l_penart_plot,X(ind),msz,c_pen,'filled','marker',mkr_pen);
        ind=find(is_cap(edges(:,2)));
        mh(3)=scatter(abs(znode(edges(ind,1)))*l_penart/l_penart_plot,X(ind),msz,c_cap,'filled','marker',mkr_cap);
        ind=find(is_paren(edges(:,2)));
        mh(4)=scatter(abs(znode(edges(ind,1)))*l_penart/l_penart_plot,X(ind),msz,c_paren,'filled','marker',mkr_paren);
    end
    % Set colormap and format
    colormap jet
    set(gca,'yscale',yscl,'fontsize',18,'clim',Xrange);
    ylim(10.^Xrange)
    lh=legend(mh,'Pial PVSs','Penetrating PVSs','Capillary PVSs','Parenchyma');
    set(lh,'fontsize',18,'location','northeast');
    % Add proper label
    if strcmp(plt_opt,'volume_flow_rate');ylabel('Volume flow rate (\muL/min)','fontsize',20);
    elseif strcmp(plt_opt,'speed');ylabel('Average flow speed (\mum/s)','fontsize',20);
    elseif strcmp(plt_opt,'Re');ylabel('Re','fontsize',20);
    elseif strcmp(plt_opt,'Pe');ylabel('Pe','fontsize',20);
    end
    xlabel('Depth below surface (\mum)','fontsize',20)
    box on
    if depthOnVert
        view([90 -90])
        set(gca, 'xdir', 'reverse')
    end
    
    % 2D plot (depth below surface w/ mean + error bars)
    f4=figure(4);
    set(gcf,'position',[0 0 900 700]);
    hold on
    % Get indices of penetrating arterioles and plot (either filled or open)
    ind=find(~is_penart(edges(:,1)) & ~is_penart(edges(:,2)) & ~is_cap(edges(:,1)) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,1)) & ~is_paren(edges(:,2)));
    [x,y]=find_stats([abs(znode(edges(ind,1)))*l_penart/l_penart_plot; X(ind)]');
    mh(1)=errorbar(x,y(:,1),y(:,1)-y(:,3),y(:,2)-y(:,1),[mkr_pial],'MarkerSize',msz_plt,'MarkerFaceColor',c_pial,'color',c_pial);
    ind=find((is_penart(edges(:,1)) | is_penart(edges(:,2))) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
    [x,y]=find_stats([abs(znode(edges(ind,1)))*l_penart/l_penart_plot; X(ind)]');
    mh(2)=errorbar(x,y(:,1),y(:,1)-y(:,3),y(:,2)-y(:,1),[mkr_pen],'MarkerSize',msz_plt,'MarkerFaceColor',c_pen,'color',c_pen);
    ind=find(is_cap(edges(:,2)));
    [x,y]=find_stats([abs(znode(edges(ind,1)))*l_penart/l_penart_plot; X(ind)]');
    mh(3)=errorbar(x,y(:,1),y(:,1)-y(:,3),y(:,2)-y(:,1),[mkr_cap],'MarkerSize',msz_plt,'MarkerFaceColor',c_cap,'color',c_cap);
    ind=find(is_paren(edges(:,2)));
    [x,y]=find_stats([abs(znode(edges(ind,1)))*l_penart/l_penart_plot; X(ind)]');
    mh(4)=errorbar(x,y(:,1),y(:,1)-y(:,3),y(:,2)-y(:,1),[mkr_paren],'MarkerSize',msz_plt,'MarkerFaceColor',c_paren,'color',c_paren);
    
    % Set colormap and format
    colormap jet
    set(gca,'yscale',yscl,'fontsize',18,'clim',Xrange);
    ylim(10.^Xrange)
    lh=legend(mh,'Pial PVSs','Penetrating PVSs','Capillary PVSs','Parenchyma');
    set(lh,'fontsize',18,'location','northeast');
    % Add proper label
    if strcmp(plt_opt,'volume_flow_rate');ylabel('Volume flow rate (\muL/min)','fontsize',20);
    elseif strcmp(plt_opt,'speed');ylabel('Average flow speed (\mum/s)','fontsize',20);
    elseif strcmp(plt_opt,'Re');ylabel('Re','fontsize',20);
    elseif strcmp(plt_opt,'Pe');ylabel('Pe','fontsize',20);
    end
    xlabel('Depth below surface (\mum)','fontsize',20)
    box on
    if depthOnVert
        view([90 -90])
        set(gca, 'xdir', 'reverse')
    end
    
    % 2D plot (four points w/ mean + error bars)
    f5=figure(5);
    set(gcf,'position',[0 0 900 700]);
    hold on
    % Get indices of penetrating arterioles and plot (either filled or open)
    ind=find(~is_penart(edges(:,1)) & ~is_penart(edges(:,2)) & ~is_cap(edges(:,1)) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,1)) & ~is_paren(edges(:,2)));
    y=[mean(X(ind)) max(X(ind)) min(X(ind))];
    mh(1)=errorbar(1+xoffset,y(1),y(1)-y(3),y(2)-y(1),'MarkerSize',msz_plt,'MarkerFaceColor',c5,'color',c5,'marker',mkr5,'linewidth',2);
    ind=find((is_penart(edges(:,1)) | is_penart(edges(:,2))) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
    y=[mean(X(ind)) max(X(ind)) min(X(ind))];
    mh(2)=errorbar(2+xoffset,y(1),y(1)-y(3),y(2)-y(1),'MarkerSize',msz_plt,'MarkerFaceColor',c5,'color',c5,'marker',mkr5,'linewidth',2);
    ind=find(is_cap(edges(:,2)));
    y=[mean(X(ind)) max(X(ind)) min(X(ind))];
    mh(3)=errorbar(3+xoffset,y(1),y(1)-y(3),y(2)-y(1),'MarkerSize',msz_plt,'MarkerFaceColor',c5,'color',c5,'marker',mkr5,'linewidth',2);
    ind=find(is_paren(edges(:,2)));
    y=[mean(X(ind)) max(X(ind)) min(X(ind))];
    mh(4)=errorbar(4+xoffset,y(1),y(1)-y(3),y(2)-y(1),'MarkerSize',msz_plt,'MarkerFaceColor',c5,'color',c5,'marker',mkr5,'linewidth',2);
    
    % Set colormap and format
    colormap jet
    set(gca,'yscale',yscl,'fontsize',18,'xtick',[1:4],'xticklabel',['       Pial PVSs';'Penetrating PVSs';'  Capillary PVSs';'      Parenchyma']);
    xtickangle(45)
    ylim(10.^Xrange)
    xlim([0.5 4.5])
    % Add proper label
    if strcmp(plt_opt,'volume_flow_rate');ylabel('Volume flow rate (\muL/min)','fontsize',20);
    elseif strcmp(plt_opt,'speed');ylabel('Average flow speed (\mum/s)','fontsize',20);
    elseif strcmp(plt_opt,'Re');ylabel('Re','fontsize',20);
    elseif strcmp(plt_opt,'Pe');ylabel('Pe','fontsize',20);
    end
    box on
    
    % Save figure
    if savefigs
        figure(f1)
        export_fig([plt_opt '_schematic_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
        figure(f2)
        export_fig([plt_opt '_2D_pialgen_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
        figure(f3)
        export_fig([plt_opt '_2D_depth_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
        figure(f4)
        export_fig([plt_opt '_2D_depth_avg_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
        figure(f5)
        export_fig([plt_opt '_2D_overall_avg_pial' num2str(dpdx_pial) '_pen' num2str(dpdx_pen) '_cap' num2str(dpdx_cap) '.pdf'],'-transparent');pause(1);
    end
end


end
