function [Qtotal,Rtotal]=branching_hexagon_model_pext(savename,p0,param)
%------------------branching_hexagon_model.m-------------------------------
% This code generates a hexagon model of perivascular spaces with the
% assumption that the flow is driven by a steady pressure gradient. There
% is a fixed influx Q from the base and then all penetrating arterioles are
% each connected to a separate large resistor, representing the parenchyma,
% which then connect to a common point which has efflux Q. 
% (This description should be updated after the geometry is finalized.)
%
% Description of input parameters:
% savename: name of the file to which the entire workspace will be saved
% p0: driving pressure difference between the inlet and ground. Units: mmHg
% param.dpdx: 1x4 vector with the pressure gradient caused by perivascular
%   pumping at the pial, penetrating, capillary, and parenchymal level. 
%   Note that this is formulated as a pressure per unit distance for the
%   the first three, but a pressure difference for the parenchyma. 
%   Units: mmHg/m for pial, penetrating, capillary; mmHg for parenchyma
% param.sc: short circuit the model by connecting the distal pial nodes to ground?
% param.C_paren: hydraulic conductivity (1/R) for flow from the penetrating nodes 
%   to the perivenous nodes, through the parenchyma. Units: mL/(min*mmHg)
% param.C_efflux: hydraulic conductivity (1/R) for flow from the start of the 
%   perivenous side to ground (beyond parenchymal or capillary flow).
%   Units: mL/(min*mmHg)
% param.g: number of generations to include in the model.
% param.K: a 1X4 vector containing the PVS-to-artery area ratio at the 
%   pial, penetrating, capillary, and parenchymal levels, respectively. The
%   last entry should be nan since this isn't defined for the parenchyma.
% param.d_pialart: diameter of pial arteries. Units: m 
% param.l_pial_art: length of pial arteries. Units: m
% param.l_pial2pen: length of pial to penetrating arteries (the horizontal 
%   offshoots from the main hexagon). Units: m
% param.d_penart: diameter of penetrating artery. Units: m
% param.l_penart: total length of penetrating artery. Units: m
% param.r_cap: radius of capillaries. Units: m
% param.l_cap: length of the capillaries. Units: m
% param.ncap_per_penart: number of capillaries branching off of each 
%   penetrating artery
% param.kappa: penetrating and capillary permeability; set to nan to model
%   the PVS as an open space. Units: m^2
%
% Example: 
%   param.g=6;
%   param.K=[1.4 1 1 nan];
%   param.d_pialart=50e-6;
%   param.l_pial_art=3*175e-6;
%   param.l_pial2pen=175e-6;
%   param.d_penart=11e-6;
%   param.l_penart=1e-3;
%   param.r_cap=3e-6;
%   param.l_cap=50e-6;
%   param.ncap_per_penart=11;
%   param.kappa=[nan nan 1.8e-14];
%   param.dpdx=[0;0;0;0];
%   param.sc=0;
%   param.C_paren=2.7379e-09;
%   param.C_efflux=1e3;
%   [Qtotal,Rtotal]=branching_hexagon_model('test_sim',0.4415,param)
%
% Written by Jeff Tithof on 25 Sept 2019
% 9 Oct 2019 - made substantial edits to make the code more readable and
%   removed commented portions that are no longer needed; in future
%   improvements, the units should be changed to avoid inaccuracies 
%   associated with very large / very small numbers
% 17 December 2019 - made updates to plot everything in 3D
% 10 April 2020 - updated to apply a constant external pressure gradient to
%   drive the flow
% 24 April 2020 - started making edits to add a "battery" in series with 
%   each individual resistor; this version is based on editing the C matrix
%   after it's been created; completed and working end of June 2020
% 1 July 2020 - added an option for "short-circuiting" the network by
%   connecting the last pial nodes to the same efflux as the capillaries
% 2 July 2020 - added an option for parenchymal flow by connecting every
%   node to the grounded efflux node with a high resistance; this
%   resistance is determined by the local area of PVS outer wall, an
%   average artery-to-venule distance of 250 micron, and a value from Jin
%   et al 2016
% 1 Oct 2020 - (KASB) changes to parameters 
%   1.) change r_penart to 5.5 microns (from 12.5)
%   2.) changed l_pial_art to be 3*175 so that each segment has a length 
%       that long
%   3.) created a new parameter, l_pial2pen, to describe the length of the 
%       horizonal transition edge, and replaced all the previous calls to 
%       get_pial_length with that parameter, which is currently set to 175 
%       microns, but could easily be adjusted.
% 12 Oct 2020 - (KASB) feed all variables into the function as parameters; 
%   Replaced the function get_pial_length with a variable
%   Changed K to be 1 X 3 vector, corresponding to K_pial, K_pen, or K_cap,
%     dependent on vessel type.
%   Change dpdx to be a 1 x 3 vector, corresponding to dpdx_pial, dpdx_pen, 
%     dpdx_cap 
%   Changed vessel_type to be a 1 X M vector, (M is the number of edges)
%   Replaced all calls to r_pen_art with d_pen_art/2
%   Make C_efflux & C_par input parameters (this way it will be easy to
%     drop in a computational calculation for either conductance)
% 18 Oct 2020 - (JT) added an example to the notes above, changed variable 
%   "is_term_cap" to "is_cap", removed "pf" input, removed two subfunctions
%   that are no longer used, added a "vessel_type" option 4 (which 
%   corresponds to parenchymal flow), added variable dp_paren, and added
%   param.kappa to model porous PVSs.
% 22 Oct 2020 - (KASB) reduced the number of input parameters by including
%   sc, dpdx, C_paren, and C_efflux in param
% 31 May 2020 - (JT) edited the hydraulic resistance calculation for
%   precapillary PVSs so that it's based on the analytical solution for
%   flow through an annulus, not the power law from Tithof et al 2019
% 1 June 2021 - (JT) simplified the code by removing the additional
%   "batteries" in series throughout the network; the network modeled in 
%   this code is only driven by one external pressure drop pext; this sped
%   the code up as well
%--------------------------------------------------------------------------

%% Set parameters
dpdx=param.dpdx;
dpdx_pial=dpdx(1); 
dpdx_pen=dpdx(2);
dpdx_cap=dpdx(3);
dp_paren=dpdx(4);
sc=param.sc;
C_paren=param.C_paren;
C_efflux=param.C_efflux;
g=param.g;
K=param.K;
d_pialart=param.d_pialart;
l_pial_art=param.l_pial_art;
l_pial2pen=param.l_pial2pen;
d_penart=param.d_penart;
l_penart=param.l_penart;
r_cap=param.r_cap;
l_cap=param.l_cap;
ncap_per_penart=param.ncap_per_penart;
kappa=param.kappa;

% Set the length of the penetrating arterioles and capillaries for the plot
% (each pial artery segment is about length 1)
l_penart_plot=5.7;
l_cap_plot=0.3;


%% Construct the variables that describe the connectivity and shape of the network

% Initialize variables
% "edges" contains the indices of the two nodes connecting that edge; 
% directionality is such that "forward" flow goes from column 1 node to 
% column 2 node
edges=[1 2]; % inlet; further down in the code, we place a "battery" before node 1
% "l_edge" is the length of each edge and should have the same number of
% elements as the number of rows in "edges"
l_edge=[l_pial_art];
% "d_artery" is the diameter of the artery for that corresponding edge; it
% should have the same number of elements as the number of rows in "edges"
d_edge=[d_pialart];
% "is_penart" is a binary variable (0 or 1) that indexes whether a given
% node is a penetrating arteriole or not; this vector has "n" entries
is_penart=[0 0];
% "is_cap" is a binary variable (0 or 1) that indexes whether a given
% node is a terminal capillary or not; this vector has "n" entries
is_cap=[0 0];
% "is_paren" is a binary variable (0 or 1) that indexes whether a given
% node corresponds to parenchymal flow or not; this vector has "n" entries
is_paren=[0 0];
% "active_nodes" are the most recently generated nodes from which branching
% will occur in the next iteration of the loop (this is a vector containing
% the indices of the nodes, not the total number of active nodes)
active_nodes=2;
% "n" is the current total number of nodes (a single number)
n=2;

% Loop over number of generations
for i=1:g
    
    %--------Draw half of hexagon for first outer-most active node---------
    prev_nodes=active_nodes;
    active_nodes=[];
    for j=prev_nodes(1)
                    
        % Add edges for half of a hexagon
        edges=[edges; j n+1; n+1 n+2; n+1 n+3; n+3 n+4; n+3 n+5;...
            n+5 n+6; n+6 n+7; n+6 n+8; n+8 n+9; n+8 n+10; n+10 n+11;...
            n+11 n+12; n+11 n+13; n+13 n+14; n+13 n+15];
        
        % Get three random edge lengths to use for the hexagon sides
        temp1=l_pial_art/3;temp2=l_pial_art/3;temp3=l_pial_art/3;

        % Update the length of all the hexagon sides and offshoots
        l_edge=[l_edge temp1 l_pial2pen temp1 l_pial2pen temp1 temp2...
            l_pial2pen temp2 l_pial2pen temp2 temp3 l_pial2pen temp3...
            l_pial2pen temp3];

        % Update the vessel diameters
        %d_edge=[d_edge repmat(d_pialart,1,15)];
        d_edge=[d_edge d_pialart d_penart d_pialart d_penart d_pialart ...
            d_pialart d_penart d_pialart d_penart d_pialart d_pialart ...
            d_penart d_pialart d_penart d_pialart];
        
        % Indicate which nodes connect to penetrating arteries/terminal
        % capillaries
        is_penart=[is_penart 0 1 0 1 0 0 1 0 1 0 0 1 0 1 0];
        is_cap=[is_cap zeros(1,15)];
        is_paren=[is_paren zeros(1,15)];
        
        % Update active nodes
        active_nodes=[active_nodes n+10 n+15];
        
        % Increment n
        n=n+15;
        
    end
    
    %------Draw interior Y-shaped region connecting adjacent hexagons------
    for j=prev_nodes(2:end-1)
        
        % Add edges for central Y-shaped portion
        edges=[edges; j n+1; n+1 n+2; n+1 n+3; n+3 n+4; n+3 n+5;...
            n+5 n+6; n+6 n+7; n+6 n+8; n+8 n+9; n+8 n;...
            n+5 n+10; n+10 n+11; n+10 n+12; n+12 n+13; n+12 n+14];

        % Get three random edge lengths to use for the hexagon sides
        temp1=l_pial_art/3;temp2=l_pial_art/3;temp3=l_pial_art/3;

        % Update the length of all the hexagon sides and offshoots
        l_edge=[l_edge temp1 l_pial2pen temp1 l_pial2pen temp1 temp2...
            l_pial2pen temp2 l_pial2pen temp2 temp3 l_pial2pen temp3...
            l_pial2pen temp3];

        % Update the vessel diameters
        %d_edge=[d_edge repmat(d_pialart,1,15)];
        d_edge=[d_edge d_pialart d_penart d_pialart d_penart d_pialart ...
            d_pialart d_penart d_pialart d_penart d_pialart d_pialart ...
            d_penart d_pialart d_penart d_pialart];
        
        % Indicate which nodes connect to penetrating arteries/terminal
        % capillaries
        is_penart=[is_penart 0 1 0 1 0 0 1 0 1 0 1 0 1 0];
        is_cap=[is_cap zeros(1,14)];
        is_paren=[is_paren zeros(1,14)];
        
        % Update active nodes
        active_nodes=[active_nodes n+14];
        
        % Increment n
        n=n+14;
        
    end
    
    %----------Add remaining half of last hexagon----------
    for j=prev_nodes(end)
        
        % Add edges for half of a hexagon
        edges=[edges; j n+1; n+1 n+2; n+1 n+3; n+3 n+4; n+3 n+5;...
            n+5 n+6; n+6 n+7; n+6 n+8; n+8 n+9; n+8 n+10; n+10 n+11;...
            n+11 n+12; n+11 n+13; n+13 n+14; n+13 n];

        % Get three random edge lengths to use for the hexagon sides
        temp1=l_pial_art/3;temp2=l_pial_art/3;temp3=l_pial_art/3;

        % Update the length of all the hexagon sides and offshoots
        l_edge=[l_edge temp1 l_pial2pen temp1 l_pial2pen temp1 temp2...
            l_pial2pen temp2 l_pial2pen temp2 temp3 l_pial2pen temp3...
            l_pial2pen temp3];

        % Update the vessel diameters
        %d_edge=[d_edge repmat(d_pialart,1,15)];
        d_edge=[d_edge d_pialart d_penart d_pialart d_penart d_pialart ...
            d_pialart d_penart d_pialart d_penart d_pialart d_pialart ...
            d_penart d_pialart d_penart d_pialart];
        
        % Indicate which nodes connect to penetrating arteries/terminal
        % capillaries
        is_penart=[is_penart 0 1 0 1 0 0 1 0 1 0 0 1 0 1];
        is_cap=[is_cap zeros(1,14)];
        is_paren=[is_paren zeros(1,14)];
        
        % Update active nodes
        active_nodes=[active_nodes n+10];
        
        % Increment n
        n=n+14;
        
    end

    
end

% Loop over penetrating arterioles and add several capillary and 
% parenchymal nodes
penart_ind=find(is_penart);
for i=1:length(penart_ind)
    
    %-------------penetrating artery nodes--------------------------
    
    % Add edges for the full penetrating artery, which capillaries will
    % branch off of
    edges=[edges;([penart_ind(i) n+1:n+ncap_per_penart-1])' (n+1:n+ncap_per_penart)'];
    
    % Set edge lengths for this penetrating artery
    l_edge=[l_edge repmat(l_penart/ncap_per_penart,1,ncap_per_penart)];
    
    % Update vessel diameter for this penetrating artery
    d_edge=[d_edge repmat(d_penart,1,ncap_per_penart)];
    
    % Indicate which nodes connect to penetrating arteries
    is_penart=[is_penart ones(1,ncap_per_penart)];
            
    % Indicate which nodes are capillaries
    is_cap=[is_cap zeros(1,ncap_per_penart)];
    
    % Indicate which nodes correspond to parenchymal flow
    is_paren=[is_paren zeros(1,ncap_per_penart)];
    
    %-----capillaries branching off the penetrating artery nodes-----------
     
    % Add each capillary
    edges=[edges;([n+1:n+ncap_per_penart])' (n+ncap_per_penart+1:n+2*ncap_per_penart)'];
    
    % Set edge lengths for this penetrating artery
    l_edge=[l_edge repmat(l_cap,1,ncap_per_penart)];
    
    % Update vessel diameter for this penetrating artery
    d_edge=[d_edge repmat(r_cap*2,1,ncap_per_penart)];
    
    % Indicate which nodes connect to penetrating arteries
    is_penart=[is_penart zeros(1,ncap_per_penart)];
    
    % Indicate which nodes are terminal capillaries
    is_cap=[is_cap ones(1,ncap_per_penart)];
    
    % Indicate which nodes correspond to parenchymal flow
    is_paren=[is_paren zeros(1,ncap_per_penart)];
    
    %-----parenchymal flow branching off the penetrating artery nodes-----------
     
    % Add each capillary
    edges=[edges;([n+1:n+ncap_per_penart])' (n+2*ncap_per_penart+1:n+3*ncap_per_penart)'];
    
    % Set edge lengths for this penetrating artery
    l_edge=[l_edge repmat(nan,1,ncap_per_penart)];
    
    % Update vessel diameter for this penetrating artery
    d_edge=[d_edge repmat(nan,1,ncap_per_penart)];
    
    % Indicate which nodes connect to penetrating arteries
    is_penart=[is_penart zeros(1,ncap_per_penart)];
    
    % Indicate which nodes are terminal capillaries
    is_cap=[is_cap zeros(1,ncap_per_penart)];
    
    % Indicate which nodes correspond to parenchymal flow
    is_paren=[is_paren ones(1,ncap_per_penart)];
    
    % Increment n
    n=n+3*ncap_per_penart;
        
end

% Assign vessel_type for each edge
% Loop over all the edges
vessel_type=zeros(size(l_edge));
for i=1:size(edges,1) 
    
    % Determine what type of vessel this is
    if is_paren(edges(i,2))
        vessel_type(i)=4;
    elseif is_cap(edges(i,2))
        vessel_type(i)=3;
    elseif is_penart(edges(i,1)) || is_penart(edges(i,2))
        vessel_type(i)=2;
    else
        vessel_type(i)=1;
    end
end

% Compute the cross sectional area of each PVS based on the PVS-to-artery
% area ratio K
PVS_area=pi*(d_edge/2).^2.*K(vessel_type);

% Compute the PVS volume of each segment
PVS_volume=PVS_area.*l_edge; 

%% Set up matrix and solve Kirchhoff's first law

% This is the version with "batteries" throughout the pial, penetrating,
% and capillary PVSs

% In the code below, we insert a separate "battery" in series with every 
% single resistor, immediately upstream of that resistor. The matrix
% structure for Cp=z looks approximately like:
% [ various sparse  |   sparse ] [p1]       [0]
% [     C vals      |    ones  ] [p2]       [0]
% [----------------------------] [...]  =   [...] 
% [     sparse      |   zeros  ] [Q12]      [dp12]
% [      ones       |          ] [Q23]      [dp23]
% [                 |          ] [...]      [...]
%
% The p vector is actually composed as:
%   p=[p1 p2 ... pn pn+1 ... pn+m-1 pn+m Q1 Q2 ... Qm-1 Qm]
%       for a total of n+2m nodes, where n is the (original) total number
%       of nodes and m is the number of edges
% Note that the node pn+m is typically grounded, so it could be removed
%   from the matrix

% Initialize a sparse matrix for constructing the circuit
%C=sparse(n+2*m+2,n+2*m+2);
C=sparse(n+2,n+2); % n nodes + 1 efflux node + 1 "battery"

% Initialize RHS vector in C*dp=z
% (z is approx. 2/3 zeros, 1/3 "battery" pressures)
z=zeros(size(C,1),1);

% Loop over all the edges
for i=1:size(edges,1) 

    % Reminder:
    % row determines node at which we're enforcing continuity
    % column determines which nodes are connected to that one

    % Compute conductance depending on vessel type and parameters
    if ~is_paren(edges(i,2))
        ctemp=compute_conductance(K,vessel_type(i),d_edge(i)/2,l_edge(i),kappa(vessel_type(i)));
    else
        ctemp=C_paren;
    end
    
    % Update C matrix
    C(edges(i,1),edges(i,2))=-ctemp; % negative upper triangle C value
    C(edges(i,2),edges(i,1))=-ctemp; % negative lower triangle C value
    C(edges(i,1),edges(i,1))=C(edges(i,1),edges(i,1))+ctemp; % positive diagonal C value
    C(edges(i,2),edges(i,2))=C(edges(i,2),edges(i,2))+ctemp; % positive diagonal C value

end

%% Construct variables that hold spatial coordinates of nodes for plotting

% Set coordinates for active node (the edge of what's just been drawn)
active_node_x=0;
active_node_y=0.5;

% Initialize variables for holding x and y coordinates of penetrating
% arterioles
xnode=[0 0];
ynode=[-0.5 0.5];
znode=[0 0];

% Loop over the number of generations
for i=1:g
    
    % Save a variable for the number of active nodes to remove after loop
    % is done
    anr=length(active_node_x);
    
    % Start with the left-most active node
    for k=1
        
        % Add left 3 sides of hexagon
        theta = [270 210 150 90];
        x=cosd(theta)+active_node_x(k);
        y=sind(theta)+active_node_y(k)+1;
        
        % Update active nodes
        active_node_x=[active_node_x x([3 4])];
        active_node_y=[active_node_y y([3 4])];
        
        % Compute vectors pointing from one node to the next
        dx=diff(x);
        dy=diff(y);
        
        % Add temporary node coordinates, sort, then add to array
        xnodetemp=[x(2:end) x(1:end-1)+dx/3 x(1:end-1)+2*dx/3];
        ynodetemp=[y(2:end) y(1:end-1)+dy/3 y(1:end-1)+2*dy/3];
        [~,ind]=sort(ynodetemp,'ascend'); % sort in ascending order
        xnode=[xnode xnodetemp(ind)]; % then add sorted nodes
        ynode=[ynode ynodetemp(ind)];
        znode=[znode zeros(1,length(ind))];
        
        % Compute perpendicular direction and add coordinates for
        % penetrating arterioles
        xhat=dx./sqrt(dx.^2+dy.^2);
        yhat=dy./sqrt(dx.^2+dy.^2);
        l=2/sqrt(3)/3;
        for j=1:length(xhat)
            perpdir=cross([0 0 1],[xhat(j) yhat(j) 0]);
            xnode=[xnode x(j)+dx(j)/3+(-1)^mod(j,2)*l*perpdir(1) x(j)+2*dx(j)/3+(-1)^mod(j,2)*l*perpdir(1)];
            ynode=[ynode y(j)+dy(j)/3+(-1)^mod(j,2)*l*perpdir(2) y(j)+2*dy(j)/3+(-1)^mod(j,2)*l*perpdir(2)];
            znode=[znode 0 0];
        end
    end
    
    % Reorganize xnode and ynode so that they're in the same order as
    % matrix construction above (no need to update znode -- all zeros)
    l=length(xnode);
    xnode(l-14:l)=xnode(l-15+[1 10 2 11 3 4 12 5 13 6 7 14 8 15 9]);
    ynode(l-14:l)=ynode(l-15+[1 10 2 11 3 4 12 5 13 6 7 14 8 15 9]);

    % Loop over every other active node in the middle (the "tops" of the
    % hexagons)
    for k=2:2:anr-1

        % Draw right and top right portion of hexagon
        theta = [330    30    90];
        x=cosd(theta)+active_node_x(k-1);
        y=sind(theta)+active_node_y(k-1)+1;
        
        % Update active nodes
        active_node_x=[active_node_x x([2])];
        active_node_y=[active_node_y y([2])];

        % Compute vectors pointing from one node to the next
        dx=diff(x);
        dy=diff(y);
        
        % Add temporary node coordinates, sort, then add to array
        xnodetemp=[x(2:end-1) x(1:end-1)+dx/3 x(1:end-1)+2*dx/3];
        ynodetemp=[y(2:end-1) y(1:end-1)+dy/3 y(1:end-1)+2*dy/3];
        [~,ind]=sort(ynodetemp,'ascend'); % sort in ascending order
        xnode=[xnode xnodetemp(ind)]; % then add sorted nodes
        ynode=[ynode ynodetemp(ind)];
        znode=[znode zeros(1,length(ind))];

        % Compute perpendicular direction and add penetrating arterioles
        xhat=dx./sqrt(dx.^2+dy.^2);
        yhat=dy./sqrt(dx.^2+dy.^2);
        l=2/sqrt(3)/3;
        for j=1:length(xhat)
            perpdir=cross([0 0 1],[xhat(j) yhat(j) 0]);
            xnode=[xnode x(j)+dx(j)/3+(-1)^mod(j+1,2)*l*perpdir(1) x(j)+2*dx(j)/3+(-1)^mod(j+1,2)*l*perpdir(1)];
            ynode=[ynode y(j)+dy(j)/3+(-1)^mod(j+1,2)*l*perpdir(2) y(j)+2*dy(j)/3+(-1)^mod(j+1,2)*l*perpdir(2)];
            znode=[znode 0 0];
        end

        % Add top left portion of next adjacent hexagon
        theta = [150 90];
        x=cosd(theta)+active_node_x(k+1);
        y=sind(theta)+active_node_y(k+1)+1;

        % Update active nodes
        active_node_x=[active_node_x x([2])];
        active_node_y=[active_node_y y([2])];

        % Draw points along hexagon edges
        dx=diff(x);
        dy=diff(y);

        % Add temporary node coordinates, sort, then add to array
        xnodetemp=[x(2) x(1:end-1)+dx/3 x(1:end-1)+2*dx/3];
        ynodetemp=[y(2) y(1:end-1)+dy/3 y(1:end-1)+2*dy/3];
        [~,ind]=sort(ynodetemp,'ascend'); % sort in ascending order
        xnode=[xnode xnodetemp(ind)]; % then add sorted nodes
        ynode=[ynode ynodetemp(ind)];
        znode=[znode zeros(1,length(ind))];

        % Compute perpendicular direction and add penetrating arterioles
        xhat=dx./sqrt(dx.^2+dy.^2);
        yhat=dy./sqrt(dx.^2+dy.^2);
        l=2/sqrt(3)/3;
        for j=1:length(xhat)
            perpdir=cross([0 0 1],[xhat(j) yhat(j) 0]);
            xnode=[xnode x(j)+dx(j)/3+(-1)^mod(j,2)*l*perpdir(1) x(j)+2*dx(j)/3+(-1)^mod(j,2)*l*perpdir(1)];
            ynode=[ynode y(j)+dy(j)/3+(-1)^mod(j,2)*l*perpdir(2) y(j)+2*dy(j)/3+(-1)^mod(j,2)*l*perpdir(2)];
            znode=[znode 0 0];
        end
        
        % Reorganize xnode and ynode so that they're in the same order as
        % matrix construction above (no need to update znode -- all zeros)
        l=length(xnode);
        xnode(l-13:l)=xnode(l-14+[1 6 2 7 3 4 8 5 9 10 13 11 14 12]);
        ynode(l-13:l)=ynode(l-14+[1 6 2 7 3 4 8 5 9 10 13 11 14 12]);

    end
    
    % Finally, add onto the right-most active node
    for k=anr
        
        % Add right 3 sides of hexagon
        theta = [270   330    30    90];
        x=cosd(theta)+active_node_x(k);
        y=sind(theta)+active_node_y(k)+1;
        
        % Update active nodes
        active_node_x=[active_node_x x([3])];
        active_node_y=[active_node_y y([3])];
        
        % Compute vectors pointing from one node to the next
        dx=diff(x);
        dy=diff(y);
        
        % Add temporary node coordinates, sort, then add to array
        xnodetemp=[x(2:end-1) x(1:end-1)+dx/3 x(1:end-1)+2*dx/3];
        ynodetemp=[y(2:end-1) y(1:end-1)+dy/3 y(1:end-1)+2*dy/3];
        [~,ind]=sort(ynodetemp,'ascend'); % sort in ascending order
        xnode=[xnode xnodetemp(ind)]; % then add sorted nodes
        ynode=[ynode ynodetemp(ind)];
        znode=[znode zeros(1,length(ind))];
        
        % Compute perpendicular direction
        xhat=dx./sqrt(dx.^2+dy.^2);
        yhat=dy./sqrt(dx.^2+dy.^2);
        l=2/sqrt(3)/3;
        for j=1:length(xhat)
            perpdir=cross([0 0 1],[xhat(j) yhat(j) 0]);
            xnode=[xnode x(j)+dx(j)/3+(-1)^mod(j,2)*l*perpdir(1) x(j)+2*dx(j)/3+(-1)^mod(j,2)*l*perpdir(1)];
            ynode=[ynode y(j)+dy(j)/3+(-1)^mod(j,2)*l*perpdir(2) y(j)+2*dy(j)/3+(-1)^mod(j,2)*l*perpdir(2)];
            znode=[znode 0 0];
        end
    end
    
    % Reorganize xnode and ynode so that they're in the same order as
    % matrix construction above (no need to update znode -- all zeros)
    l=length(xnode);
    xnode(l-13:l)=xnode(l-14+[1 9 2 10 3 4 11 5 12 6 7 13 8 14]);
    ynode(l-13:l)=ynode(l-14+[1 9 2 10 3 4 11 5 12 6 7 13 8 14]);
    
    % Remove the active nodes corresponding to this current loop iteration
    active_node_x(1:anr)=[];
    active_node_y(1:anr)=[];
end

% Loop over penetrating arterioles and add nodes below the brain's surface
penart_ind=find(is_penart);
for i=1:length(penart_ind)
    
    % Only add coordinates for penetrating arteries with z==0
    if znode(penart_ind(i))==0
        %-------------penetrating artery nodes--------------------------
        xnode=[xnode repmat(xnode(penart_ind(i)),1,ncap_per_penart)];
        ynode=[ynode repmat(ynode(penart_ind(i)),1,ncap_per_penart)];
        znode=[znode -linspace(l_penart_plot/ncap_per_penart,l_penart_plot,ncap_per_penart)];

        %-----capillaries branching off the penetrating artery nodes----------- 
        %randang=linspace(0,2*pi,ncap_per_penart+1);%rand(1,ncap_per_penart)*2*pi;
        randang=0:110*pi/180:110*pi/180*ncap_per_penart;%rand(1,ncap_per_penart)*2*pi;
        randang=randang(1:end-1);
        xnode=[xnode xnode(penart_ind(i))+cos(randang)*l_cap_plot];
        ynode=[ynode ynode(penart_ind(i))+sin(randang)*l_cap_plot];
        znode=[znode -linspace(l_penart_plot/ncap_per_penart,l_penart_plot,ncap_per_penart)];
        
        %-----parenchymal nodes branching off the penetrating artery nodes-----------
        %randang=linspace(0,2*pi,ncap_per_penart+1);%rand(1,ncap_per_penart)*2*pi;
        randang=0:110*pi/180:110*pi/180*ncap_per_penart;%rand(1,ncap_per_penart)*2*pi;
        randang=randang(1:end-1)+pi;
        xnode=[xnode xnode(penart_ind(i))+cos(randang)*l_cap_plot];
        ynode=[ynode ynode(penart_ind(i))+sin(randang)*l_cap_plot];
        znode=[znode -linspace(l_penart_plot/ncap_per_penart,l_penart_plot,ncap_per_penart)];
    end
        
end

%% Now that geometry is drawn, finish computation

% Connect capillaries to common low-resistance efflux node (node n+1)
ind=find(is_cap);
for i=1:length(ind)
    C(ind(i),n+1)=-C_efflux;
    C(n+1,ind(i))=-C_efflux;
    C(ind(i),ind(i))=C(ind(i),ind(i))+C_efflux;
    C(n+1,n+1)=C(n+1,n+1)+C_efflux;
end

% Connect parenchymal nodes to common low-resistance efflux node (node n+1)
ind=find(is_paren);
for i=1:length(ind)
    C(ind(i),n+1)=-C_efflux;
    C(n+1,ind(i))=-C_efflux;
    C(ind(i),ind(i))=C(ind(i),ind(i))+C_efflux;
    C(n+1,n+1)=C(n+1,n+1)+C_efflux;
end


%---------------This part is for a prescribed pressure----------------

% Now connect a pressure source (due to the choroid plexus) between the 
% last node and first node
C(1,end)=-1; % this specifies the voltage into the first node
C(end,1)=1; % this specifies voltage of this battery
C(end,n+1)=-1;
C(n+1,end)=-1;

% Set the last entry in z, which is the external pressure driving the flow.
z(end)=p0;

% Ground the efflux node
C(:,n+1)=[]; 
C(n+1,:)=[]; 
z(n+1)=[];
% there are now n+1 nodes

% Resulting pressures should be scaled up by 1e10
% diagonalize [C|z] to get [p;Q]
% z is order 1
% C goes from 1e-6 to 1e5
% pressures go from 1e-15 to 1e-1
% Q goes from 1e-19 1e-4

% Now, finally, compute the pressures and total volume flow rate by 
% computing the reduced-row echelon form of [C|z]
% tic;
%R=frref([C(:,[1:end-1]) z]);
R=frref([C z]);
% toc
pressures=full(R(1:n,end))'; % pressure at each of the original nodes; units of mmHg
Cind=sub2ind(size(C),edges(:,1),edges(:,2));
Q=C(Cind).*(pressures(edges(:,2))-pressures(edges(:,1)))'; % volume flow rate through each edge (except the one between ground and node 1)
Qtotal=full(R(end,end)); % volume flow into node 1 (which equals total volumetric flow)
Rtotal=p0/Qtotal;
% disp(['Total volumetric flow is ' sprintf('%1.10f',Qtotal) ' mL/min']);
%disp(['Total hydraulic resistance is ' sprintf('%1.10f',p0/Qtotal) ' mmHg*min/mL']);

% Faghih & Sharp found 1.14 mmHg/ml/min resistance for paraarterial flow
% Some quotes:
% 0.13 ml/min through the paraarterial tree. This lowest required pressure 
% difference was 0.15 mmHg.
% The total resistance of the paravenous tree was equal to 1.75 × 10^3 
% mmHg/ml/min, about three orders of magnitude smaller than that of the 
% paraarterial tree (Fig. 5), which can be expected based on the larger 
% gaps and larger vessel diameters compared to the paraarterial channels. 
% A more consistent generation-to-generation increase in resistance is also
% evident. The required pressure difference to drive 0.13 ml/min of flow 
% through the paravenous tree was calculated to be 0.00023 mmHg.


%% Save files, if a savename is given
if ~isempty(savename)
    save([savename '.mat']);
end

end

function C=compute_conductance(K,t,r,l,kappa)
% This function takes the following inputs:
%   K - ratio of PVS-to-artery cross sectional areas; usually 1.4
%   t - type of vessel; 1 - pial, 2 - penetrating, 3 - capillary
%   r - artery radius (m)
%   l - length of segment (m)
% and returns the conductance (inverse of hydraulic resistance) in units
% of mL / (mmHg*min)

% Set viscosity (kg / m s); rho=993 kg/m^3 and nu=7e-7 m^2/s
mu=993*7e-7;
K=K(t);

% Conversion factor for kg/(m^5 s) to mmHg*min/mL / m
cf=1.25e-10;

% Calculate conductance for each possible type of vessel
% Expression for R come from Tithof et al, Fluids Barriers CNS 2019
% Assuming r is in meters, R is in units of kg/(m^5 s), ie it's the
% hydraulic resistance per unit length
switch t
    case 1 % pial artery (optimal elliptical annulus)
        if isnan(kappa)
            R=cf*mu*6.67*K^(-1.96)/r^4; % (units mL / (mmHg*min*m))
            C=1/(R*l); % (units mmHg*min/mL)
        else
            r2=r*sqrt(K+1); % radius of outer wall for circular annulus
            R=cf*mu/(pi*kappa*(r2^2-r^2)); % (units mL / (mmHg*min*m))
            C=1/(R*l); % (units mmHg*min/mL)
        end
    case 2 % penetrating artery (tangent eccentric annulus)
        if isnan(kappa)
            R=cf*mu*8.91*K^(-2.78)/r^4; % (units mL / (mmHg*min*m))
            C=1/(R*l); % (units mmHg*min/mL)
        else
            r2=r*sqrt(K+1); % radius of outer wall for circular annulus
            R=cf*mu/(pi*kappa*(r2^2-r^2)); % (units mL / (mmHg*min*m))
            C=1/(R*l); % (units mmHg*min/mL)
        end
    case 3 % capillary basement membrane (concentric annulus)
        if isnan(kappa)
            r2=r*sqrt(K+1); % radius of outer wall for circular annulus
            R=cf*8*mu/pi/(r2^4-r^4-(r2^2-r^2)^2/log(r2/r)); % (units mL / (mmHg*min*m))
            C=1/(R*l); % (units mmHg*min/mL)
        else
            r2=r*sqrt(K+1); % radius of outer wall for circular annulus
            R=cf*mu/(pi*kappa*(r2^2-r^2)); % (units mL / (mmHg*min*m))
            C=1/(R*l); % (units mmHg*min/mL)
        end
end

end

