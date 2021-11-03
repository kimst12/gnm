% code the iteratively solves the model to find what driving pressure pext
% is required so that the medial pial velocity is 18.7 microns/s. Requires
% the following variables be defined: (example values provided)
% sleep_or_awake='sleep';
% paren_type='high_res';
% pen_perm=nan;
% cap_perm=1.8e-14;
% cap_ar=0.07;


guessP=[0 80];
fun = @(p0) g(p0,pen_perm,cap_perm,cap_ar,paren_type,sleep_or_awake);
options = optimset('Display','iter','TolX',1e-6); % show iterations
pext = fzero(fun,guessP,options);

function z=g(p0,pen_perm,cap_perm,cap_ar,paren_type,sleep_or_awake)
    filename='findP_results';
    
    [param] = setParam(pen_perm,cap_perm,cap_ar,paren_type,sleep_or_awake);

    %% Run the hydraulic network model code
    [~,~]=branching_hexagon_model_pext(filename,p0,param);

    load(filename)

    cf=1.667e-8; % 1 mL/min = 1.667e-8 m^3/s
    all_velocities=abs(Q')./PVS_area*cf; % This will return "nan" for parenchymal nodes
    pial_inds=find(~is_penart(edges(:,1)) & ~is_penart(edges(:,2)) & ~is_cap(edges(:,1)) & ~is_cap(edges(:,2)) & ~is_paren(edges(:,2)));
    vpial=median(all_velocities(pial_inds));
    disp(['vpial is ' num2str(vpial*1e6) ' microns/s'])
    disp(['p0 is ' num2str(p0) ' mmHg'])

    z=18.7e-6-vpial; 
end
