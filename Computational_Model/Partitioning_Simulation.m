%----------------------------------------------------------------------------
% This code was developed to generate Figure 6 and Supplementary Figure 12
% of the following paper:
%
%
%
%----------------------------------------------------------------------------
% Tatsat Banerjee, Satomi Matsuoka, Debojyoti Biswas, 
% Yuchuan Miao, Dhiman Sankar Pal, Yoichiro Kamimura, Masahiro Ueda, 
% Peter N. Devreotes, Pablo A. Iglesias. 
% 
% "A dynamic partitioning mechanism polarizes membrane protein distribution".
% 
% 
% Nature Communications, 14, 7909 (2023). DOI:10.1038/s41467-023-43615-2. 
% 
% PMID: 38036511. 
%
%
%----------------------------------------------------------------------------


%% function Partitioning_Simulation

    clear
    close all
    clc

    format short
    warning('off','all')

%% adding paths

    path_to_urdme = fullfile(pwd,'/urdme_local/');
    addpath(genpath([path_to_urdme 'urdme/']));
%%
tic
   clear umod Umod1 param
%% Time  

    param.dt = 0.1;

    param.T_vec = [30 30 30];
    for i  = 1:length(param.T_vec)
        param.tspan{i} = 0:param.dt:param.T_vec(i);     
    end
    param.Time = 0:param.dt:sum(param.T_vec);

%
    
    param.Na = 6.022e23; %Avogadro's Number

%% build the geometry

    param.L = 20;% 20 micron
    param.x1 = 0; param.x2 = param.L;
    param.y1 = 0; param.y2 = param.x2;
    param.dx = 0.25;

    S1 = [3 4 param.x1 param.x2 param.x2 param.x1 ...
              param.y1 param.y1 param.y2 param.y2]'; % Square, 20x20um
    gd = S1;
    sf = 'S1';
    ns = char('S1')';
    G  = decsg(gd,sf,ns);

% create the mesh
    [P,E,T] = initmesh(G,'Hmax',param.dx);
    param.P = P;
    param.E = E;
    param.T = T;


%% Assemble the diffusion part

% diffusion constants are in micron^2/sec

param.Mspecies = 17;

param.D_F = 0.15;     % diffusion constant of F 
param.D_R = 0.09;     % diffusion constant of R
param.D_B = 0.075;    % diffusion constant of B

% S1c: PP_c, S1: PP_u, BS1: PP_b
param.D_S1c = 0.75;        % diffusion constant of cytosolic PP
param.D_S1  = 0.45;        % diffusion constant of unbound PP
param.D_BS1 = 0.05;        % diffusion constant of B:PP
param.D_FS1 = 0.05;        % diffusion constant of F:PP

% S2: LP_u, BS2: LP_b, 
param.D_S2  = 0.45;  % diffusion constant of unbound LP
param.D_BS2 = 0.05;  % diffusion constant of B:LP
param.D_FS2 = 0.05;  % diffusion constant of F:LP


% S3c: PP_c, S3: photoconverted PP_u, BS3: photoconverted PP_b
param.D_S3c = 0.75;       
param.D_S3  = 0.45;                    
param.D_BS3 = 0.05; 
param.D_FS3 = 0.05; 

% S4: photoconverted LP_u, BS4: photoconverted LP_b, 
param.D_S4  = 0.45;   
param.D_BS4 = 0.05;  
param.D_FS4 = 0.05;  

D_vector = [param.D_F,   param.D_B,   param.D_R,...
            param.D_S1c, param.D_S1,  param.D_BS1, param.D_FS1,...
            param.D_S2,  param.D_BS2, param.D_FS2,...
            param.D_S3c, param.D_S3,  param.D_BS3, param.D_FS3,...
            param.D_S4,  param.D_BS4, param.D_FS4];

Diffusion_constants = num2cell(D_vector);
umod = pde2urdme(param.P,param.T,Diffusion_constants);

umod.pde.E = param.E;
% umod.private.Mspecies = length(Diffusion_constants);



%% reaction parameters

    % Vary these parameters to achieve intended number of molecules and time response


    % response scaling parameters
    alpha =  0.05;     
    beta = 0.05;      
    gamma = 0.05;     
    epsilon = 2;   
    
 
    % parameters for F dynamics
    param.a1 =  0.0830;  
    param.a2 =  8.3300;  
    param.a3 = 93.7000;  
    param.a4 =  144;  
    param.a5 =  40.000;  
    param.ub =  0.1470;  
    
    % parameters for B dynamics
    param.b1 =  10; 
    param.b2 =  0.0100;
    param.b3 =  10;
     
    
    % parameters for R dynamics
    param.c1 =  0.01;
    param.c2 = 0.64; 
    
%--------------------------------------------------------------------------
    
  % parameters for LP dynamics
    
    param.e4 =   72; 
    fracn = 0.3414;
    mean_B_level = 20;
    param.e3 =  6*(fracn/(1-fracn))*mean_B_level; 

    param.e5 =   0.05;
    param.e6 =   10;


    % parameters for PP dynamics


    param.d1 =  1;
    param.d2 =  8;    
    param.d3 =  1.6*(fracn/(1-fracn))*mean_B_level;
    param.d4 =  0.8;
%--------------------------------------------------------------------------    
    
 
    % Scaling parameters
    param.alpha = alpha;
    param.beta  = beta;
    param.gamma = gamma;
    param.epsilon = epsilon;
       
%%
    param.run_no = 1;
    sd_tmp = ones(1,size(param.P,2));
    umod.sd = sd_tmp;

    

%% load initial data

init_cond_tmp = load("Initial_Condition_Data.mat");
init_cond_tmp = init_cond_tmp.init_cond;

init_cond = [init_cond_tmp(1:6,:); ones(1,size(init_cond_tmp,2));...
             init_cond_tmp(7:8,:);ones(1,size(init_cond_tmp,2));...
             init_cond_tmp(9:11,:);zeros(1,size(init_cond_tmp,2));...
             init_cond_tmp(12:13,:);ones(1,size(init_cond_tmp,2))];

field_names = {'F','B','R', 'S1c','S1','BS1','FS1',...
                                  'S2','BS2','FS2',...
                            'S3c','S3','BS3','FS3',...
                                  'S4','BS4','FS4'};


%%

    umod.private.param = param;

    Umod = cell(1,length(param.T_vec));
    processed_data = cell(1,length(param.T_vec));
    seed_val = 2144; % fixed for reporducibility


    for k = 1:length(param.tspan)

            Umod{k} = umod;
            Umod{k}.private.run_no = k;


            if k == 1

                 Umod{k}.private.init_cond = init_cond;
                 sd_tmp = ones(1,size(param.P,2));
                 Umod{k}.sd = sd_tmp;

            else

                Umod{k}.private.init_cond      = create_init_cond(processed_data{k-1},field_names);
                
%                 if k == 2
%                    load("PC_indices.mat"); % load spatial indices for the photoconversion in the intended region
%                    Node_index_pc = unique(Node_index_pc);
%                    sd_tmp = ones(1,size(param.P,2));
%                    sd_tmp(Node_index_pc) = 2;
%                    Umod{k}.sd = sd_tmp;
                %else
                   sd_tmp = ones(1,size(param.P,2));
                   Umod{k}.sd = sd_tmp;
                %end
            end

            if isfield(param, 'tspan') && any(param.tspan{k})
                Umod{k}.tspan = param.tspan{k};
            end  

            Umod{k}         = reaction_model(Umod{k},param);

            Umod{k} = urdme(Umod{k},'seed',seed_val,'propensities','reaction_model','report',0);


            processed_data{k} = data_processing(Umod{k},field_names);
        
     end 
          
%% Combining the data
close all
clc

        plot_data = 1;
        comb_data = combining_data(processed_data,umod,field_names,param,plot_data);  

elapsed_time = toc;
fprintf('Elapsed time = %f in min \n',elapsed_time/60)


%%
[Node_index_front, Node_index_back] = front_back(comb_data,param);


%%

figure('color','white')
set(gca,'linewidth',1.5,'fontsize',14)
hold on
h(1) = plot(comb_data.tspan(1:10:end),(comb_data.sum_F_molecule(1:10:end)/max(comb_data.sum_F_molecule)));
h(2) = plot(comb_data.tspan(1:10:end),comb_data.sum_B_molecule(1:10:end)/max(comb_data.sum_B_molecule));
h(3) = plot(comb_data.tspan(1:10:end),comb_data.sum_R_molecule(1:10:end)/max(comb_data.sum_R_molecule));
ylim([0 1.1])
legend(h,{'F','B','R'},'Location','best','EdgeColor','none')
xlabel('Time (s)')
ylabel('Norm. total conc.')


S1_membrane = comb_data.sum_S1_molecule + comb_data.sum_BS1_molecule + comb_data.sum_FS1_molecule;
S2_membrane = comb_data.sum_S2_molecule + comb_data.sum_BS2_molecule + comb_data.sum_FS2_molecule;
S1_membrane = S1_membrane/max(S1_membrane);
S2_membrane = S2_membrane/max(S2_membrane);

figure('color','white')
set(gca,'linewidth',1.5,'fontsize',14)
hold on
h(1) = plot(comb_data.tspan(1:10:end),(S1_membrane(1:10:end)));
h(2) = plot(comb_data.tspan(1:10:end),(S2_membrane(1:10:end)));
ylim([0.5 1.1])
legend(h,{'PP','LP'},'Location','best','EdgeColor','none')
xlabel('Time (s)')
ylabel('Norm. total conc.')

%%

    X = 10;
    Y = 10;
    
    p = [X',Y'];
    spatial_points = dsearchn(P',p);
    index1 = 1;

    
    F_smooth = smoothdata(comb_data.F_molecule(spatial_points,index1:end),'movmean',11);
    B_smooth = smoothdata(comb_data.B_molecule(spatial_points,index1:end),'movmean',11);
    R_smooth = smoothdata(comb_data.R_molecule(spatial_points,index1:end),'movmean',11);

    S1c_smooth = smoothdata(comb_data.S1c_molecule(spatial_points,index1:end),'movmean',11);
    S1_smooth = smoothdata(comb_data.S1_molecule(spatial_points,index1:end),'movmean',11);
    BS1_smooth = smoothdata(comb_data.BS1_molecule(spatial_points,index1:end),'movmean',11);
    FS1_smooth = smoothdata(comb_data.FS1_molecule(spatial_points,index1:end),'movmean',11);

    S2_smooth = smoothdata(comb_data.S2_molecule(spatial_points,index1:end),'movmean',11);
    BS2_smooth = smoothdata(comb_data.BS2_molecule(spatial_points,index1:end),'movmean',11);
    FS2_smooth = smoothdata(comb_data.FS2_molecule(spatial_points,index1:end),'movmean',11);

    S1_total = S1_smooth + BS1_smooth + FS1_smooth + S1c_smooth;
    S2_total = S2_smooth + BS2_smooth + FS2_smooth;
    
    figure('color','white')    
    hold on
    set(gca,'linewidth',1.5,'fontsize',14)
    h(1) = plot(comb_data.tspan,F_smooth(:,:)/max(F_smooth(:)),'g-','linewidth',1);
    h(2) = plot(comb_data.tspan,B_smooth(:,:)/max(B_smooth(:)),'r-','linewidth',1);
    h(3) = plot(comb_data.tspan,R_smooth(:,:)/max(R_smooth(:)),'b-','linewidth',1);
    ylim([0 1])
    set(gca,'TickDir','out');
    legend(h,{'F','B','R'},'Location','best','EdgeColor','none')
    xlabel('Time (s)')
    ylabel('Norm. total conc. at node a')

    figure('color','white')    
    hold on
    set(gca,'linewidth',1.5,'fontsize',14)
    h(1) = plot(comb_data.tspan,S2_smooth(:,:)/max(S2_total(:)),'g-','linewidth',1);
    h(2) = plot(comb_data.tspan,(BS2_smooth(:,:)+FS2_smooth(:,:))/max(S2_total(:)),'r-','linewidth',1);
    ylim([0 1])
    set(gca,'TickDir','out');
    legend(h,{'LP unbound at membrane','LP bound on membrane'},'Location','best','EdgeColor','none')
    xlabel('Time (s)')
    ylabel('Norm. total conc. at node a')
    

    figure('color','white')    
    hold on
    set(gca,'linewidth',1.5,'fontsize',14)
    plot(comb_data.tspan,S1c_smooth(:,:)/max(S1_total(:)),'b-','linewidth',1)
    plot(comb_data.tspan,S1_smooth(:,:)/max(S1_total(:)),'g-','linewidth',1)   
    plot(comb_data.tspan,(BS1_smooth(:,:)+FS1_smooth(:,:))/max(S1_total(:)),'r-','linewidth',1)
    ylim([0 1])
    set(gca,'TickDir','out');
    legend(h,{'PP at cytosol','PP unbound at membrane','PP bound on membrane'},'Location','best','EdgeColor','none')
    xlabel('Time (s)')
    ylabel('Norm. total conc. at node a')

%%
    data_visualization(comb_data,field_names,param)

%%  Auxilliary Fuctions 

function umod = reaction_model(umod,param)

% transitions and rates

Na = 6.022e23; % Avogradro's Number
k  = 1e-21*Na; % volume to molecule conversion factor
Ts = 0.2; % time scaling for FBR dyanamics


%--------------------------------------------------------------------------
% F Dynamics
%--------------------------------------------------------------------------
r1  = 'F > a1*F > @';
r2  = 'F > (a2/vol)*F*R > @';
r3 = '@ > (a3/(a4*a4/vol/vol*B*B+1))*(a5*vol-F) > F';
r4 = '@ > ub*(a5*vol-F) > F';
%--------------------------------------------------------------------------
a1  = param.a1*Ts;
a2  = param.a2/k*Ts;
a3  = param.a3*Ts;
a4  = param.a4/k;
a5  = param.a5*k;
ub  = param.ub*Ts;
%--------------------------------------------------------------------------
% R Dynamics
%--------------------------------------------------------------------------
r5  = 'R > c1*R > @';
r6  = '@ > c2*F > R';
%--------------------------------------------------------------------------
c1 = param.c1*Ts;
c2 = param.c2*Ts;
%--------------------------------------------------------------------------
% B Dynamics
%--------------------------------------------------------------------------
r7  = '@ > b1*vol       > B';
r8  = 'B > (b2)*B       > @';
r9  = 'B > (b3/vol)*F*B > @';
%--------------------------------------------------------------------------
b1 = param.b1*k*Ts;
b2 = param.b2*Ts;
b3 = param.b3/k*Ts;
%--------------------------------------------------------------------------
% PP Dynamics
%--------------------------------------------------------------------------
r10  = 'S1c > d1*S1c > S1';
r11  = 'S1 > (d2)*S1    > S1c';
r12  = 'S1 > (d3/vol)*B*S1 > BS1';
r13  = 'BS1 > (d4)*BS1 > S1';
r14  = 'S1 > 0.1*(d3/vol)*F*S1 > FS1';% 0.01
r15  = 'FS1 > 0.1*(d4)*FS1 > S1';
%--------------------------------------------------------------------------
d1 = param.d1;
d2 = param.d2;
d3 = param.d3/k;
d4 = param.d4;
%--------------------------------------------------------------------------
% LP Dynamics
%--------------------------------------------------------------------------
r16  = 'S2 > 2*(e3/vol)*B*S2 > BS2';
r17  = 'BS2 > 1.7*(e4)*BS2 > S2';
r18  = 'S2 > (e5/vol)*F*S2 > FS2';
r19  = 'FS2 > 2.6*(e6)*FS2 > S2';
%--------------------------------------------------------------------------
e3 = param.e3/k;
e4 = param.e4;
e5 = param.e5/k;
e6 = param.e6;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Photoconversion of PP
%--------------------------------------------------------------------------
r20  = 'S1c > sd == 2 ? 2*kon*S1c : 0.0 > S3c';
r21  = 'S1  > sd == 2 ? 2*kon*S1  : 0.0 > S3';
r22  = 'BS1 > sd == 2 ? 2*kon*BS1 : 0.0 > BS3';
r23  = 'FS1 > sd == 2 ? 2*kon*FS1 : 0.0 > FS3';
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Photoconversion of LP
%--------------------------------------------------------------------------
r24  = 'S2  > sd == 2 ? 2*kon*S2  : 0.0 > S4';
r25  = 'BS2 > sd == 2 ? 2*kon*BS2 : 0.0 > BS4';
r26  = 'FS2 > sd == 2 ? 2*kon*FS2 : 0.0 > FS4';
%--------------------------------------------------------------------------

kon = 0;

%--------------------------------------------------------------------------
% S3:photoconverted PP Dynamics
%--------------------------------------------------------------------------
r27  = 'S3c > d1*S3c > S3';
r28  = 'S3 > (d2)*S3    > S3c';
r29  = 'S3 > (d3/vol)*B*S3 > BS3';
r30  = 'BS3 > (d4)*BS3 > S3';
r31  = 'S3 > 0.1*(d3/vol)*F*S3 > FS3';
r32  = 'FS3 > 0.1*(d4)*FS3 > S3';
%--------------------------------------------------------------------------
% S4:photoconverted LP Dynamics
%--------------------------------------------------------------------------
r33  = 'S4 > 2*(e3/vol)*B*S4 > BS4';
r34  = 'BS4 > 1.5*(e4)*BS4 > S4';
r35  = 'S4 > (e5/vol)*F*S4 > FS4';
r36  = 'FS4 > 2.4*(e6)*FS4 > S4';



%--------------------------------------------------------------------------
% Adjustments
%--------------------------------------------------------------------------
alpha = param.alpha;
beta  = param.beta;
gamma = param.gamma;
epsilon = param.epsilon;

a2 = a2/gamma;
a4 = a4/beta;
a5 = a5*alpha;
b1 = b1*beta;
b3 = b3/alpha;
c1 = c1*epsilon;
c2 = (c2*gamma/alpha)*epsilon;


%------------------------------------------------------------------------
    [~,umod.N,umod.G] = ...
        rparse(...
        {r1 r2 r3 r4 r5 r6 r7 r8 r9 ...
         r10 r11 r12 r13 ...
         r14 r15 ...
         r16 r17 r18 r19 r20 ...
         r21 r22 r23 r24 ...
         r25 r26 ...
         r27 r28 r29 r30 r31 r32 ...
         r33 r34 r35 r36},...
        {'F' 'B' 'R' ...
         'S1c' 'S1' 'BS1' 'FS1'...
         'S2' 'BS2' 'FS2'...
         'S3c' 'S3' 'BS3' 'FS3'...
         'S4' 'BS4' 'FS4'}, ...
        {'a1' a1 'a2' a2 'a3' a3 'a4' a4 'a5' a5 'ub' ub ...
         'b1' b1 'b2' b2 'b3' b3 ...
         'c1' c1 'c2' c2 ...
         'd1' d1 'd2' d2 'd3' d3 'd4' d4 ...
         'e3' e3 'e4' e4 'e5' e5 'e6' e6...
         'kon' kon}, ...
         'reaction_model.c');


    if isfield(umod.private, 'init_cond')    
        umod.u0 = umod.private.init_cond;
    else
        umod.u0 = zeros(size(umod.N,1),numel(umod.vol));
    end   


    % simulation time interval
    if ~isfield(umod,'tspan')
        umod.tspan = 0:0.1:1;
    end



end

function processed_data = data_processing(umod,field_names)

        Na = 6.022e23; %Avogadro's Number
        n_var = umod.private.param.Mspecies;
        T_sampling = 1;
        
        
        processed_data.Mspecies       = n_var;
        
        
        
        
        if umod.private.run_no == 1
            
            tspan = umod.tspan;
            processed_data.tspan = tspan;
            
            for k = 1:length(field_names)
                
                field_name_tmp1 = strcat(field_names{k},'_molecule');
                processed_data.(field_name_tmp1)   = umod.U(k:n_var:end,1:T_sampling:end);
                field_name_tmp2 = strcat('sum_',field_names{k},'_molecule');
                processed_data.(field_name_tmp2) = sum(processed_data.(field_name_tmp1));
                field_name_tmp3 = strcat(field_names{k},'_molecule_tf');
                processed_data.(field_name_tmp3) = umod.U(k:n_var:end,end);

                field_name_tmp4 = strcat(field_names{k},'_conc');
                processed_data.(field_name_tmp4)   = diag(1./umod.vol)*umod.U(k:n_var:end,1:T_sampling:end)/Na*1e21;

                 
                
            end

        else

            tspan = umod.tspan(2:end);
            processed_data.tspan = tspan;

            for k = 1:length(field_names)
                
                field_name_tmp1 = strcat(field_names{k},'_molecule');
                processed_data.(field_name_tmp1)   = umod.U(k:n_var:end,2:T_sampling:end);
                field_name_tmp2 = strcat('sum_',field_names{k},'_molecule');
                processed_data.(field_name_tmp2) = sum(processed_data.(field_name_tmp1));
                field_name_tmp3 = strcat(field_names{k},'_molecule_tf');
                processed_data.(field_name_tmp3) = umod.U(k:n_var:end,end);
                field_name_tmp4 = strcat(field_names{k},'_conc');
                processed_data.(field_name_tmp4)   = diag(1./umod.vol)*umod.U(k:n_var:end,2:T_sampling:end)/Na*1e21;

                
            end      
        end
     

end

function init_cond = create_init_cond(processed_data,field_names)

   
    init_cond = [];
    for k = 1:length(field_names)
        field_name_tmp3 = strcat(field_names{k},'_molecule_tf');
        init_cond = [init_cond; (processed_data.(field_name_tmp3))'];

    end


end

function comb_data = combining_data(processed_data,umod,field_names,param,plot_data)



    comb_data.Mspecies = processed_data{1}.Mspecies;
   
    t_max = 0;
    for k = 1:length(processed_data)
        t_max = t_max + processed_data{k}.tspan(end);
    end
    comb_data.tspan = 0:umod.private.param.dt:t_max;

    for  K =  1:length(field_names)
    
        field_name_tmp1 = strcat(field_names{K},'_molecule');
        comb_data.(field_name_tmp1) = [];
        for k = 1:length(processed_data)
            comb_data.(field_name_tmp1) = [comb_data.(field_name_tmp1), processed_data{k}.(field_name_tmp1)];
        end


         field_name_tmp4 = strcat(field_names{K},'_conc');
         comb_data.(field_name_tmp4) = [];
         for k = 1:length(processed_data)
            comb_data.(field_name_tmp4) = [comb_data.(field_name_tmp4), processed_data{k}.(field_name_tmp4)];
        end

    end

   
    
    for  K =  1:length(field_names)
        
        field_name_tmp1 = strcat(field_names{K},'_molecule');
        field_name_tmp2 = strcat('sum_',field_names{K},'_molecule');
        field_name_tmp3 = strcat(field_names{K},'_molecule_tf');
        
        comb_data.(field_name_tmp2) = sum(comb_data.(field_name_tmp1));
        comb_data.(field_name_tmp3) = (comb_data.(field_name_tmp1 )(:,end))';  
       
    end
   

end

function data_visualization(comb_data,field_names,param)

    time_stamp = char(datetime);
    time_stamp_for_file = [time_stamp(1:2),'_',time_stamp(4:6),'_',time_stamp(8:11),'_',...
                                time_stamp(13:14),'_',time_stamp(16:17),'_',time_stamp(19:20)];
                  
    Data_folder = fullfile('Fig6new',time_stamp_for_file);
    if ~exist(Data_folder,'dir')
        mkdir(Data_folder)
    end
 

    n = length(comb_data.F_conc(:));
    nlow  = round(0.02*n);
    nhigh = round(0.98*n);
    
    Ftmpsorted = sort(comb_data.F_conc(:));
    F_low  = Ftmpsorted(nlow);
    F_high = Ftmpsorted(nhigh);
%--------------------------------------------------------------------------           

    Bmin = min(min(comb_data.B_conc));
    Rmin = min(min(comb_data.R_conc));

    S1cmin = min(min(comb_data.S1_conc));
    S1min = min(min(comb_data.S1_conc));
    BS1min = min(min(comb_data.BS1_conc));
    FS1min = min(min(comb_data.FS1_conc));

    S2min = min(min(comb_data.S2_conc));
    BS2min = min(min(comb_data.BS2_conc));
    FS2min = min(min(comb_data.FS2_conc));

    S3cmin = min(min(comb_data.S3c_conc));
    S3min = min(min(comb_data.S3_conc));
    BS3min = min(min(comb_data.BS3_conc));
    FS3min = min(min(comb_data.FS3_conc));


    S4min = min(min(comb_data.S4_conc));
    BS4min = min(min(comb_data.BS4_conc));
    FS4min = min(min(comb_data.FS4_conc));
%--------------------------------------------------------------------------           

    Bmax = max(max(comb_data.B_conc));
    Rmax = max(max(comb_data.R_conc));

    S1cmax = max(max(comb_data.S1_conc));
    S1max = max(max(comb_data.S1_conc));
    BS1max = max(max(comb_data.BS1_conc));
    FS1max = max(max(comb_data.FS1_conc));

    S2max = max(max(comb_data.S2_conc));
    BS2max = max(max(comb_data.BS2_conc));
    FS2max = max(max(comb_data.FS2_conc));

    S3cmax = max(max(comb_data.S3c_conc));
    S3max = max(max(comb_data.S3_conc));
    BS3max = max(max(comb_data.BS3_conc));
    FS3max = max(max(comb_data.FS3_conc));

    S4max = max(max(comb_data.S4_conc));
    BS4max = max(max(comb_data.BS4_conc));
    FS4max = max(max(comb_data.FS4_conc));
    
%%
        tspan = 1:1:size(comb_data.F_molecule,2);  

 
        P = param.P;
        T = param.T;
        l = 400;
        x = linspace(param.x1,param.x2,l);
        y = linspace(param.y1,param.y2,l);
        
   


        for t = 1:10:tspan(end) % every 10th frame
                            
            fprintf('%f\n',t)
                            
 
  %--------------------------------------------------------------------------            
           
            Ft = squeeze(comb_data.F_conc(:,t));
            Bt = squeeze(comb_data.B_conc(:,t));
            Rt = squeeze(comb_data.R_conc(:,t));

            S1ct = squeeze(comb_data.S1c_conc(:,t));
            S1t = squeeze(comb_data.S1_conc(:,t));
            BS1t = squeeze(comb_data.BS1_conc(:,t));
            FS1t = squeeze(comb_data.FS1_conc(:,t));

            S2t = squeeze(comb_data.S2_conc(:,t));
            BS2t = squeeze(comb_data.BS2_conc(:,t));
            FS2t = squeeze(comb_data.FS2_conc(:,t));

            S3ct = squeeze(comb_data.S3c_conc(:,t));
            S3t = squeeze(comb_data.S3_conc(:,t));
            BS3t = squeeze(comb_data.BS3_conc(:,t));
            FS3t = squeeze(comb_data.FS3_conc(:,t));

            S4t = squeeze(comb_data.S4_conc(:,t));
            BS4t = squeeze(comb_data.BS4_conc(:,t));
            FS4t = squeeze(comb_data.FS4_conc(:,t));

            S1_Tt = S1t + BS1t + FS1t;
            S2_Tt = S2t + BS2t + FS2t;
            S3_Tt = S3t + BS3t + FS3t;
            S4_Tt = S4t + BS4t + FS4t;
%--------------------------------------------------------------------------            
            [rgbFront,tn,ta2,ta3] = tri2grid(P,T,Ft,x,y);           
            rgbFront = (rgbFront-F_low)/(F_high-F_low);

                                           
            [rgbBack,  tn,ta2,ta3] = tri2grid(P,T,Bt,tn,ta2,ta3);
            rgbBack = (rgbBack-Bmin)/(Bmax-Bmin);

            [rgbRef,  tn,ta2,ta3] = tri2grid(P,T,Rt,tn,ta2,ta3); 
            rgbRef  = (rgbRef-Rmin)/(Rmax-Rmin);
             
            [rgbS1c ,tn,ta2,ta3] = tri2grid(P,T,S1ct,tn,ta2,ta3);
            rgbS1c = (rgbS1c-S1cmin)/(S1cmax-S1cmin);

            [rgbS1 ,tn,ta2,ta3] = tri2grid(P,T,S1t,tn,ta2,ta3);
            rgbS1 = (rgbS1-S1min)/(S1max-S1min);


            [rgbBS1  ,tn,ta2,ta3] = tri2grid(P,T,BS1t,tn,ta2,ta3);        
            rgbBS1  = (rgbBS1-BS1min)/(BS1max-BS1min);

            [rgbFS1  ,tn,ta2,ta3] = tri2grid(P,T,FS1t,tn,ta2,ta3);        
            rgbFS1  = (rgbFS1-FS1min)/(FS1max-FS1min);
            

            [rgbS2 ,tn,ta2,ta3] = tri2grid(P,T,S2t,tn,ta2,ta3);
            rgbS2 = (rgbS2-S2min)/(S2max-S2min);

            [rgbBS2  ,tn,ta2,ta3] = tri2grid(P,T,BS2t,tn,ta2,ta3);
            rgbBS2  = (rgbBS2-BS2min)/(BS2max-BS2min);

            [rgbFS2  ,tn,ta2,ta3] = tri2grid(P,T,FS2t,tn,ta2,ta3);
            rgbFS2  = (rgbFS2-FS2min)/(FS2max-FS2min);

            [rgbS3c ,tn,ta2,ta3] = tri2grid(P,T,S3ct,tn,ta2,ta3);
            rgbS3c = (rgbS3c-S3cmin)/(S3cmax-S3cmin);

            [rgbS3 ,tn,ta2,ta3] = tri2grid(P,T,S3t,tn,ta2,ta3);
            rgbS3 = (rgbS3-S3min)/(S3max-S3min);


            [rgbBS3  ,tn,ta2,ta3] = tri2grid(P,T,BS3t,tn,ta2,ta3);        
            rgbBS3  = (rgbBS3-BS3min)/(BS3max-BS3min);

            [rgbFS3  ,tn,ta2,ta3] = tri2grid(P,T,FS3t,tn,ta2,ta3);        
            rgbFS3  = (rgbFS3-FS3min)/(FS3max-FS3min);
            

            [rgbS4 ,tn,ta2,ta3] = tri2grid(P,T,S4t,tn,ta2,ta3);
            rgbS4 = (rgbS4-S4min)/(S4max-S4min);

            [rgbBS4  ,tn,ta2,ta3] = tri2grid(P,T,BS4t,tn,ta2,ta3);
            rgbBS4  = (rgbBS4-BS4min)/(BS4max-BS4min);

            [rgbFS4  ,tn,ta2,ta3] = tri2grid(P,T,FS4t,tn,ta2,ta3);
            rgbFS4  = (rgbFS4-FS4min)/(FS4max-FS4min);
            
            
            [rgbS1_T ,tn,ta2,ta3] = tri2grid(P,T,S1_Tt,tn,ta2,ta3);           
            [rgbS2_T ,tn,ta2,ta3] = tri2grid(P,T,S2_Tt,tn,ta2,ta3);
            [rgbS3_T ,tn,ta2,ta3] = tri2grid(P,T,S3_Tt,tn,ta2,ta3);            
            [rgbS4_T ,tn,ta2,ta3] = tri2grid(P,T,S4_Tt,tn,ta2,ta3);

            rgbS1_T = (rgbS1_T-S1min)/(S1max-S1min);
            rgbS2_T = (rgbS2_T-S2min)/(S2max-S2min);
            rgbS3_T = (rgbS3_T-S3min)/(S3max-S3min);
            rgbS4_T = (rgbS4_T-S4min)/(S4max-S4min);
%--------------------------------------------------------------------------   

            image_F = imadjust(rgbFront, [0 1]);
            image_B = imadjust(rgbBack, [0 1]);
            image_R = imadjust(rgbRef, [0 1]);
      
            image_S1c = imadjust(rgbS1c, [0 1]);
            image_S1  = imadjust(rgbS1, [0 1]);
            image_BS1 = imadjust(rgbBS1, [0 1]);
            image_FS1 = imadjust(rgbFS1, [0 1]);
            image_S1T = imadjust(rgbS1_T, [0 1]);

            
            image_S2 = imadjust(rgbS2, [0 1]);
            image_BS2 = imadjust(rgbBS2, [0 1]);
            image_FS2 = imadjust(rgbFS2, [0 1]);
            image_S2T = imadjust(rgbS2_T, [0 1]);


            image_S3T = min(1,max(0,rgbS3_T));
            image_S4T = min(1,max(0,rgbS4_T));
%--------------------------------------------------------------------------            
% Uncomment for saving the data in .tiff format   



        imwrite(image_F, fullfile(Data_folder,'F.tif'),'WriteMode','append','Compression','none');
        imwrite(image_B, fullfile(Data_folder,'B.tif'),'WriteMode','append','Compression','none');
        imwrite(image_R, fullfile(Data_folder,'R.tif'),'WriteMode','append','Compression','none');
        imwrite(image_S1c, fullfile(Data_folder,'S1c.tif'),'WriteMode','append','Compression','none');
        imwrite(image_S1T, fullfile(Data_folder,'S1T.tif'),'WriteMode','append','Compression','none');
        imwrite(image_S2T, fullfile(Data_folder,'S2T.tif'),'WriteMode','append','Compression','none');
        imwrite(image_S3T, fullfile(Data_folder,'S3T.tif'),'WriteMode','append','Compression','none');
        imwrite(image_S4T, fullfile(Data_folder,'S4T.tif'),'WriteMode','append','Compression','none');
%--------------------------------------------------------------------------            
        end
end

function [Node_index_front, Node_index_back] = front_back(comb_data,param)

    Bmin = min(min(comb_data.B_conc));
    Fmin = min(min(comb_data.F_conc));
    Bmax = max(max(comb_data.B_conc));
    Fmax = max(max(comb_data.F_conc));

    l = 200;
    x = linspace(param.x1,param.x2,l);
    y = linspace(param.y1,param.y2,l);
    [xx,yy] = meshgrid(x);

    N = (length(comb_data.tspan)-1)*param.dt;
    Node_index_front = cell(N,1);
    Node_index_back = cell(N,1);
    
    
    tic
    parfor t = 1:N+1 % every 10th frame
        idx = 1+(t-1)/param.dt;
        Bt = squeeze(comb_data.B_conc(:,idx));
        [rgbBack,~,~,~] = tri2grid(param.P,param.T,Bt,x,y);
        rgbBack = (rgbBack-Bmin)/(Bmax-Bmin);

        Ft = squeeze(comb_data.F_conc(:,idx));
        [rgbFront,~,~,~] = tri2grid(param.P,param.T,Ft,x,y);
        rgbFront = (rgbFront-Fmin)/(Fmax-Fmin);
%--------------------------------------------------------------------------      

        Front_tmp = 0*rgbFront;
        Front_tmp(rgbFront>0.2) = 1;
        Front_tmp = imfill(Front_tmp,'holes');
        Front_tmp = bwareaopen(Front_tmp,5);
        index_front = find(Front_tmp == 1);
   
        x_front = xx(index_front);
        y_front = yy(index_front);
                        
        Back_tmp = 0*rgbBack;
        Back_tmp(rgbBack>0.2) = 1;
        Back_tmp = imfill(Back_tmp,'holes');
        Back_tmp = bwareaopen(Back_tmp,5);
        index_back = find(Back_tmp == 1);
   
        x_back = xx(index_back);
        y_back = yy(index_back);


%--------------------------------------------------------------------------

%--------------------------------------------------------------------------    
        p_f = [x_front,20-y_front];
        Node_index_front_tmp = dsearchn(param.P',p_f);
        Node_index_front{t} = unique(Node_index_front_tmp);

        p_b = [x_back,20-y_back];
        Node_index_back_tmp = dsearchn(param.P',p_b);
        Node_index_back{t} = unique(Node_index_back_tmp);
%--------------------------------------------------------------------------

   
    end
    elapsed_time = toc;
    fprintf('Elapsed time = %f in min \n',elapsed_time/60)
end