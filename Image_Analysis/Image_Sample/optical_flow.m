%----------------------------------------------------------------------------
% Intensity and optical flow analysis of photoconvertion microscopy
% based protein tracking assay
%
%----------------------------------------------------------------------------
% This code was developed to generate Figure 4, Supplementary Figure 9, and
% Supplementary figure 10 of the following paper:
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



clear
close all
clc

format short
warning('off','all'); 
  
  

Path_to_Load = pwd;
file_name = 'Image_Sample';

full_path_to_load = fullfile(Path_to_Load,file_name);

idcs   = strfind(Path_to_Load,'\');
newdir = Path_to_Load(1:idcs(end)-1);


% Define the root folder to save different results 
Path_to_Save = strcat(newdir,'\Save_Folder'); 
     
%% Select initial area of interest or use an existing area of interest

% Change to 0, to select new areas of interest (masks) to initialize segmentaion. 
% Keep 1, if you already have existing masks for the cell and background (and you want to use those).

user_input_temp=1;    

     
%% Save the channel data separately to mat files for quicker access    

savedata_as_MATRIX(full_path_to_load,Path_to_Save);


     
%%  Load separate channels as a matrix 

    full_path_green_chanel_data = strcat(fullfile(Path_to_Save,file_name),'\Green_chanel.mat');
    full_path_red_chanel_data = strcat(fullfile(Path_to_Save,file_name),'\Red_chanel.mat');

    n_frames = 16; % # of frames upto which the analysis has to be perfomed
    
    %  loading the data
    D_g  = load(full_path_green_chanel_data);
    Green_chanel = D_g.Green_chanel;

    D_r  = load(full_path_red_chanel_data);
    Red_chanel = D_r.Red_chanel;

%%   Create required directories


    full_Path_to_Save = [Path_to_Save,'\',file_name];
    
    
    if ~exist(fullfile(full_Path_to_Save, 'quiver_figures'),'dir') == 1
        mkdir(fullfile(full_Path_to_Save, 'quiver_figures'))
    end
    
    if ~exist(fullfile(full_Path_to_Save, 'mask_figures'),'dir') == 1
    mkdir(fullfile(full_Path_to_Save, 'mask_figures'))
    end

  
%%  Initialize optical flow objects
  
  opticFlow_mask_r = opticalFlowHS; % Estimate optical flow of red photoconverted regions using Horn-Schunck method
  opticFlow_mask_sw = opticalFlowHS; % Estimate optical flow of shadow waves using Horn-Schunck method

  
%%   Iterate over required frames
  
    for frame_index = 2:1:n_frames

        fprintf('frame = %d\n',frame_index)  
        image_frame_green = Green_chanel(:,:,frame_index);
        image_frame_red = Red_chanel(:,:,frame_index);
        
        
     
%-------------------------------------------------------------------------------------------------------------------
%% Generate and save or load the regions of interest (masks) for the cell and background 

% If the code is being executed for the first time (i.e. no predefinded 
% area of interest exists) please use the following section
% code to generate it (first input should be set as: user_input_temp == 0). 

% Select or load pre-selected areas of interest for the cell

       if user_input_temp == 0
        initial_cell_mask = user_defined_cell_mask(imtophat(image_frame_green, strel('disk',100)));

        
        while any(initial_cell_mask(:)) == 0 || logical(sum(initial_cell_mask(:))<20000)
            initial_cell_mask = user_defined_cell_mask(imtophat(image_frame_green, strel('disk',100)));
        end 
               
       
        MASK.initial_cell_mask_full{frame_index}=initial_cell_mask;
        save([full_Path_to_Save, '\MASK.mat'],'MASK');
            
        
        
%       Load the predefined maks from the directory 
       else
              
            load([full_Path_to_Save, '\MASK.mat']); 
            initial_cell_mask= MASK.initial_cell_mask_full{frame_index};
            
       end


       
% Select or load pre-selected areas of interest for the background

       if user_input_temp == 0
            mask_bg = user_defined_bg(filt_image_g);
        
        while any(mask_bg(:)) == 0 || logical(sum(mask_bg(:))>20000)
            mask_bg = user_defined_bg(filt_image_g);
        end
           
        
%       Save masks of background        
        
        MASK.mask_bg_full{frame_index}=mask_bg;
        save([full_Path_to_Save, '\MASK.mat'],'MASK');

%       Load masks of background  
           
       else 
           mask_bg= MASK.mask_bg_full{frame_index};
       end 
%--------------------------------------------------------------------------------------------------------
       

%% Segmentations of the back-state of the cell, photoconverted region, and the whole cell
        
        % Segmentation of back-state regions of the cell (from green channel) 
        
        G_thresh_1 = 1000; 
        [filt_image_g,I_g1] = preprocessing_image(image_frame_green, G_thresh_1,initial_cell_mask);

        [mask_g,boundary_abs_g, image_mod, mask] = segmentation_green(I_g1, frame_index);
               
        
        
        % Segmentation of photoconversion area of the cell (from red channel)
        
        R_thresh = 1100; 
        [filt_image_r,I_r] = preprocessing_image(image_frame_red, R_thresh,initial_cell_mask);
        [mask_r,boundary_abs_r,area_r,centroid_mask_r] = segmentation_red(I_r, frame_index); 
        
        
        % Segmenting whole cell (from green channel)
        
        G_thresh_2 = 400;
        [~,I_g2] = preprocessing_image(image_frame_green, G_thresh_2,initial_cell_mask);
        [mask_cell, image_mod, ~] = segmentation_cell(I_g2, frame_index);
        
        
%         imshow(mask_r)        
%         title(['Segmented red mask, Frame - ' num2str(frame_index)])
%         hold on
%                
%         xcentroid=centroid_mask_r(:,1);
%         ycentroid=centroid_mask_r(:,2);
%         plot(xcentroid, ycentroid, 'r+', 'MarkerSize', 15, 'LineWidth', 2);       
%         
%         hold off
        

%% Find the minimum bounding circle of the photoconverted region

        if sum(mask_r(:)) >0
            
            [center,radius] = minboundcircle(boundary_abs_r(:,1),boundary_abs_r(:,2));
            angle = linspace(0, 2*pi, 1000);
            x = radius * cos(angle) + center(1);
            y = radius * sin(angle) + center(2);
            
            
            L = max(0.2*radius,20);
            x_prime = (radius+L) * cos(angle) + center(1);
            y_prime = (radius+L) * sin(angle) + center(2);
            
            
            
            mask_circle = poly2mask(x_prime,y_prime,size(mask_r, 1), size(mask_r, 2));


        else

             x = NaN;
             y = NaN;
             x_prime = NaN;
             y_prime = NaN;
             
             mask_circle = zeros(size(mask_r, 1), size(mask_r, 2));
        end
            
%         imshow(mask_circle) 
%         title(['\textbf{mask circle, Frame -} ' num2str(frame_index)],'interpreter','Latex')



     % Total back-state regions of the cell = mask_g + mask_r
     
        mask_back = or(mask_g,mask_r);
        Back_area = sum(mask_back(:));
        
     % Shadow Wave = mask_cell - mask_back  
     
        [mask_sw, ~, sw_mask_near_comb, ~, boundary_abs_sw] = shadow_wave_mask(mask_cell,mask_back,mask_circle);
        
%       imshowpair(mask_back,mask_sw,'montage')
%       title(['\textbf{Shadow wave and back region masks -} ' num2str(frame_index)],'interpreter','Latex')

        
        
        boundary_abs_sw(cellfun('isempty',boundary_abs_sw)) = []; % remove the empty cells from a vector of cells
        
%% Calculate intensities          
        
        [mask_r_I_avg, mask_r_I_tot,...
         cell_Ig_avg_wrt_bg, cell_Ig_tot_wrt_bg,...
         cell_Ir_avg_wrt_bg, cell_Ir_tot_wrt_bg] = intensity_calculation(mask_r, mask_cell, mask_bg, filt_image_g, filt_image_r);
     

     
%% Compute the optical flow vectors and save the quiver figures


          red_mask_flow = estimateFlow(opticFlow_mask_r, logical(mask_r));
          sw_flow = estimateFlow(opticFlow_mask_sw, logical(sw_mask_near_comb));     
     
          
          if sum(mask_r(:)) >0
              [eroded_mask_r, centroid_eroded_mask_r] = mask_erosion(mask_r);
          else
              eroded_mask_r = zeros(size(mask_r,1),size(mask_r,2));
              centroid_eroded_mask_r = [NaN, NaN];
          end
        
          
          
            h = figure;
            movegui(h);
            hViewPanel = uipanel(h,'Position',[0 0 1 1],'Title','Plot of Optical Flow Vectors');
            hPlot = axes(hViewPanel);
         
            

              
          if frame_index==3 % The first frame of photoconversion 

              imshowpair(mask_r,sw_mask_near_comb,'ColorChannels',[1 0 2])
              hold on

              newcolors = [1 1 1
                  0 1 0];
              
              colororder(newcolors)
              
              plot(red_mask_flow,'DecimationFactor',[50 50],'ScaleFactor',1,'Parent',hPlot);

%               q_red_mask_flow = findobj(gca,'type','Quiver');
%               q_red_mask_flow.Color = 'w';
               

              plot(sw_flow,'DecimationFactor',[8 8],'ScaleFactor',60,'Parent',hPlot);

              hold off
          
   
              
          else % Subsequent frames of moving photoconvered regions
              imshowpair(mask_r,sw_mask_near_comb,'ColorChannels',[1 0 2])
              hold on
                
              newcolors = [1 1 1
                  0 1 0];
         
              colororder(newcolors)
              
              plot(red_mask_flow,'DecimationFactor',[8 8],'ScaleFactor',60,'Parent',hPlot);


              plot(sw_flow,'DecimationFactor',[8 8],'ScaleFactor',60,'Parent',hPlot);

              hold off
             
                    
          end
          
          

        saveas(gcf,fullfile(full_Path_to_Save, 'quiver_figures',strcat('frame_',num2str(frame_index),'.png')))
        saveas(gcf,fullfile(full_Path_to_Save, 'quiver_figures',strcat('frame_',num2str(frame_index),'.fig')))  
          
          
%% Remaining calculations and saving the figures with combined segmentation

          Resultant_Vx_mask_r = sum(sum(red_mask_flow.Vx.*single(eroded_mask_r)));
          Resultant_Vy_mask_r = sum(sum(red_mask_flow.Vy.*single(eroded_mask_r)));
          
          R_magnitude_mask_r = sqrt(Resultant_Vx_mask_r^2+Resultant_Vy_mask_r^2);
          R_angle_mask_r = atan2(Resultant_Vy_mask_r,Resultant_Vx_mask_r);
          
          
          Resultant_Vx_sw_mask_comb = sum(sum(sw_flow.Vx.*single(sw_mask_near_comb)));
          Resultant_Vy_sw_mask_comb = sum(sum(sw_flow.Vy.*single(sw_mask_near_comb)));
          
          R_magnitude_sw_mask_comb = sqrt(Resultant_Vx_sw_mask_comb^2+Resultant_Vy_sw_mask_comb^2);
          R_angle_mask_sw_mask_comb = atan2(Resultant_Vy_sw_mask_comb,Resultant_Vx_sw_mask_comb);
          
          
          mask_r_vector = [R_magnitude_mask_r*cos(R_angle_mask_r), R_magnitude_mask_r*sin(R_angle_mask_r)];
        
          
          mask_sw_vector_double = [R_magnitude_sw_mask_comb*cos(R_angle_mask_sw_mask_comb), R_magnitude_sw_mask_comb*sin(R_angle_mask_sw_mask_comb)];
          
          
          unit_mask_r_vector=mask_r_vector/norm(mask_r_vector);

          
          mask_sw_vector_projected_double = dot(mask_sw_vector_double,mask_r_vector/norm(mask_r_vector)^2)*mask_r_vector;   
        
          
          unit_mask_sw_vector_double=mask_sw_vector_double/norm(mask_sw_vector_double);
          
          dot_product_double=dot(unit_mask_sw_vector_double,unit_mask_r_vector);
          
          
          
          
          
%         figure('color','white')  
%         title(['\textbf{Shadow wave and back region masks -} ' num2str(frame_index)],'interpreter','Latex')

          
        figure('color','white')
        hold on
        imshowpair(mask_back,mask_sw,'ColorChannels',[0 1 2])

        hold on

        plot(boundary_abs_r(:,1), boundary_abs_r(:,2),'k','linewidth',2)

        %---------------------------------------------------------------%
        
        plot(x,y,'m--','linewidth',2)
        plot(x_prime,y_prime,'m--','linewidth',2)

        
        title(['\textbf{Frame $\#$ (combined segmented image): } ' num2str(frame_index)],'interpreter','Latex')

        
        saveas(gcf,fullfile(full_Path_to_Save, 'mask_figures',strcat('frame_',num2str(frame_index),'.png')))
        saveas(gcf,fullfile(full_Path_to_Save, 'mask_figures',strcat('frame_',num2str(frame_index),'.emf')))
            
        %%  Storing in structures        

        
        Data.Mask_r_I_tot{frame_index-1} = mask_r_I_tot;
        Data.Angle_between_SW_R_Vectors_double(frame_index-1) = acosd(dot_product_double);
        
        
    end


save([full_Path_to_Save, '\DATA.mat'],'Data');


%-------------------------------------------------------------------------%


function cell_mask = user_defined_cell_mask(filt_image_g)
% user specific selection of region of interest to help in segmenting the
% cell


        figure('color','white')
        hold on
        title('Draw AROUND cell to select an initial area of interest')
        cell_mask = roipoly(imadjust(filt_image_g));
        close all
end


function mask_bg = user_defined_bg(filt_image_g)
% user specific background selection


        figure('color','white')
        hold on
        title('Draw a BACKGROUND polygon','color','red')
        mask_bg = roipoly(imadjust(filt_image_g));
        
        close all
end
%-------------------------------------------------------------------------% 





  