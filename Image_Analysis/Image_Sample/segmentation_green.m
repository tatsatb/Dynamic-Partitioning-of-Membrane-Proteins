function [MASK_g,boundary_abs_2, image_mod, mask] = segmentation_green(I, frame_index)


%%      Processing by dilation, erosion, filling holes and eliminating isolated components

        image_mod = bwareaopen(I,20); 
%         figure; 
%         imshow(image_mod);
        
        
%       image_mod = imclearborder(image_mod);
%       image_mod = imfill(image_mod,'holes');
%       image_mod = bwareaopen(I,20);
        
        % To connect the lines
        for j = 1:4
            image_mod = imdilate(image_mod, strel('disk',1,4));
%           figure; 
%           imshow(image_mod)
        end

 
    
       %  To erode the extras
       %  for k = 1:2
%             image_mod = imerode(image_mod,strel('disk',1));
%         end



%%      Extracting the boundary

        stats = regionprops(image_mod,'Area','BoundingBox','Centroid','ConvexHull','Orientation','Perimeter');
        
        area  = zeros(1,length(stats));
        
        
        for stat_i = 1:length(stats)
            area(stat_i) = stats(stat_i).Area;
        end

        
        [~,index1] = max(area);
        BB = stats(index1).BoundingBox + [-1 -1 2 2];
        TMP_BW = imcrop(image_mod,BB);
        
        

        [B,~] = bwboundaries(TMP_BW,'holes');


        area_enclosed = zeros(1,length(B));
        
        for stat_j = 1:length(B)
            TMP_boundary = B{stat_j};  
            area_enclosed(stat_j) = polyarea(TMP_boundary(:,2),TMP_boundary(:,1));
        end
        
        
        [~,index2] = max(area_enclosed);
        TMP_boundary = B{index2};    
        
        
        
        % mask for the bright green region
        mask = poly2mask(TMP_boundary(:,2)+BB(1),TMP_boundary(:,1)+BB(2),size(I,1),size(I,2));
        boundary_abs = [TMP_boundary(:,2)+BB(1),TMP_boundary(:,1)+BB(2)];
    
        mask_g = mask;



        image_mod_2 = ~image_mod.*mask_g;
     
        
        for tempo=1:2
        image_mod_2 = imerode(image_mod_2,strel('disk',1,8)); 
        end
        
        
        image_mod_2 = bwareaopen(image_mod_2,250); 
        
%         figure; 
%         imshow(image_mod_2);


%       image_mod_2 = imclearborder(image_mod_2);
        image_mod_2 = imfill(image_mod_2,'holes');
    
        
%      To connect the lines
        for j = 1:3
            image_mod_2 = imdilate(image_mod_2, strel('disk',1,4));
        end
        
        image_mod_2 = imfill(image_mod_2,'holes');
        
%       To erode the extras
%       for k = 1:2
%           I_tmp2 = imerode(I_tmp2,strel('disk',1));
%       end
     

        stats_2 = regionprops(image_mod_2,'Area','BoundingBox','Centroid','ConvexHull','Orientation','Perimeter');
    
        mask_g_2 = zeros(size(image_mod_2,1),size(image_mod_2,2));
        area_2  = zeros(1,length(stats_2));
        c_2 = 1;
        mask_2 = cell(1,10);
        boundary_abs_2 = cell(1,10);
        
        
        for stat_i2 = 1:length(stats_2)
            
            area_2(stat_i2) = stats_2(stat_i2).Area;
            
            if area_2(stat_i2) > 300
                
                BB_2 = stats_2(stat_i2).BoundingBox + [-1 -1 2 2];
                TMP_BW_2 = imcrop(image_mod_2,BB_2);

                [B_2,~] = bwboundaries(TMP_BW_2);

                area_enclosed_2 = zeros(1,length(B_2));
                
                for stat_j = 1:length(B_2)
                    TMP_boundary_2 = B_2{stat_j};  
                    area_enclosed_2(stat_j) = polyarea(TMP_boundary_2(:,2),TMP_boundary_2(:,1));
                end
                [~,index2_2] = max(area_enclosed_2);
                TMP_boundary_2 = B_2{index2_2};    
    
               
                mask_2{c_2} = poly2mask(TMP_boundary_2(:,2)+BB_2(1),TMP_boundary_2(:,1)+BB_2(2),size(image_mod_2,1),size(image_mod_2,2));
                boundary_abs_2{c_2} = [TMP_boundary_2(:,2)+BB_2(1),TMP_boundary_2(:,1)+BB_2(2)];
                
                mask_g_2 = or(mask_g_2,mask_2{c_2});
                c_2 = c_2+1;
            end
                
        end     
         
%---------------------------------------------------------------------------------------------------------------------------------------    

        boundary_abs_2(cellfun('isempty',boundary_abs_2)) = []; % remove the empty cells from a vector of cells
        
        if iscell(boundary_abs) == 1
            boundary_abs(cellfun('isempty',boundary_abs)) = []; % remove the empty cells from a vector of cells
            
            for j = 1:length(boundary_abs)
            
                boundary_abs_2{end+1} = boundary_abs{j};
            
            end

        else    
            boundary_abs_2{end+1} = boundary_abs;
        end

      
%---------------------------------------------------------------------------------------------------------------------------------------    

            MASK_g = mask_g.*~mask_g_2;
            
%           figure; 
%           imshow(MASK_g);
            
            
            for temp=1:6
            MASK_g = imerode(MASK_g,strel('disk',1,8)); 
%           figure; 
%           imshow(MASK_g);
            end


end