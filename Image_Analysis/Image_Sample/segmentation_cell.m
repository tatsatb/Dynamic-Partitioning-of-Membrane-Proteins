function [mask_cell, image_mod, mask_bg] = segmentation_cell(I, frame_index)

        
%%      Processing by dilation, erosion, filling holes and eliminating isolated components

        image_mod = bwareaopen(I,20); % 20
%       image_mod = imclearborder(image_mod);
        image_mod = imfill(image_mod,'holes');
                
        for j = 1:2
            image_mod = imdilate(image_mod, strel('disk',1,4));
        end
        
        image_mod = imfill(image_mod,'holes');
         
    
%       To erode the extras
        for k = 1:2
            image_mod = imerode(image_mod,strel('disk',1));
        end
 
%       To fill back 
%        for m = 1:3
%           image_mod = imdilate(image_mod, strel('disk',1,4));
%        end    
        

%%      Extracting the boundary

        stats = regionprops(image_mod,'Area','BoundingBox','Centroid','ConvexHull','Orientation','Perimeter');
    
        area  = zeros(1,length(stats));
        for stat_i = 1:length(stats)
            area(stat_i) = stats(stat_i).Area;
        end

        [~,index1] = max(area);
        BB = stats(index1).BoundingBox + [-1 -1 2 2];
        TMP_BW = imcrop(image_mod,BB);

        [B,~] = bwboundaries(TMP_BW,'noholes');

        area_enclosed = zeros(1,length(B));
        for stat_j = 1:length(B)
            TMP_boundary = B{stat_j};  
            area_enclosed(stat_j) = polyarea(TMP_boundary(:,2),TMP_boundary(:,1));
        end
        [~,index2] = max(area_enclosed);
        TMP_boundary = B{index2};    
    
        % mask for the bright green region
        mask_cell = poly2mask(TMP_boundary(:,2)+BB(1),TMP_boundary(:,1)+BB(2),size(I,1),size(I,2));
        boundary_abs = [TMP_boundary(:,2)+BB(1),TMP_boundary(:,1)+BB(2)];


        % Extracting background around the cell
        convex_hull_tmp = {stats.ConvexHull};
        convex_hull = convex_hull_tmp{index1};
    
        additional_L = 20;
        cell_bg_dim = [BB(1)-additional_L, BB(2)-additional_L, BB(3)+2*additional_L, BB(4)+2*additional_L];
    
        backgroud_mask_tmp = poly2mask( [cell_bg_dim(1),cell_bg_dim(1)+cell_bg_dim(3),cell_bg_dim(1)+cell_bg_dim(3),cell_bg_dim(1),cell_bg_dim(1)],...
                                        [cell_bg_dim(2),cell_bg_dim(2),cell_bg_dim(2)+cell_bg_dim(4),cell_bg_dim(2)+cell_bg_dim(4),cell_bg_dim(2)],...
                                        size(I,1),size(I,2));
                                
        convex_mask = poly2mask(convex_hull(:,1),convex_hull(:,2),size(I,1),size(I,2));
                          

    
       image_mod_2 = bwareaopen(I,20); % 20
%      image_mod_2 = imclearborder(image_mod_2);
       image_mod_2 = imfill(image_mod_2,'holes');
        
       

        for j = 1
            image_mod_2 = imdilate(image_mod_2, strel('disk',1,4));
        end
        image_mod_2 = imfill(image_mod_2,'holes');
        
%       To erode the extras
%       for k = 1:2
%           I_tmp2 = imerode(I_tmp2,strel('disk',1));
%       end
     

        stats_2 = regionprops(image_mod_2,'Area','BoundingBox','Centroid','ConvexHull','Orientation','Perimeter');
    
        
        convex_hull_2_tmp = {stats_2.ConvexHull};
        
        convex_mask_comb = zeros(size(image_mod_2,1),size(image_mod_2,2));
        
        
        for i = 1:length(convex_hull_2_tmp)
            
            convex_hull = convex_hull_2_tmp{i};
            convex_mask_tmp = poly2mask(convex_hull(:,1),convex_hull(:,2),size(image_mod,1),size(image_mod,2));
            convex_mask_comb = or(convex_mask_comb,convex_mask_tmp);
        end

        mask_bg = and(backgroud_mask_tmp,~convex_mask_comb); 
end