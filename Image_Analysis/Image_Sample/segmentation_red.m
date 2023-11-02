function [mask,boundary_abs,area_r,centroid_mask_r] = segmentation_red(I, frame_index)
        
%%      Processing by dilation, erosion, filling holes and eliminating isolated components

        image_mod = bwareaopen(I,20); % 20
        image_mod = imclearborder(image_mod);
        image_mod = imfill(image_mod,'holes');
              

        % To connect the lines
        for j = 1:2
            image_mod = imdilate(image_mod, strel('disk',1,4));
        end
        
        image_mod = imfill(image_mod,'holes');
        
 
    
%         To erode the extras
%         for k = 1
%             image_mod = imerode(image_mod,strel('disk',1));
%         end
        

 

%%      Extracting the boundary

        stats = regionprops(image_mod,'Area','BoundingBox','Centroid','ConvexHull','Orientation','Perimeter');
    
        area  = zeros(1,length(stats));
        for stat_i = 1:length(stats)
            area(stat_i) = stats(stat_i).Area;
        end

        [max_area,index1] = max(area);
        
        if max_area > 250
            
            BB = stats(index1).BoundingBox + [-1 -1 2 2];
            TMP_BW = imcrop(image_mod,BB);

            [B,~] = bwboundaries(TMP_BW,'noholes');

            area_enclosed = zeros(1,length(B));
            for stat_j = 1:length(B)
                TMP_boundary = B{stat_j};  
                area_enclosed(stat_j) = polyarea(TMP_boundary(:,2),TMP_boundary(:,1));
            end
            [area_r,index2] = max(area_enclosed);
            
            if any(area_r) == 0
                area_r = NaN;
            end
            
            TMP_boundary = B{index2};    
    
            % mask for the bright red region
            mask = poly2mask(TMP_boundary(:,2)+BB(1),TMP_boundary(:,1)+BB(2),size(I,1),size(I,2));
            boundary_abs = [TMP_boundary(:,2)+BB(1),TMP_boundary(:,1)+BB(2)];
            centroid_mask_r = stats(index1).Centroid;
            
        else 
            mask = zeros(size(I,1),size(I,2));
            TMP_boundary = NaN(1,2);
            BB = NaN(1,4);
            boundary_abs = [TMP_boundary(:,2)+BB(1),TMP_boundary(:,1)+BB(2)];
            area_r = NaN;
            centroid_mask_r = [NaN, NaN];
        end


end