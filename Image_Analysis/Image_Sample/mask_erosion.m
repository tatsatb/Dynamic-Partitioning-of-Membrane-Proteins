function [eroded_mask, centroid_eroded_mask_r] = mask_erosion(mask)


        eroded_mask_tmp = mask;
        
        for k = 1
            eroded_mask_tmp = imerode(eroded_mask_tmp,strel('disk',1));
        end
        
        eroded_mask_tmp = bwareaopen(eroded_mask_tmp,20);
        
        
        stats = regionprops(eroded_mask_tmp,'Area','BoundingBox','Centroid','ConvexHull','Orientation','Perimeter');
    
        area  = zeros(1,length(stats));
        for stat_i = 1:length(stats)
            area(stat_i) = stats(stat_i).Area;
        end

        [~,index1] = max(area);
        BB = stats(index1).BoundingBox + [-1 -1 2 2];
        TMP_BW = imcrop(eroded_mask_tmp, BB);

        [B,~] = bwboundaries(TMP_BW,'noholes');

        area_enclosed = zeros(1,length(B));
        for stat_j = 1:length(B)
            TMP_boundary = B{stat_j};  
            area_enclosed(stat_j) = polyarea(TMP_boundary(:,2),TMP_boundary(:,1));
        end
        [~,index2] = max(area_enclosed);
        TMP_boundary = B{index2};    
    

        eroded_mask  = poly2mask(TMP_boundary(:,2)+BB(1),TMP_boundary(:,1)+BB(2),size(mask,1),size(mask,2));
        centroid_eroded_mask_r = stats(index1).Centroid;
        
end