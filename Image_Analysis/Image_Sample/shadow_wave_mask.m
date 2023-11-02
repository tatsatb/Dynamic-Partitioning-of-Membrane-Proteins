function [mask_SW, mask, sw_mask_near_comb, sw_mask_near, boundary_abs] = shadow_wave_mask(mask_cell,mask_back,mask_circle)

    mask_sw = (mask_cell.*~mask_back);
    
    
    mask_sw = bwareaopen(mask_sw,20);
    
    
    % To erode the extras
    for k = 2
    mask_sw = imerode(mask_sw,strel('disk',1));
    end

    mask_sw = bwareaopen(mask_sw,40); % 20
    
    
    %%  Extracting the boundary and segmenting the shadow wave mask

        stats = regionprops(mask_sw,'Area','BoundingBox','Centroid','ConvexHull','Orientation','Perimeter');
    
        mask_SW = zeros(size(mask_sw,1),size(mask_sw,2));
        
        
        area  = zeros(1,length(stats));
        c = 1;
        mask = cell(1,10);
        mask_holes = cell(1,10);
        boundary_abs = cell(1,10);
        for stat_i = 1:length(stats)
            
            area(stat_i) = stats(stat_i).Area;
            
            if area(stat_i) > 150
                
                BB = stats(stat_i).BoundingBox + [-1 -1 2 2];
                TMP_BW = imcrop(mask_sw,BB);
                
                [B,L] = bwboundaries(TMP_BW,'noholes');
                    
                [B_holes,L_holes] = bwboundaries(TMP_BW,'holes');
                
                
                if isequal(B,B_holes) ~= 1 && length(find(L_holes==2))>500
                    
                        
                        holes_binary=im2bw(L_holes,1);
                        
                        [BX,~] = bwboundaries(holes_binary,'noholes');
                        
                         
                        area_enclosed_holes = zeros(1,length(BX));
                        
                        
                        for stat_holes = 1:length(BX)
                            TMP_boundary_holes = BX{stat_holes};
                                                        
                            area_enclosed_holes(stat_holes) = polyarea(TMP_boundary_holes(:,2),TMP_boundary_holes(:,1));
                    
                        end
                        
                        
                        [~,index2_holes] = max(area_enclosed_holes);
                        
                        TMP_boundary_holes = BX{index2_holes}; 
                        
                        mask_holes{c}=poly2mask(TMP_boundary_holes(:,2)+BB(1),TMP_boundary_holes(:,1)+BB(2),size(mask_sw,1),size(mask_sw,2));
                        


                else
                    
                    mask_holes{c}=0;
                        
                end
                

                
            
                
                               
                area_enclosed = zeros(1,length(B));
                
                
                for stat_j = 1:length(B)
                    TMP_boundary = B{stat_j};  
                    
                    
                    area_enclosed(stat_j) = polyarea(TMP_boundary(:,2),TMP_boundary(:,1));
                    
                    
                end
                [~,index2] = max(area_enclosed);
                TMP_boundary = B{index2};    

                
                mask{c} = poly2mask(TMP_boundary(:,2)+BB(1),TMP_boundary(:,1)+BB(2),size(mask_sw,1),size(mask_sw,2));
                

                boundary_abs{c} = [TMP_boundary(:,2)+BB(1),TMP_boundary(:,1)+BB(2)];
                
                if length(find(mask{c} - mask_holes{c}))<1000
                
                    mask_holes{c}=0;
                    
                end
                
                
                mask{c} = mask{c} - mask_holes{c}; 
                
            
                
                mask_SW = or(mask_SW,mask{c});
                c = c+1;
            end
                
        end


        mask = 0;  
        sw_mask_near = 0;        
        sw_mask_near_comb = and(mask_SW,mask_circle);
        
        
end