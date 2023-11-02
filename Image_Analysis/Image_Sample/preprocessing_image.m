function [filt_image,I] = preprocessing_image(image_frame, thresh, initial_cell_mask)
        
        
        % Modify the parameters to get proper segmentation
        
        % Solving nonuniform background and background subtraction  
          bg = imopen(image_frame,strel('disk',150));
          image1 = image_frame - bg;

        % morphological top-hat filtering 
          tophatFiltered_cell = imtophat(image1, strel('disk',100));

        % gaussian smoothing, if required
          
          %filt_image = imgaussfilt(tophatFiltered_cell,1);
          filt_image = tophatFiltered_cell;             
            
        I  = zeros(size(filt_image,1),size(filt_image,2));
        I(filt_image.*uint16(initial_cell_mask)>thresh) = 1;

end