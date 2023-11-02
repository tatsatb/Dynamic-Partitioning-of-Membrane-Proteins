function [mask_r_I_avg, mask_r_I_tot,...
          cell_Ig_avg_wrt_bg, cell_Ig_tot_wrt_bg,...
          cell_Ir_avg_wrt_bg, cell_Ir_tot_wrt_bg] = intensity_calculation(mask_r, mask_cell, mask_bg, filt_image_g, filt_image_r)
      
        if sum(mask_r(:)) == 0
            mask_r_I_tot = NaN;
        else
            mask_r_I_tot_tmp = filt_image_r.*uint16(mask_r);
            mask_r_I_tot = sum(mask_r_I_tot_tmp(:));
        end
        mask_r_I_avg = mask_r_I_tot/sum(mask_r(:));

        cell_Ir_tot_tmp = filt_image_r.*uint16(mask_cell);
        cell_Ir_tot = sum(cell_Ir_tot_tmp(:));      
        cell_Ir_avg = cell_Ir_tot/sum(mask_cell(:));

        cell_Ig_tot_tmp = filt_image_g.*uint16(mask_cell);
        cell_Ig_tot = sum(cell_Ig_tot_tmp(:));      
        cell_Ig_avg = cell_Ig_tot/sum(mask_cell(:));

        cell_Ir_bg_tot_tmp = filt_image_r.*uint16(mask_bg);
        cell_Ir_bg_tot = sum(cell_Ir_bg_tot_tmp(:));      
        cell_Ir_bg_avg = cell_Ir_bg_tot/sum(mask_bg(:));

        cell_Ig_bg_tot_tmp = filt_image_g.*uint16(mask_bg);
        cell_Ig_bg_tot = sum(cell_Ig_bg_tot_tmp(:));      
        cell_Ig_bg_avg = cell_Ig_bg_tot/sum(mask_bg(:));


        cell_Ig_tot_wrt_bg = cell_Ig_tot - cell_Ig_bg_tot;
        cell_Ir_tot_wrt_bg = cell_Ir_tot - cell_Ir_bg_tot;
        cell_Ig_avg_wrt_bg = cell_Ig_avg - cell_Ig_bg_avg;
        cell_Ir_avg_wrt_bg = cell_Ir_avg - cell_Ir_bg_avg;
        
end