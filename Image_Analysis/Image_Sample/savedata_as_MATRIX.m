function savedata_as_MATRIX(full_path_to_load, Path_to_Save)

% This function saves each channel of the multichannel image data 
% as separate matrices. 
 

  index_backslash = strfind(full_path_to_load,'\');
  file_name = full_path_to_load(index_backslash(end)+1:end);
  
  
  full_Path_to_Save = [Path_to_Save,'\',file_name];
  
  if exist('fname','dir') == 0
     mkdir(full_Path_to_Save);
  end
  

  tiffile = strcat(full_path_to_load,'.tif');
     

     
 info = imfinfo(tiffile);
 nimages = length(info);
 cols = info(1).Width;
 rows = info(1).Height;
 I_tmp = uint16(zeros(rows,cols,nimages));
     
     
   %% Saving  chanel data separately

  for i = 1:1:nimages
     tmp = imread(tiffile,i);    
     I_tmp(:,:,i) = tmp;
  end

 Green_chanel = I_tmp(:,:,1:2:end);
 Red_chanel   = I_tmp(:,:,2:2:end);
 

 save([full_Path_to_Save, '\Green_chanel.mat'],'Green_chanel');
 save([full_Path_to_Save, '\Red_chanel.mat'],'Red_chanel');
 
end
