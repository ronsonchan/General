% output functino that output case comps
function [] = output(str,n,final_state)
    if(strcmp(str,'X'))
       [meas_X,~,~] = meas(n,'XYZ',final_state) 
    elseif(strcmp(str,'Y'))
       [~,meas_Y,~] = meas(n,'XYZ',final_state) 
    elseif(strcmp(str,'Z'))
       [~,~,means_Z] = meas(n,'XYZ',final_state)
    elseif(strcmp(str,'XY'))
       [meas_X,meas_Y,~] = meas(n,'XYZ',final_state)
    elseif(strcmp(str,'YZ'))
       [~,meas_Y,meas_Z] = meas(n,'XYZ',final_state)
    elseif(strcmp(str,'XZ'))
       [meas_X,~,meas_Z] = meas(n,'XZ',final_state)
    elseif(strcmp(str,'XYZ'))
       [meas_X,meas_Y,meas_Z] = meas(n,'XYZ',final_state)
    else
        error('Specified axis does not exist')
    end
end