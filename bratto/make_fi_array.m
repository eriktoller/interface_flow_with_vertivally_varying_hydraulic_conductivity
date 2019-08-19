function [fi_array,Phi_head] = make_fi_array(func_Phi_head,fi_start,fi_end,inc)
fi_array = linspace(fi_start,fi_end,((fi_end-fi_start)/inc));
Phi_head = func_Phi_head(fi_array);
end

