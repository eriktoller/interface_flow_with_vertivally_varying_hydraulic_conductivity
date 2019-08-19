function [fi] = fi_solver(x,func_Phi_pos,Phi_head, fi_array)
fi = zeros(length(x),1);
for ii = 1:length(x)
    Phi_pos = func_Phi_pos(x(ii));
    diff = abs(Phi_head-Phi_pos);
    if isnan(Phi_pos)
        fi(ii) = NaN;
    else
        fi(ii) = fi_array(diff==min(diff));
    end
end
end

