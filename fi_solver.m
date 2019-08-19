function [fi] = fi_solver(z,func_Phi_pos,Phi_head,fi_array)
fi = zeros(length(z),1);
    for ii = 1:length(z)
        Phi_pos = func_Phi_pos(z(ii));
        diff = abs(Phi_head-Phi_pos);
        if isnan(Phi_pos)
            fi(ii) = NaN;
        else
            fi(ii) = fi_array(diff==min(diff));
        end
    end
end

