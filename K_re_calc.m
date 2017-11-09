function Kre = K_re_calc (C2,c,nu)
    
    Re_c = (C2*c)/nu;
    e = 5E-6;
    Re_r = 100 * (c/e);
    
    if Re_c <= 1E5
        Kre = sqrt ((1E5/Re_c));
        
    elseif Re_c > 5E5
        Kre = 1 + ((((log10(5E5)/log10(Re_r))^2.58) - 1) *(1 - (5E5/Re_c)));
        
    else
        Kre = 1;
    end
end
