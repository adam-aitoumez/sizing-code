function ycl = Y_CLcalc(alpha1, alpha2, Lc, sc, dL, ishr)
    
    alpha1 = abs(alpha1);
    alpha2 = abs(alpha2);

    if ishr == 1
        Fb = 0.36;
    else
        Fb = 0.47;
    end

    alpha_m = abs( atand( 2 / (cotd(alpha2) + cotd(alpha1)) )); 
    CL = abs( 2*sc*(cotd(alpha1) - cotd(alpha2)) * sind(alpha_m) ); %coefficient of lift
    Z = abs((CL/sc)^2 * ((sind(alpha2)^2)/(sind(alpha_m)^3))); %Ainley loading parameter

    ycl = Fb * Z * (1/Lc) * (dL*Lc)^0.78;
end