function yp = Ypcalc(alpha1,alpha2,M1,M2,C2,c,nu,yp1,yp2)

    alpha1= abs(alpha1);
    alpha2= abs(alpha2);

    Kmod = 1;
    Kinc = 1;
    Km = K_m_calc(M1);
    Kp = K_p_calc(M1,M2);
    Kre = K_re_calc (C2,c,nu);
    Kte = 1;
    

    zeta = (90 - abs(alpha1))/(90 - abs(alpha2));
    


    yp = Kmod*Kinc*Km*Kp*Kre*Kte*(yp1 + (yp2-yp1)*zeta^2);

end