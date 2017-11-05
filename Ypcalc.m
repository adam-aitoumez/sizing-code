function yp = Ypcalc(alpha1,alpha2,sc,Kmod,Kinc,KM,M1,M2,C_exit,c,nu)

    zeta = (90 - abs(alpha1))/(90 - abs(alpha2));

    yp1 = Y_p1calc(sc, alpha2); yp2 = Y_p2calc(sc, alpha2);

    Kp = KPcalc(M1,M2);
%    KRe = KRecalc(C_exit, c, nu);
    yp = Kmod*Kinc*KM*Kp*KRe*(yp1 + (yp2-yp1)*zeta^2);

end