function Ys = Y_Scalc(alpha1,alpha2,beta1,Lc,sc)

    alpha1 = abs(alpha1);
    alpha2 = abs(alpha2);
    
    alpham = abs ( atand ( 2/ ((1/tand(alpha2)) + (1/tand(alpha1)) ))) ;

   
    CL = abs ( 2 * sc * ((1/tan(alpha1*pi/180)) - (1/tan(alpha2*pi/180)) )* sin (alpham *pi/180));

    
    Z = abs((CL / (sc))^2* ((sin (alpha2*pi/180))^2/(sin(alpham*pi/180))^3));

    
    if Lc >= 2
        Far = 1/Lc;
    else
        Far = 0.5 * (2* (1/Lc))^.7;
    end
    
   
    Ysbar = abs(0.0334 * Far * Z * sind(alpha2) / sind(beta1));
    
%     if (nargin == 8)
%         Ks = 1 - (1 - KP)*(bzL^2 / (1 + bzL^2));
%         Ys = Ks*KRe*(Ysbar^2)/(1 + 7.5*Ysbar^2);
%     else
%         Ys = Ysbar;
%     end
    Ys = Ysbar;
end