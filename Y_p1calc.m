function Yp1 = Y_p1calc(sc, alpha2)

    alpha2 = abs(alpha2);

    if alpha2 <= 30
        sc_min = 0.46 + abs(alpha2)/77;
    else
        sc_min = 0.614 + abs(alpha2)/130;
    end
    
    X = sc - sc_min;
    
    preA = 27 - abs(alpha2);
    
    if alpha2 <= 27
        A = 0.025 + preA/530;
    else
        A = 0.025 + preA/3085;
    end
    
    B = 0.1583 - abs(alpha2)/1640;
    C = 0.08*((abs(alpha2)/30)^2 - 1);
    n = 1 + abs(alpha2)/30;
    
    if alpha2 <= 30
        Yp1 = A + B*X^2 + C*X.^3;
    else
        Yp1 = A + B*abs(X).^n;
    end
    
end 
