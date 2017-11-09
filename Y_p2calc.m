function yp2 = Y_p2calc(sc, alpha2)

    alpha2 = abs(alpha2);

    scMin = 0.224 + (1.575*(alpha2/90)) - ((alpha2 /90)^2);


    X = sc - scMin;
    A = 0.242 - (alpha2/151) + ((alpha2/127)^2);
    
    if alpha2 <= 30
            B = 0.3 + ((30 - alpha2)/50);
    elseif alpha2 > 30
            B = 0.3 + ((30 - alpha2)/275);   
    end
    
    
    C = 0.88 - (alpha2/42.4) + ((alpha2/72.8)^2);
    yp2 = A + (B.*(X.^2)) - (C.*(X.^3));


end

