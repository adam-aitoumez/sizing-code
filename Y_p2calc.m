function yp2 = Y_p2calc(sc, alpha2)

    scStator = sc;

    scMinYp2Nozzle = scRotor(alpha2);


    X_Yp2_Nozzle = scStator - scMinYp2Nozzle;
    A_Yp2_Nozzle = 0.242 - (alpha2/151) + ((alpha2/127)^2);
    
    if alpha2 <= 30
            B_Yp2_Nozzle = 0.3 + ((30 - alpha2)/50);
    elseif alpha2 > 30
            B_Yp2_Nozzle = 0.3 + ((30 - alpha2)/275);   
    end 
    

    C_Yp2_Nozzle = 0.88 - (alpha2/42.4) + ((alpha2/72.8)^2);
    yp2 = A_Yp2_Nozzle + (B_Yp2_Nozzle.*(X_Yp2_Nozzle.^2)) + (C_Yp2_Nozzle.*(X_Yp2_Nozzle.^3));


end

function scRotorMin = scRotor(alpha2)

    x = abs(alpha2);
    scRotorMin =  0.224 + (1.575*(x/90)) - ((x /90)^2);

end