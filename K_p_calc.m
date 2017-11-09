function Kp = K_p_calc (M1,M2)
    M1_Bar = (M1 + 0.566 - abs(0.566 - M1))/2;
    M2_Bar = (M2 + 1.0 - abs(M2 - 1.0))/2;

    X = (2 * M1_Bar)/ (M1_Bar + M2_Bar + abs(M2_Bar - M1_Bar));

    K1 = 1 - (0.625*(M2_Bar - 0.2 + abs(M2_Bar - 0.2)));

    Kp = 1 - ((1 - K1)*(X^2));
    
end
