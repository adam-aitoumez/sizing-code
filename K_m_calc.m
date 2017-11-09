function Km = K_m_calc (M1)

     if  M1 <= 1.0
            Km = 1;
    elseif M1 > 1
            Km = 1/((-.06553003*M1^2)+(0.092484848*M1)+(0.96965756));  
     end
end
