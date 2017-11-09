function Mach2 = Mach2_calc (gam, M2);

Mach2 = (1 + 0.5 * (gam-1)*M2^2)^(gam/(gam-1));


end
