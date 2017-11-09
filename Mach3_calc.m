function Mach3 = Mach3_calc (gam,M3);

Mach3 = (1+0.5*(gam-1)*M3^2)^(gam/(gam-1));

end
