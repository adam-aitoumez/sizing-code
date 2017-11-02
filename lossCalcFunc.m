function delta = lossCalcFunc(lossEst, REF)
gam = REF(1);
P01 = REF(2);
M2 = REF(3);
Mw2 = REF(4);
Mw3 = REF(5);
P3 = REF(6);

eta = 1 - lossEst;
P02 = P01*eta;
P2 = P02/(1 + 0.5*(gam-1)*M2^2)^(gam/(gam-1));
Pw2 = P2*(1 + 0.5*(gam-1)*Mw2^2)^(gam/(gam-1));
Pw3 = eta*Pw2;
P3ver = Pw3/(1 + 0.5*(gam-1)*Mw3^2)^(gam/(gam-1));

delta = P3ver - P3;
end