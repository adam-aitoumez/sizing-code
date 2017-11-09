%% AME 415 Turbine Sizing Code
% by Loren and Kaylyn
%% Cleaning
clear; clc;
%% Constants
g0 = 9.807; % [m/s^2]

m2ft = 3.2801;
hp2w = 745.699872;

%% Inputs
T01 = 500;  % [K]
P01 = 13842074;   % [Pa]
P3 = 11011992;    % [Pa]
mdot = 13.61;   % [kg/s]
M1 = 0.05;

c_p = 15636.11;   % [J/kg/K]
gam = 1.3994;
R = 4123.311;
z = 1;
Zwiff = 0.8;

power = 6215;   % [HP]
N = 50000;   % [RPM]

%% Inputs for Loop 1

eta_TT = 0.90;

Nb_stator = 90; Nb_rotor = 100;

BL = 0.9; 

phi = 0.4;  psi = 1.6;
Rc = 0.4;
KlossN = 0.98; KlossR = 0.98;

%% Calculations


eta_TT_old = 0;

i=0;
ii=0;

 while (abs(eta_TT_old-eta_TT)>0.001)
    
    i=i+1;
    
    
    UoverC0 = sqrt(1 / (phi^2 - 4*(Rc - 1)/eta_TT));

    T3is = T01*(P3/P01)^((gam-1)/gam);
    C0 = sqrt(2*c_p*(T01 - T3is));
    U = C0*UoverC0;

    Ca = U*phi;

    dh0is = (C0^2 - Ca^2)/2;
    W = eta_TT*dh0is;
    power_est = mdot*W/hp2w;

    alpha1 = 90;
    alpha2 = atand(U*Ca/W);

    C2 = Ca/sind(alpha2);
    C2u = Ca/tand(alpha2);
    C2a = Ca;
    W2u = C2u - U;
    W2a = Ca;
    W2 = sqrt(W2u^2 + W2a^2);
    alpha2p = atand(W2a/W2u);

    W3a = Ca;
    W3u = -U;
    W3 = sqrt(W3a^2 + W3u^2);
    C3a = Ca;
    C3u = 0;
    C3 = Ca;
    alpha3 = 90;
    alpha3p = atand(W3a/W3u);

    T2 = T01 - (C2^2 / 2 / c_p);
    a2 = sqrt(gam*R*T2);
    M2 = C2/a2;

    Mw2 = W2/a2;
    Tw2 = T2*(1 + 0.5*(gam-1)*Mw2^2);
    T3 = Tw2 - (W3^2 / 2 / c_p);
    a3 = sqrt(gam*R*T3);
    Mw3 = W3/a3;
    M3 = C3/a3;


    P02 = P01*KlossN;
    P2 = P02/(1 + 0.5*(gam-1)*M2^2)^(gam/(gam-1));
    Pw2 = P2*(1 + 0.5*(gam-1)*Mw2^2)^(gam/(gam-1));
    Pw3 = Pw2*KlossR;
    P3veri = Pw3/(1+ 0.5*(gam-1)*Mw3^2)^(gam/(gam-1));

    rho2 = P2/z/R/T2;
    rho3 = P3/z/R/T3;
    rhor = (rho2+rho3)/2;

    A2 = mdot/rho2/C2a;
    A3 = mdot/rho3/C3a;

    rmean = 30*U/pi/N;

    L2 = A2/2/pi/rmean/BL;
    L3 = A3/2/pi/rmean/BL;

    %% Stator

    sc0_stator = 0.427 + alpha2/58 - (alpha2/93)^2;
    sBz_stator = Zwiff/(2*(sind(alpha2)^2)*(cotd(alpha2) - cotd(alpha1)));

    s_stator = 2*pi*rmean/Nb_stator;

    Bz_stator = s_stator/sBz_stator; c_stator = s_stator/sc0_stator;
    betaS_stator = asind(Bz_stator/c_stator);

    %% Rotor
    sc0_rotor = 0.427 + abs(alpha3p)/58 - (alpha3p/93)^2;
    sc1_rotor = 0.224 + (1.575 - abs(alpha3p)/90)*(abs(alpha3p)/90);

    zeta = (90 - abs(alpha2p))/(90 - abs((alpha3p)));

    sc_opt = sc0_rotor + (sc1_rotor - sc0_rotor)*abs(zeta)*zeta;

    sBz_rotor = Zwiff/(2*(sind(alpha3p)^2)*(cotd(alpha2p) - cotd(alpha3p)));

    s_rotor = 2*pi*rmean/Nb_rotor;

    Bz_rotor = s_rotor/sBz_rotor; c_rotor = s_rotor/sc_opt;
    betaS_rotor = asind(Bz_rotor/c_rotor);
    
    


    eta_TT_old = eta_TT;
    
    
    %% 3D Design

    r_tip_stator = rmean + 0.5*L2;  r_tip_rotor = rmean + 0.5*L3;
    r_hub_stator = rmean - 0.5*L2;  r_hub_rotor = rmean - 0.5*L3;

    U_tip = r_tip_stator*pi*N/30;
    U_hub = r_hub_stator*pi*N/30;

    C2u_tip = rmean*C2u/r_tip_stator;   C2a_tip = Ca;
    C2u_hub = rmean*C2u/r_hub_stator;   C2a_hub = Ca;
    C2_tip = sqrt(C2u_tip^2 + C2a_tip^2);
    C2_hub = sqrt(C2u_hub^2 + C2a_hub^2);

    W2u_tip = C2u_tip - U_tip;  W2a_tip = C2a_tip;
    W2_tip = sqrt(W2u_tip^2 + W2a_tip^2);
    W2u_hub = C2u_hub - U_hub;  W2a_hub = C2a_hub;
    W2_hub = sqrt(W2u_hub^2 + W2a_hub^2);

    alpha2_tip = asind(C2a_tip/C2_tip);  alpha2p_tip = atand(W2a_tip/W2u_tip);
    alpha2_hub = asind(C2a_hub/C2_hub);  alpha2p_hub = atand(W2a_hub/W2u_hub);

    C3u_tip = 0;    C3a_tip = Ca;
    C3u_hub = 0;    C3a_hub = Ca;
    C3_tip = sqrt(C3u_tip^2 + C3a_tip^2);
    C3_hub = sqrt(C3u_hub^2 + C3a_hub^2);

    W3u_tip = C3u_tip - U_tip;  W3a_tip = C3a_tip;
    W3_tip = sqrt(W3u_tip^2 + W3a_tip^2);
    W3u_hub = C3u_hub - U_hub;  W3a_hub = C3a_hub;
    W3_hub = sqrt(W3u_hub^2 + W3a_hub^2);

    alpha3_tip = asind(C3a_tip/C3_tip);  alpha3p_tip = atand(W3a_tip/W3u_tip);
    alpha3_hub = asind(C3a_hub/C3_hub);  alpha3p_hub = atand(W3a_hub/W3u_hub);

    tipVals = [C2a_tip, C2u_tip, C2_tip, W2a_tip, W2u_tip, W2_tip, ...
                    alpha2_tip, alpha2p_tip, C3a_tip, C3u_tip, C3_tip, ...
                    W3a_tip, W3u_tip, W3_tip, alpha3_tip, alpha3p_tip];
    hubVals = [C2a_hub, C2u_hub, C2_hub, W2a_hub, W2u_hub, W2_hub, ...
                    alpha2_hub, alpha2p_hub, C3a_hub, C3u_hub, C3_hub, ...
                    W3a_hub, W3u_hub, W3_hub, alpha3_hub, alpha3p_hub];             

    meanVals = [C2a, C2u, C2, W2a, W2u, W2, alpha2, alpha2p, ...
                    C3a, C3u, C3, W3a, W3u, W3, alpha3, alpha3p];


    %% Y Loss Stator
        %Yp1
        Yp1_stator = Y_p1calc(sc0_stator,alpha2);


        %Yp2
        Yp2_stator = Y_p2calc(sc0_stator,alpha2);

        %Yp
        mu_stator = 8.756E-6*((293.85+72)/(T2+72))*(T2/293.85)^(3/2);
        nu_stator = mu_stator/rho2;
        Yp_stator = Ypcalc(alpha1,alpha2,M1,M2,C2,c_stator,nu_stator,Yp1_stator,Yp2_stator);

        %Ys Stator
        Lc_stator = L2/c_stator;
        Ys_stator = Y_Scalc(alpha1,alpha2,alpha1,Lc_stator,sc0_stator);

        % Ycl Stator

        Ycl_stator = 0;

        % Yex Stator

        Yex_stator = 0;


        %Y Total
        Y_N = Yp_stator + Ys_stator + Ycl_stator + Yex_stator;


    %% Y Loss Rotor
        %Yp1
        Yp1_rotor = Y_p1calc(sc_opt,alpha3p);


        %Yp2
         Yp2_rotor = Y_p2calc(sc_opt,alpha3p);


        %Yp
        mu_rotor = 8.756E-6*((293.85+72)/(T3+72))*(T3/293.85)^(3/2);
        nu_rotor = mu_rotor/rhor;
        Yp_rotor = Ypcalc(alpha2p,alpha3p,Mw2,Mw3,W3,c_rotor,nu_rotor,Yp1_rotor,Yp2_rotor);


        % Ys Rotor

        Lc_rotor = L3/c_rotor;
        Ys_rotor = Y_Scalc(alpha2p,alpha3p,alpha2p,Lc_rotor,sc_opt);

        % Ycl Rotor

        dL_rotor = 0.00025 / L3;
        Ycl_rotor = Y_CLcalc(alpha2p,alpha3p,Lc_rotor,sc_opt,dL_rotor,1);


       % Yex Rotor

       Yex_rotor = 0;

       % Y Total

        Y_R = Yp_rotor + Ys_rotor + Ycl_rotor + Yex_rotor;

        %% First Pass

        M2_func = Mach2_calc(gam,M2);

        P2_it = P01/((M2_func*(Y_N+1))-Y_N);

        P02_it = P2_it * M2_func;

        Mw3_func = Machw3_calc (gam,Mw3);

        M3_func = Mach3_calc (gam,M3);

        P3_it = Pw2/((Mw3_func*(Y_R+1))-Y_R);

        Pw3_it = P3_it * Mw3_func;

        P03_it = P3_it * M3_func;

        Kloss_stator = P02_it/P01;

        Kloss_rotor = Pw3_it/Pw2;
        
        KlossN = Kloss_stator;
        KlossR = Kloss_rotor;

        while (abs (P3veri - P3) >100)
    
            ii=ii+1;
           
            
            C2_it = C2*(P3_it/P3);

            Ca_it = C2_it*sind(alpha2);

            C2 = C2_it;
            Ca = Ca_it;

            C2u = Ca/tand(alpha2);
            C2a = Ca;
            W2u = C2u - U;
            W2a = Ca;
            W2 = sqrt(W2u^2 + W2a^2);
            alpha2p = atand(W2a/W2u);

            W3a = Ca;
            W3u = -U;
            W3 = sqrt(W3a^2 + W3u^2);
            C3a = Ca;
            C3u = 0;
            C3 = Ca;
            alpha3 = 90;
            alpha3p = atand(W3a/W3u);

            T2 = T01 - (C2^2 / 2 / c_p);
            a2 = sqrt(gam*R*T2);
            M2 = C2/a2;

            Mw2 = W2/a2;
            Tw2 = T2*(1 + 0.5*(gam-1)*Mw2^2);
            T3 = Tw2 - (W3^2 / 2 / c_p);
            a3 = sqrt(gam*R*T3);
            Mw3 = W3/a3;
            M3 = C3/a3;


            P02 = P01*Kloss_stator;
            P2 = P02/(1 + 0.5*(gam-1)*M2^2)^(gam/(gam-1));
            Pw2 = P2*(1 + 0.5*(gam-1)*Mw2^2)^(gam/(gam-1));
            Pw3 = Pw2*Kloss_rotor;
            P3veri = Pw3/(1+ 0.5*(gam-1)*Mw3^2)^(gam/(gam-1));

            rho2 = P2/z/R/T2;
            rho3 = P3/z/R/T3;
            rhor = (rho2+rho3)/2;

            A2 = mdot/rho2/C2a;
            A3 = mdot/rho3/C3a;

            rmean = 30*U/pi/N;

            L2 = A2/2/pi/rmean/BL;
            L3 = A3/2/pi/rmean/BL;

            %% Stator

            sc0_stator = 0.427 + alpha2/58 - (alpha2/93)^2;
            sBz_stator = Zwiff/(2*(sind(alpha2)^2)*(cotd(alpha2) - cotd(alpha1)));

            s_stator = 2*pi*rmean/Nb_stator;

            Bz_stator = s_stator/sBz_stator; c_stator = s_stator/sc0_stator;
            betaS_stator = asind(Bz_stator/c_stator);

            %% Rotor
            sc0_rotor = 0.427 + abs(alpha3p)/58 - (alpha3p/93)^2;
            sc1_rotor = 0.224 + (1.575 - abs(alpha3p)/90)*(abs(alpha3p)/90);

            zeta = (90 - abs(alpha2p))/(90 - abs((alpha3p)));

            sc_opt = sc0_rotor + (sc1_rotor - sc0_rotor)*abs(zeta)*zeta;

            sBz_rotor = Zwiff/(2*(sind(alpha3p)^2)*(cotd(alpha2p) - cotd(alpha3p)));

            s_rotor = 2*pi*rmean/Nb_rotor;

            Bz_rotor = s_rotor/sBz_rotor; c_rotor = s_rotor/sc_opt;
            betaS_rotor = asind(Bz_rotor/c_rotor);
            
            
            %% 3D Design

            r_tip_stator = rmean + 0.5*L2;  r_tip_rotor = rmean + 0.5*L3;
            r_hub_stator = rmean - 0.5*L2;  r_hub_rotor = rmean - 0.5*L3;

            U_tip = r_tip_stator*pi*N/30;
            U_hub = r_hub_stator*pi*N/30;

            C2u_tip = rmean*C2u/r_tip_stator;   C2a_tip = Ca;
            C2u_hub = rmean*C2u/r_hub_stator;   C2a_hub = Ca;
            C2_tip = sqrt(C2u_tip^2 + C2a_tip^2);
            C2_hub = sqrt(C2u_hub^2 + C2a_hub^2);

            W2u_tip = C2u_tip - U_tip;  W2a_tip = C2a_tip;
            W2_tip = sqrt(W2u_tip^2 + W2a_tip^2);
            W2u_hub = C2u_hub - U_hub;  W2a_hub = C2a_hub;
            W2_hub = sqrt(W2u_hub^2 + W2a_hub^2);

            alpha2_tip = asind(C2a_tip/C2_tip);  alpha2p_tip = atand(W2a_tip/W2u_tip);
            alpha2_hub = asind(C2a_hub/C2_hub);  alpha2p_hub = atand(W2a_hub/W2u_hub);

            C3u_tip = 0;    C3a_tip = Ca;
            C3u_hub = 0;    C3a_hub = Ca;
            C3_tip = sqrt(C3u_tip^2 + C3a_tip^2);
            C3_hub = sqrt(C3u_hub^2 + C3a_hub^2);

            W3u_tip = C3u_tip - U_tip;  W3a_tip = C3a_tip;
            W3_tip = sqrt(W3u_tip^2 + W3a_tip^2);
            W3u_hub = C3u_hub - U_hub;  W3a_hub = C3a_hub;
            W3_hub = sqrt(W3u_hub^2 + W3a_hub^2);

            alpha3_tip = asind(C3a_tip/C3_tip);  alpha3p_tip = atand(W3a_tip/W3u_tip);
            alpha3_hub = asind(C3a_hub/C3_hub);  alpha3p_hub = atand(W3a_hub/W3u_hub);

            tipVals = [C2a_tip, C2u_tip, C2_tip, W2a_tip, W2u_tip, W2_tip, ...
                            alpha2_tip, alpha2p_tip, C3a_tip, C3u_tip, C3_tip, ...
                            W3a_tip, W3u_tip, W3_tip, alpha3_tip, alpha3p_tip];
            hubVals = [C2a_hub, C2u_hub, C2_hub, W2a_hub, W2u_hub, W2_hub, ...
                            alpha2_hub, alpha2p_hub, C3a_hub, C3u_hub, C3_hub, ...
                            W3a_hub, W3u_hub, W3_hub, alpha3_hub, alpha3p_hub];             

            meanVals = [C2a, C2u, C2, W2a, W2u, W2, alpha2, alpha2p, ...
                            C3a, C3u, C3, W3a, W3u, W3, alpha3, alpha3p];




            M2_func = Mach2_calc(gam,M2);

            P2_it = P01/((M2_func*(Y_N+1))-Y_N);

            P02_it = P2_it * M2_func;

            Mw3_func = Machw3_calc (gam,Mw3);

            M3_func = Mach3_calc (gam,M3);

            P3_it = Pw2/((Mw3_func*(Y_R+1))-Y_R);

            Pw3_it = P3_it * Mw3_func;

            P03_it = P3_it * M3_func;

            Kloss_stator = P02_it/P01;

            Kloss_rotor = Pw3_it/Pw2;
            
       end
    
   T03 = T3*(1+0.5*(gam-1)*M3^2);
   eta_TT = (1-(T03/T01))/(1-(P03_it/P01)^((gam-1)/gam));
   
end
   
