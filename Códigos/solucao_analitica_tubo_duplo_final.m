%% Cálculo para determinação dos diâmetro dos tubos do trocador de calor tubo duplo %%
clc
clear all;
close all;
% format bank;
format long;
%Q_ponto_q = mq_ponto * Cp_quente(Tq_entrada - Tq_saida) - eq 1 para o fluido
%quente
%Q_ponto_f mf_ponto * Cp_f(Tf_entrada - Tf_saida) - eq 2 para o fluido
%frio
% Q = U*A*DTML;

%% Dados gerais
    %DADOS:
     vaz_vol_mosto = 3.333e-5 ; %21 l/min =  3.33e-5 m^3/s
     rho_mosto = 997; % 1086.73 %997; % kg/m^31100;% 
     vaz_mas_mosto = vaz_vol_mosto * rho_mosto;
     Cp_quente = 4.184; %KJ/Kg K
     Cp_frio = 4.184; %KJ/Kg K
     Tq_entrada = 100; % em °C
     Tq_saida = 25; % em °C
     
     Q_ponto_q = vaz_mas_mosto * Cp_quente * (Tq_entrada-Tq_saida); % em KW
     Q_ponto_q_joulepersec = vaz_mas_mosto * Cp_quente * (Tq_entrada-Tq_saida)*1000; % em j/s
     
     vaz_vol_agua = 1.667e-4 ; %10 L/min =  1.667e-4 m^3/s
     rho_agua = 997; %kg/m^3
     vaz_mas_agua = vaz_vol_agua* rho_agua;
     Tf_entrada = 20; % °C
     Tf_saida_K =   (Tf_entrada + 273 ) + (Q_ponto_q/vaz_mas_agua/Cp_frio);
     
    Tf_saida =  Tf_saida_K - 273;

    DT1 = Tq_entrada - Tf_saida; %delta T
    DT2 = Tq_saida - Tf_entrada;
    DTML = (DT1-DT2)/(log(DT1/DT2));
    %% Dados dos tubos escolhidos:
   %Tubo interno
        espessura_tubo_mosto = 1.24e-3; %espessura 
        Dia_int_tubo_mosto = 10.29e-3 - (2*espessura_tubo_mosto); %em metros, diâmetro interno do tubo interno
        Dia_ext_tubo_mosto = 10.29e-3 ;  %em metros, diâmetro externo do tubo interno
   %Tubo externo:
        espessura_tubo_agua = 1.65e-3;

        Dia_int_tubo_agua = 17.15e-3 - (2*espessura_tubo_agua); %em metros, diâmetro interno do tubo externo
        Dia_ext_tubo_agua = 17.15e-3 ; %em metros, diâmetro externo do tubo externo

        k_tubo = 14.9;% coeficiente de condutividade termica do material proposto - AÇO INOX AISI 304. a 300K

    
    %% Análise para o tubo interno
    mi_mosto = 4.5125e-4 ;% em N*s/m^2
    Re_i = (4*vaz_mas_mosto)/(pi*Dia_int_tubo_mosto*mi_mosto); %admensional
    Pr_mosto = 2.89; %a 62.5°C
    f_i = (0.790*log(Re_i)-1.64)^-2; % admensional
    Nud_i = (f_i/8*(Re_i-1000)*Pr_mosto)/(1+12.7*((Pr_mosto^(2/3)-1)*(f_i/8)^0.5)); % Em w/m^2 °C
    d_espira = 0.25; %m
    lambda = d_espira / Dia_int_tubo_mosto;
     
    f_c_i = (0.0075/sqrt(lambda)) + f_i;
    
    Nud_i_c = Nud_i * (1 + (3.4*(1/lambda)));
    k_mosto = 0.659;
    h_i = Nud_i_c * k_mosto / Dia_int_tubo_mosto;
    
    %% para o tubo externo: 
%     D_ext = 15.875e-3;
%     d_int_e = 13.875e-3; 
%     d_ext_e = 9.525e-3;
%     
    mi_agua = 0.849e-3; % N*s/m^2
    Re_e = (4*vaz_mas_agua)/(pi*(Dia_int_tubo_agua+Dia_ext_tubo_mosto)*mi_agua); %admensional
    
    Pr_agua = 5.75; %a 27.5°C
    f_e = (0.790*log(Re_e)-1.64)^-2; % admensional
    Nud_e = (f_e/8*(Re_e-1000)*Pr_agua)/(1+12.7*((Pr_agua^(2/3)-1)*(f_e/8)^0.5)); % Em w/m^2 °C
    % no caso da regiao anular, utiliza-se o diametro hidraulico
    diam_hidr=(Dia_int_tubo_agua^2-Dia_ext_tubo_mosto^2)/Dia_ext_tubo_mosto; % calculo do diametro hidraulico
    lambda_e = d_espira / diam_hidr;
     
    f_c_e = (0.0075/sqrt(lambda_e)) + f_e;
    
    Nud_i_e = Nud_i * (1 + (3.4*(1/lambda_e)));
    k_agua = 0.659;
    h_e = Nud_i_e * k_agua / diam_hidr;
    
    %% calculo do coeficiente global de calor 
    Rdi = 1.76e-4; % em w/m^2 K
    U = ((Dia_ext_tubo_mosto/(h_i*Dia_int_tubo_agua))+((Dia_ext_tubo_agua)*(log(Dia_ext_tubo_agua/Dia_int_tubo_agua))/(2*k_tubo))+(Rdi)+(1/h_e))^-1; % em w/m^2 K
    %% calculando a área final de troca de calor por Q=U*A*DTML
    
    A= (Q_ponto_q_joulepersec)/(DTML*U); % em m^2
    
    % A=pi*r*L
    L = A/(pi*Dia_ext_tubo_mosto) % em metos
    
    %% Calculos de eficiência
    Cmosto = vaz_mas_mosto*Cp_quente; %kW/K
    Cagua = vaz_mas_agua*Cp_frio; %kW/K
    
     if Cmosto < Cagua
        Cmin = Cmosto;%kW/K
        Cmax = Cagua;
     else
        Cmin = Cagua; %kW/K
        Cmax = Cmosto;
     end
         
    DTmax = Tq_entrada - Tf_entrada;
     
    Qmax = Cmin* DTmax; %em kW
    
    eff = Q_ponto_q/Qmax; % eficiência do trocador de calor
    Cmin_em_wattporkelvin = Cmin *1000;
    
    NTU = U*A/Cmin_em_wattporkelvin % numero de unidades de transferência eq 11-32 çengel
    
    c = Cmin/Cmax;
    %calculos de eficiência para um trocador de calor tubo duplo para
    %contrafluxo, tabela 11-4 do çengel
    if c < 1
        
        Eff2 = (1-exp(-NTU*(1-c)))/(1-c*exp(-NTU*(1-c)))
        
    else c == 1  
        
        Eff2 = NTU/(NTU+1)
    
    end
     
     %% calculo de n de espiras:
     % L = N*SQRT(C^2+P^2)
     p = 18e-3;
     circ = 2*pi*12.5e-2
     
     n = L/(sqrt(circ^2 + p^2))
     h = n*p
     
     %% calculos gerais
     A_entrada_agua = (pi/4)*(Dia_int_tubo_agua^2 - Dia_ext_tubo_mosto^2); % em m^2
     
     A_entrada_mosto = (pi/4)* (Dia_int_tubo_mosto^2); % em m^2
     
     
     %% calculo das velocidades dos escoamentos
     %V = u * A >> u = V/A velocidade media = vazao sobre area
     u_mosto = vaz_vol_mosto/A_entrada_mosto
     u_agua = vaz_vol_agua/A_entrada_agua
      
     % perda de carga: f = deltaP*Dt/2*rho*velmédia^2*comprimento
     % deltap = fator de atrito* 2 * rho * velmedia^2 * L / Dt 
     
     deltaPagua_perda = (f_c_e*(1/Dia_int_tubo_agua)*rho_agua*((u_agua^2)/2))
     
     deltaPmosto_perda = (f_c_e*(1/Dia_int_tubo_mosto)*rho_mosto*((u_mosto^2)/2))
     % 
     deltaP_mosto = rho_mosto*9.81*h + deltaPmosto_perda
     deltaP_agua = rho_agua*9.81*h + deltaPagua_perda
     