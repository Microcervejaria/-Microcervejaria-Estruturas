%% calcular o diâmetro do tubo do reservatorio %%
%% Dados
Q = 1.667e-4; %m^3/s, vazão volumetrica de água na entrada do trocador de calor
um = 2.46; %m/s velocidade média do escoamento
% Q = A*um
A = Q/um; % área necessária para o escoamento nesta vazão e velocidade

%A = pi *(D/2_^2

D = sqrt(A/pi)*2 % diametro em metros, consequente da área dimensionada

D_mm = D*1000 %diametro em mm, ou seja, um tubo de 1/4 de polegada está dentro dos limites
