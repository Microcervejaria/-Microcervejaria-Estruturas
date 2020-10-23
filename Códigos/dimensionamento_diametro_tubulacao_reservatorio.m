%% calcular o di�metro do tubo do reservatorio %%
%% Dados
Q = 1.667e-4; %m^3/s, vaz�o volumetrica de �gua na entrada do trocador de calor
um = 2.46; %m/s velocidade m�dia do escoamento
% Q = A*um
A = Q/um; % �rea necess�ria para o escoamento nesta vaz�o e velocidade

%A = pi *(D/2_^2

D = sqrt(A/pi)*2 % diametro em metros, consequente da �rea dimensionada

D_mm = D*1000 %diametro em mm, ou seja, um tubo de 1/4 de polegada est� dentro dos limites
