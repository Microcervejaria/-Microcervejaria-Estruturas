% Grupo MICROCERVEJARIA AUTOMATIZADA
% Projeto Integrador 2
% Universidade de Brasília - Faculdade do Gama

% CÁLCULO DA PERDA DE PRESSÃO ATRAVÉS DA TUBULAÇÃO
%% INPUTS

clear all;close all;clc
format long

vazao_L_m = [26
    26
    26
    2
    26
    2]; % vazao em L/min
vazao = vazao_L_m/(60*1000); % vazao em m3/s

D = [(33.40-2*2.77)/1000
   (33.40-2*2.77)/1000
   (33.40-2*2.77)/1000
   (33.40-2*2.77)/1000
   (33.40-2*2.77)/1000
   (10.20-2*1.24)/1000]; % uniform hydraulic (INTERNAL) diameter of the pipe (m)

v = vazao./(pi*((D/2).^2)); % velocidade

% comprimento total do cano (m) em dado ciclo
L = [2.159 % de 1 para 2
	 2.029 % recirculacao
	 1.665 % de 2 para 3
	 2.094 % resfriamento1
	 2.180 % whirlpool
     0];% resfriamento2
 
dh = [0.853
    0.440
    0.564
    0.160
    0.160
    0]; % diferenca de altura de bombeamento

rho = [980
    1.1*980
    1.1*980
    1.1*997
    1.1*997
    1.1*997]; % fluid density (kg/m3);
% rho1 = 997 agua
% rho2+ = 997*1.1 IPA mais densa

e = 0.00002304; % roughness of the material (Stainless steel)

mu = [0.000431
	0.001
	0.001
	0.002
	0.002
    0.002]; % fluid dynamic (or absolute) viscosity
% mu1 = 0.000431 agua 65 Celsius Table A.8 fox and mcdonalds
% mu2 = 0.001 mosto malteado 65 Celsius
% mu3 = 0.001 mosto malteado 65 Celsius
% mu4 = 0.002 mosto lupulado 25 Celsius
% mu5 = 0.002 mosto lupulado 25 Celsius
% mu6 = 0.002 mosto lupulado 25 Celsius

n_tee_threaded = [3
    4
    4
    3
    4
    0]; % TEEs em dado ciclo

n_long_elbow_threaded = [4
    3
    2
    6
    3
    0]; % cotovelos 90deg long radius em dado ciclo

n_valve_diaf = [3
    3
    3
    3
    3
    0]; % valvulas em dado ciclo

n_contraction = [0
    0
    0
    0
    0
    1]; % contracoes em dado ciclo

n_expansion = [0
    0
    0
    0
    0
    1]; % expansoes em dado ciclo

n_largexit = [0
    0
    1
    1
    1
    0]; % quando o fluido entra num reservatorio ja cheio de liquido

%% Calculos
g = 9.81; % gravity (m/s^2)

Re = rho.*v.*D./mu % Reynolds number
disp('Re_laminar < 2300; Re_turbulent > 2900')

flaminar = 64./Re

f_re3000 = (1./(-1.8*log10( (((e./D)/3.7).^1.11) + (6.9./Re) ))).^2

f = f_re3000;
f(4) = flaminar(4);
% f(6) = flaminar(6);
delta_p_friction = rho.*L.*f.*v.^2./(D*2) % pressure loss in horizontal pipe (Pa)
% e_D = e./D

% r_D = r/D
% disp('Se r/D > 0.15, entao delta_p_entrance eh negligivel')

K_largexit = 1;
K_tee_thread = 1;
K_long_elbow_threaded = 0.2;
K_valve_diaf = 2.3;
K_contraction = 0.5;
K_expansion = 0.85;

K = 0.5 + n_largexit*K_largexit + n_contraction*K_contraction + n_expansion*K_expansion + n_tee_threaded*K_tee_thread + n_long_elbow_threaded*K_long_elbow_threaded + n_valve_diaf*K_valve_diaf;
delta_p_minor = (rho.*K.*v.^2)/2

delta_p_height = rho.*dh.*g

delta_p_total = delta_p_friction + delta_p_minor + delta_p_height % PASCAL

% Potencia_entregue = vazao.*delta_p_total; % WATT

delta_p_total(4) = delta_p_total(4)+delta_p_total(6)+(8.155*1006.28)
% Potencia_entregue(4) = Potencia_entregue(4)+Potencia_entregue(6)

D_bomba = (21.34-2*2.11)/1000;
v_bomba = vazao./(pi*((D_bomba/2).^2));
K_bomba = 0.4 + 0.6;
delta_p_bomba = (rho.*K_bomba.*v_bomba.^2)/2;
delta_p_total = delta_p_total + delta_p_bomba