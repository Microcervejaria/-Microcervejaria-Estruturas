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
L = [0.4+0.8+0.05+0.465+0.278+(2*pi*0.0254/4)+0.062+0.1362+0.0254 % de 1 para 2
	 0.4+0.8+0.05+0.278+(2*pi*0.0254/4)+0.062+0.1362+0.0254 % recirculacao
	 0.4+0.4+0.278+0.1362+0.0254+0.142+0.0254+0.0381 % de 2 para 3
	 (2*pi*0.0254)+0.11+0.0885+0.0516+0.286+0.2485+0.4+0.0381 % resfriamento1
	 0.4+0.4+(2*pi*0.0254/4)+0.2885+0.1362+0.0254+0.142+0.0254+0.0381 % whirlpool
     0];% resfriamento2

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

n_tee_threaded = [4
    4
    4
    3
    3
    0]; % TEEs em dado ciclo

n_long_elbow_threaded = [2
    1
    0
    4
    1
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

dh = [0.85
    0.45
    0.72
    0.32
    0.32
    0]; % diferenca de altura de bombeamento

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
% K_entrance = ;
% delta_p_entrance = (rho*K_entrance*v^2)/2

K_tee_thread = 2;
K_long_elbow_threaded = 0.7;
K_valve_diaf = 2.3;
K_contraction = 0.5;
K_expansion = 0.85;

K = n_contraction*K_contraction + n_expansion*K_expansion + n_tee_threaded*K_tee_thread + n_long_elbow_threaded*K_long_elbow_threaded + n_valve_diaf*K_valve_diaf;
delta_p_minor = (rho.*K.*v.^2)/2

delta_p_height = rho.*dh.*g

delta_p_total = delta_p_friction + delta_p_minor + delta_p_height; % PASCAL

% Potencia_entregue = vazao.*delta_p_total; % WATT

delta_p_total(4) = delta_p_total(4)+delta_p_total(6)
% Potencia_entregue(4) = Potencia_entregue(4)+Potencia_entregue(6)