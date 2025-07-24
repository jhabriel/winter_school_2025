% adv_diff_column_corrected.m
% Adveccion-difusion 1D en columna vertical (esquema explicito upwind + central)
% Flujo descendente, especie quimica ingresando en la parte superior

clear; clc; close all;

%% Parametros fisicos
L      = 1.0;          % [m] longitud columna
phi    = 0.35;         % [-] porosidad
q      = -1.0e-4;      % [m/s] descarga (negativo = descenso)
v      = q/phi;        % [m/s] velocidad lineal (negativo)
v_abs  = abs(v);       % [m/s] velocidad positiva medida hacia abajo
cin    = 1.0;          % concentracion de entrada en z=L
D_star = 1e-6;         % [m^2/s] difusividad molecular
D_eff  = D_star/phi;   % [m^2/s] difusividad efectiva
T      = 2e03;         % [s] tiempo total de simulacion

%% Parametros numericos
Nz    = 101;
dz    = L/(Nz-1);
z     = linspace(0, L, Nz)';   % z=0 abajo, z=L arriba
Nt    = 200;
dt    = T/Nt;
C     = v_abs*dt/dz;
alpha = D_eff * dt / dz^2;
fprintf('dt=%.2e s, Nt=%d, Courant=%.3f, alpha=%.3f\n', dt, Nt, C, alpha);

%% Inicializacion
c_num = zeros(Nz,1);

%% Grafica inicial
figure(1); clf;
h_num = plot(c_num, z, 'r-', 'LineWidth',1.4); hold on;
h_exa = plot(c_num, z, 'k--','LineWidth',1.5);
grid on; xlim([0, cin]); ylim([0, L]);
xlabel('Concentracion','FontSize',14);
ylabel('z [m]','FontSize',14);
title('Adveccion–Difusion 1D – Analitico vs Numerico','FontSize',14);
legend({'Numerico','Analitico'}, 'Location','South');

%% Bucle temporal
for n = 1:Nt
    t = n*dt;

    % BC en la entrada (z=L)
    c_num(end) = cin;
    c_old = c_num;

    % Upwind (v<0 usa c(i+1)-c(i))
    adv = zeros(Nz,1);
    adv(1:Nz-1) = (-v * dt / dz) .* (c_old(2:Nz) - c_old(1:Nz-1));

    % Central differences para difusion
    diff = zeros(Nz,1);
    diff(2:Nz-1) = alpha * (c_old(3:Nz) - 2*c_old(2:Nz-1) + c_old(1:Nz-2));

    % Actualizacion explicita
    c_num(1:Nz-1) = c_old(1:Nz-1) + adv(1:Nz-1) + diff(1:Nz-1);

    % --- Solucion analitica corregida ---
    % definimos x como distancia desde la entrada (z=L) hacia abajo
    x = L - z;
    arg1 = (x - v_abs*t) ./ (2*sqrt(D_eff*t));
    arg2 = (x + v_abs*t) ./ (2*sqrt(D_eff*t));
    c_exa = 0.5 * cin * ( erfc(arg1) + exp(v .* x ./ D_eff) .* erfc(arg2) );
    disp(c_exa);


    % opcional: asegurar que en region ya invadida c_exa=cin
    %c_exa(x <= v_abs*t) = cin;

    % Actualiza grafica
    set(h_num, 'XData', c_num, 'YData', z);
    set(h_exa, 'XData', c_exa, 'YData', z);
    drawnow; pause(0.1)
end

