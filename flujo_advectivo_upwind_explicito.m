% Adveccion 1D en una columna vertical (esquema explicito upwind)
% Flujo descendente, especie quimica ingresando en la parte superior
clear; clc;

%% Parametros fisicos
L   = 1.0;          % [m] longitud columna
phi = 0.35;         % [-] porosidad
q   = -1.0e-4;      % [m/s] descarga especifica (negativo = descenso)
v   = q/phi;        % [m/s] velocidad lineal
cin = 1.0;          % concentracion de entrada en z = L
T   = 8e03;         % [s] tiempo final de simulacion

%% Parametros numericos
Nz = 101;                   % nodos espaciales
dz = L/(Nz-1);              % [m] paso espacial
Nt = 400;                   % nodos temporales
dt = T / Nt;                % [s] paso temporal
C = abs(v)*dt/dz;           % numero de CFL
fprintf('Explicito upwind: dt=%.2e s, Nt=%d, Courant=%.3f \n',dt,Nt,C);

%% Mallado e inicializacion
z = linspace(0,L,Nz)';      % z=0 abajo, z=L arriba
c = zeros(Nz,1);            % IC: columna limpia

%% Prepara figura
figure(1); clf;
h_num = plot(c,z,'r-','LineWidth',1.4);
hold on;
h_exa = plot(c,z,'k-','LineWidth',1.5);
grid on; xlim([0, 1]); ylim([0, 1]);
xlabel('Concentracion', 'FontSize', 14);
ylabel('z [m]', 'FontSize', 14);
title('Adveccion 1D – Explicito upwind', 'FontSize', 14);
legend({'Numerico', 'Exacto'}, 'Location','South');

%% Bucle temporal
for n = 1:Nt
    t = n*dt;

    % BC de entrada (z=L => ultimo nodo)
    c(end) = cin;

    % esquema explicito upwind (v<0 -> usa i+1)
    c(1:Nz-1) = c(1:Nz-1) + (-v*dt/dz)*(c(2:Nz) - c(1:Nz-1));

    % --- Solucion exacta instantanea
    zf = L - abs(v)*t;      % posicion del frente
    c_exact = zeros(Nz,1);
    c_exact(z > zf) = cin;  % arriba del frente = cin

    % --- Actualiza grafica
    set(h_num,'XData',c,'YData',z);
    set(h_exa,'XData',c_exact,'YData',z);
    drawnow; pause(0.1);
end


% La discrepancia se debe a la difusion numerica introducida por el upwind
% de primer orden. Mejores resultados se obtienen con otros metodos mas
% precisos tales como Lax–Wendroff, MacCormack, QUICK, MUSCL‑TVD
