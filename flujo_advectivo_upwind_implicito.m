% Adveccion 1D en una columna vertical (esquema implicito upwind)
% Flujo descendente, especie quimica ingresando en la parte superior
clear; clc;

%% Parámetros físicos
L   = 1.0;          % [m]  longitud columna
phi = 0.35;         % [-]  porosidad
q   = -1.0e-4;      % [m/s] descarga especifica (negativo = descenso)
v   = q/phi;        % [m/s] velocidad lineal
cin = 1.0;          % concentracion de entrada en z = L
T   = 8e3;          % [s]  tiempo final de simulación

%% Parámetros numéricos
Nz = 101;                   % nodos espaciales
dz = L/(Nz-1);              % [m]  paso espacial
Nt = 400;                   % nodos temporales
dt = T / Nt;                % [s]  paso temporal
C  = abs(v)*dt/dz;          % numero de Courant (solo informativo)
alpha = -v*dt/dz;           % >0 con v<0  (coef. del esquema implicito)

fprintf(['Implicito upwind: dt=%.2e s, Nt=%d, Courant=%.3f, ' ...
         'alpha=%.3f (metodo incondicionalmente estable)\n'], ...
         dt, Nt, C, alpha);

% Adveccion 1D en una columna vertical (esquema IMPLÍCITO upwind)
% Flujo descendente, especie química ingresando en la parte superior
clear; clc;

%% Parametros fisicos
L   = 1.0;          % [m]  longitud columna
phi = 0.35;         % [-]  porosidad
q   = -1.0e-4;      % [m/s] descarga específica (negativo = descenso)
v   = q/phi;        % [m/s] velocidad lineal
cin = 1.0;          % concentracion de entrada en z = L
T   = 8e3;          % [s]  tiempo final de simulación

%% Parametros numericos
Nz = 101;                   % nodos espaciales
dz = L/(Nz-1);              % [m]  paso espacial
Nt = 400;                   % nodos temporales
dt = T / Nt;                % [s]  paso temporal
C  = abs(v)*dt/dz;          % numero de Courant (solo informativo)
alpha = -v*dt/dz;           % >0 con v<0  (coef. del esquema implicito)

fprintf(['Implicito upwind: dt=%.2e s, Nt=%d, Courant=%.3f, ' ...
         'alpha=%.3f (metodo incondicionalmente estable)\n'], ...
         dt,Nt,C,alpha);

%% Mallado e inicialización
z = linspace(0,L,Nz)';      % z=0 abajo, z=L arriba
c = zeros(Nz,1);            % IC: columna limpia

%% Matriz bidiagonal constante para Backward‑Euler upwind
% (1+alpha) c_i^{n+1} - alpha c_{i+1}^{\,n+1} = c_i^{n}
% Notar que la matriz es constante, por lo que no se necesita incluirla
% dentro del bucle for
unos = ones(Nz-1,1);
A = spdiags([(1+alpha)*unos, -alpha*unos], [0 1], Nz-1, Nz-1);

%% Prepara figura
figure(1); clf;
h_num = plot(c,z,'b-','LineWidth',1.4);
hold on;
h_exa = plot(c,z,'k-','LineWidth',1.5);
grid on; xlim([0, 1]); ylim([0, 1]);
xlabel('Concentracion', 'FontSize', 14);
ylabel('z [m]',           'FontSize', 14);
title('Adveccion 1D – Implicito upwind', 'FontSize', 14);
legend({'Numerico','Exacto'}, 'Location','South');

%% Bucle temporal
for n = 1:Nt
    t = n*dt;

    % --- RHS del sistema lineal
    b          = c(1:Nz-1);          % valores antiguos
    b(end)     = b(end) + alpha*cin; % anhade BC en z = L

    % --- Resuelve sistema bidiagonal inferior
    c_new      = A \ b;
    c(1:Nz-1)  = c_new;
    c(end)     = cin;                % BC de entrada

    % --- Solucion exacta instantnea
    zf = L - abs(v)*t;               % posicion del frente
    c_exact = zeros(Nz, 1);
    c_exact(z > zf) = cin;

    % --- Actualiza grafica
    set(h_num,'XData',c,'YData',z);
    set(h_exa,'XData',c_exact,'YData',z);
    drawnow; pause(0.1);
end

% El esquema implicito no explota aunque C>1, pero presenta mas difusion
% numerica si dt es grande. Mejores resultados se obtienen empleando esquemas
% de orden superior.




