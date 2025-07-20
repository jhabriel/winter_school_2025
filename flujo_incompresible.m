% Flujo incompresible en una columna vertical

% Datos
h0 = 0;  % [m] carga hidraulica en la base de la columna
hL = 5;  % [m] carga hidraulica en la parte superior de la columna
L = 1;   % [m] longitud de la columna

% Mallado
N = 10;  % numero de nodos (incluyendo externos)
dz = L / (N - 1);  % [m] paso espacial
nodos = 0:dz:L; % [m] vector de nodos

% Solucion exacta definida como funcion anonima
h_exact = @(z) ((hL - h0) / L) * z + h0;

% Ensamblaje de la matriz y del vector constante
unos = ones(N-2, 1);

A = spdiags([-1*unos, 2*unos, -1*unos], [-1, 0, 1], N-2, N-2);
% full(A) % muestra la matriz completa
% spy(A)  % muestra la esparsidad de la matriz

b = zeros(N-2, 1);
b(1) = h0;
b(end) = hL;

% Resolviendo el sistema lineal
h_int = A\b;
h = [h0; h_int; hL];

% Graficando
figure();
plot(h_exact(nodos), nodos, '-k', 'LineWidth', 2);
hold on;
plot(h, nodos, '.r', 'MarkerSize', 15);
xlabel('Carga hidraulica [m]', 'FontSize', 16);
ylabel('Elevacion [m]', 'FontSize', 16);
legend({'Solucion exacta', 'Diferencias Finitas'}, ...
        'Location', 'North', ...
        'FontSize', 14);





