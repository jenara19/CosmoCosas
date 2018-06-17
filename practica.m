% APARTADO 2
% La función distancia luminosidad es la integral de la función g
g = @(qz, omm, om) (1./sqrt(omm.*(1 + qz).^3 + (1 - omm).*(1 + qz).^(3 + 3.*om)));

zmax = 1.5;
distluz = @(omm, om, z) (1 + z).*(integral(@(qz) heaviside(z-qz).*g(qz, omm, om), 0, zmax,'ArrayValued',true,'RelTol',1e-2));


% APARTADO 3
z = linspace(0, 1.5, 100);

dist_apartadotres = distluz(0.3, -1, z);

figure
plot(z, dist_apartadotres);

valor_pedido = distluz(0.3, -1, 0.5);
disp(fprintf('El valor pedido de distancia en apartado 3 es: %f', valor_pedido)) ;


% APARTADO 4
load('Union2.mat');

% Función Mm a la que hay que proporcionar el vector de distancias
% calculado
Mm = @(omm, om) (sum((m - 5.*log10(distluz(omm, om, z)))./(dm.^2))/sum(1./dm.^2));
disp(Mm(0.3, -1));

% Función chi cuadrado marginal, sustituyo lo que de Mm
chicuadrado_marg = @(omm, om) (sum(((m - 5.*log10(distluz(omm, om, z)) - Mm(omm, om)).^2)./(dm.^2)));

chi_pedida = chicuadrado_marg(0.3, -1);
disp(fprintf('El valor pedido de chi cuadrado en apartado 4 es: %f', chi_pedida));


% APARTADO 5
omm_pruebas = linspace(0, 1, 100);
y = zeros(1, 100);
for i=1:length(omm_pruebas)
        y(i) = chicuadrado_marg(omm_pruebas(i), -1);
        
end
figure
plot(omm_pruebas, y);

% Resultado: omega que minimiza: 0.2727
% En ese punto, chicuadrado vale 541.8
omm_mejor = 0.2727;
chi_mejor = 541.8;

% APARTADO 6
% En 1 sigma, chicuadrado tiene que valer 542.8

hold on
plot(omm_pruebas, 542.8*ones(1, length(omm_pruebas)));

% Para ese valor, omega M vale 0.2929
error = 0.2929 - 0.2727;
disp(fprintf('El error es: %f', error));


% APARTADO 7
omm_sitter = 1;
om_sitter = -1; % esta da igual lo que valga
chi_sitter = chicuadrado_marg(omm_sitter, om_sitter);

omm_mejor = 0.2727;
% comparar los valores de chi

% Representar los módulos de distancia
moddist_mejor = 5*log10(distluz(omm_mejor, -1, z)) + Mm(omm_mejor, -1);
moddist_sitter = 5*log10(distluz(omm_sitter, om_sitter, z)) + Mm(omm_sitter, om_sitter);
moddist_exp = m;

figure
errorbar(z, moddist_exp, dm, 'g');
hold on
plot(z, moddist_sitter, 'b');
hold on
plot(z, moddist_mejor, 'r');

% Calculamos chi cuadrado para los parámetros de Sitter, restamos, y
% el valor de sigma es la raíz de esa diferencia
delta_chi = abs(chi_sitter - chi_mejor);
no_sigmas = sqrt(delta_chi);

disp(fprintf('El valor de sigma del apartado 7 es: %f', no_sigmas));

% APARTADO 8
Hzeroinv = 4285;
% Calculamos la distancia luminosidad, en Mpc (mult por TAL)
dist_apocho = Hzeroinv*distluz(omm_mejor, -1, 1);
disp(fprintf('Distancia luminosidad del apartado 8: %f', dist_apocho));

% Ahora calculamos la distancia comóvil, que es dividir por 1+z=2
dist_com = dist_apocho/2;
disp(fprintf('Distancia comóvil es: %f', dist_com));


% APARTADO 9
% x = linspace(-10, 1, 100);
% r = zeros(1, length(x));
% for i=1:length(x)
%     r(i) = chicuadrado_marg(omm_mejor, x(i));
% end
% 
% figure
% plot(x, r);

N = 100;
% x son valores de omega grande
% y son valores de omega pequeña
x = linspace(0, 5, N);
y = linspace(-2, 2, N);
[X, Y] = meshgrid(x, y);

Z = zeros(length(y),length(x));

for i=1:length(x)
    for j=1:length(y)
        Z(i, j) = chicuadrado_marg(x(i), y(j)); 
    end
end

[zmin, I] = min(Z(:));
[I_row, I_col] = ind2sub(size(Z), I);

xmin = x(I_col);
ymin = y(I_row);

v = zmin + [1 2 3];

figure
contour(X,Y,Z,v)
colormap winter
hold on
plot(xmin,ymin,'o')

xlabel('x')
ylabel('y')

