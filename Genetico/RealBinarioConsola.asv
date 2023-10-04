clc;
clear; 
datos = [
    -3 3 1;
    -3 3 1;
];

n = 10;

poblacion_total = [];
real_total = [];

for i = 1:size(datos, 1)
    fila = datos(i, :);

    li = fila(1);
    ls = fila(2);
    press = fila(3);

    Nb = ceil(log((ls - li) * 10^press) + 0.9);

    poblacion = round(custom_random(Nb, n));

    s = 2^Nb - 1;

    real = zeros(n, 1);
    for j = 1:n
        suma_fila = 0;
        for k = 1:Nb
            bit = poblacion(j, k);
            valor_bit = bit * 2^(Nb - k);
            suma_fila = suma_fila + valor_bit;
        end
        real(j) = li + (suma_fila / s) * (ls - li);
    end

    if isempty(poblacion_total)
        poblacion_total = poblacion;
        real_total = real;
    else
        poblacion_total = [poblacion_total, poblacion];
        real_total = [real_total, real];
    end

    fprintf('Limite inferior: %d, Limite Superior: %d, Precisión: %d\n', li, ls, press);
    fprintf('Nb: %d\n',Nb);
    disp('Poblacion:');
    disp(poblacion);

    disp('Valores en reales:');
    disp(real);
    disp('--------------------------------------');
    pause(0.5);
end
disp('Todas las poblaciones:');
disp(poblacion_total);
disp('Todos los reales:');
disp(real_total);

x = real_total(:, 1);
y = real_total(:, 2);

fx = 3 * (1 - x).^2 .* exp(-x.^2 - (y + 1).^2) + ...
    10 * (x / 5 - x.^3 - y.^5) .* exp(-x.^2 - y.^2) - ...
    1/3 * exp(-((x + 1).^2) - y.^2);

disp('F(x):');
disp(fx)

normalizado = abs(fx/(min(fx)+1));
disp('Normalizado:');
disp(normalizado)

E = sum(normalizado);
fprintf("E: %f",E)

p = normalizado / E;
disp('P:');
disp(p);

q = zeros(size(p));
q(1) = p(1);

for i = 2:length(p)
    q(i) = q(i - 1) + p(i);
end

disp('Q: ');
disp(q);


num_iteraciones = n / 2; 
resultados_finales = [];

for i = 1:num_iteraciones

    indice_padre1 = 0;
    indice_padre2 = 0;
    while indice_padre1 == indice_padre2
        numero_aleatorio1 = rand();
        numero_aleatorio2 = rand();
        indice_padre1 = find(numero_aleatorio1 <= q, 1);
        indice_padre2 = find(numero_aleatorio2 <= q, 1);
    end
    
    fila_padre1 = poblacion_total(indice_padre1, :);
    fila_padre2 = poblacion_total(indice_padre2, :);
    
    fprintf('Fila del padre 1 (índice %d):\n', indice_padre1);
    disp(fila_padre1);
    
    fprintf('Fila del padre 2 (índice %d):\n', indice_padre2);
    disp(fila_padre2);
    
    mitad1_padre1 = fila_padre1(1:floor(end/2));
    mitad2_padre1 = fila_padre1(floor(end/2)+1:end);
    
    mitad1_padre2 = fila_padre2(1:floor(end/2));
    mitad2_padre2 = fila_padre2(floor(end/2)+1:end);
    
    hijo1 = [mitad1_padre1, mitad2_padre2];
    hijo2 = [mitad1_padre2, mitad2_padre1];
    
    fprintf('Hijo 1:\n');
    disp(hijo1);
    
    fprintf('Hijo 2:\n');
    disp(hijo2);
    
    mut1 = hijo1;
    mut2 = hijo2;
    for k = 1:length(hijo1)
        if rand() < 0.5
            mut1(k) = ~hijo1(k);
            mut2(k) = ~hijo2(k);
        end
    end
    fprintf('Hijo 1 mutado:\n');
    disp(mut1);
    fprintf('Hijo 2 mutado:\n');
    disp(mut2);
    
    resultados_finales = [resultados_finales; mut1; mut2];
end

fprintf('Resultados finales:\n');
disp(resultados_finales);


