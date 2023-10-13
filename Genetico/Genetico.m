clc;
clear;

datos = [
    -3 3 3;
    -3 3 3;
];
individuos = 50;
iteraciones = 40;

poblacion_total = [];
real_total = [];
real_mutado_con = [];
fx_primera_generacion = [];
for i = 1:size(datos, 1)
    fila = datos(i, :);

    li = fila(1);
    ls = fila(2);
    press = fila(3);

    Nb = ceil(log((ls - li) * 10^press) + 0.9);

    poblacion = round(custom_random(Nb, individuos));

    real = binarioAReal(poblacion,individuos,li,ls);
    % disp(real);

    if isempty(poblacion_total)
        poblacion_total = poblacion;
        real_total = real;
    else
        poblacion_total = [poblacion_total, poblacion];
        real_total = [real_total, real];
    end

    % fprintf('Limite inferior: %d, Limite Superior: %d, Precisión: %d\n', li, ls, press);
    % fprintf('Nb: %d\n',Nb);
    % disp('Poblacion:');
    % disp(poblacion);
    % disp('Valores en reales:');
    % disp(real);
    % disp('--------------------------------------');
    pause(0.5);
end
% disp('Primera generacion:');
% disp(poblacion_total);
% disp('Real de la primera generacion');
% disp(real_total);
real_total_primera = real_total;

for a = 1:iteraciones
    x = real_total(:, 1);
    y = real_total(:, 2);
    
    fx = 3 * (1 - x).^2 .* exp(-x.^2 - (y + 1).^2) + ...
        10 * (x / 5 - x.^3 - y.^5) .* exp(-x.^2 - y.^2) - ...
        1/3 * exp(-((x + 1).^2) - y.^2);
    
    %disp('F(x):');
    % disp(fx)
     if a == 1
        x = real(:, 1);
        y = real(:, 1);
        fx_primera_generacion = 3 * (1 - x).^2 .* exp(-x.^2 - (y + 1).^2) + ...
            10 * (x / 5 - x.^3 - y.^5) .* exp(-x.^2 - y.^2) - ...
            1/3 * exp(-((x + 1).^2) - y.^2);
    end
    normalizado = abs(fx/(min(fx)+1));
    % disp('Normalizado:');
    % disp(normalizado)
    
    E = sum(normalizado);
    % fprintf("E: %f\n",E)
    
    p = normalizado / E;
    % disp('P:');
    % disp(p);
    
    q = zeros(size(p));
    q(1) = p(1);
    
    for i = 2:length(p)
        q(i) = q(i - 1) + p(i);
    end
    
    % disp('Q: ');
    % disp(q);
    
    
    num_iteraciones = floor(individuos/2); 
    resultados_finales = [];
    
    if rem(individuos, 2) == 1  
        num_iteraciones = num_iteraciones + 1;
    end
    
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
        
        % fprintf('Fila del padre 1 (índice %d):\n', indice_padre1);
        % disp(fila_padre1);
        
        % fprintf('Fila del padre 2 (índice %d):\n', indice_padre2);
        % disp(fila_padre2);
        
        mitad1_padre1 = fila_padre1(1:floor(end/2));
        mitad2_padre1 = fila_padre1(floor(end/2)+1:end);
        
        mitad1_padre2 = fila_padre2(1:floor(end/2));
        mitad2_padre2 = fila_padre2(floor(end/2)+1:end);
        
        hijo1 = [mitad1_padre1, mitad2_padre2];
        hijo2 = [mitad1_padre2, mitad2_padre1];
        
        % fprintf('Hijo 1:\n');
        % disp(hijo1);
        
        % fprintf('Hijo 2:\n');
        % disp(hijo2);
        
        mut1 = hijo1;
        mut2 = hijo2;
        for k = 1:length(hijo1)
            if rand() < 0.5
                mut1(k) = ~hijo1(k);
                mut2(k) = ~hijo2(k);
            end
        end
        % fprintf('Hijo 1 mutado:\n');
        % disp(mut1);
        % fprintf('Hijo 2 mutado:\n');
        % disp(mut2);
        
        resultados_finales = [resultados_finales; mut1; mut2];
        while size(resultados_finales, 1) > size(poblacion_total, 1)
            if rand() < 0.5
                resultados_finales(end,:) = [];
            else
                resultados_finales(end-1,:) = [];
            end
        end
    end
    
    % fprintf('Resultados finales:\n');
    % disp(resultados_finales);

    real_mutado = binarioAReal(resultados_finales,individuos,li,ls);
    % fprintf('Real mutado\n');
    % disp(real_mutado);

    poblacion_total = resultados_finales;
    num_cols = size(resultados_finales, 2);
    mitad_cols = num_cols / 2;
    
    resultados_finales_1 = resultados_finales(:, 1:mitad_cols);
    resultados_finales_2 = resultados_finales(:, (mitad_cols + 1):end);
    % fprintf('Resultados finales (Primera mitad):\n');
    % disp(resultados_finales_1);
    % fprintf('Resultados finales (Segunda mitad):\n');
    % disp(resultados_finales_2);
    
    Nb = size(resultados_finales_1, 2);  

    real_resultados_finales_1 = binarioAReal(resultados_finales_1,individuos,li,ls);
    
    real_resultados_finales_2 = binarioAReal(resultados_finales_2,individuos,li,ls);
    
    % fprintf('Valores reales para resultados_finales_1:\n');
    % disp(real_resultados_finales_1);
    
    % fprintf('Valores reales para resultados_finales_2:\n');
    % disp(real_resultados_finales_2);
    
    real_total =[real_resultados_finales_1,real_resultados_finales_2];
    % fprintf('Generacion %d\n', a);
    % fprintf('Real total:\n');
    % disp(real_total)
    % fprintf('-------------------------------------------\n');
    real_mutado_con = [real_mutado_con,real_mutado];

end
% disp('Fx de la primera generación:');
% disp(fx_primera_generacion);
% disp('-------------------------------------------');
% fprintf('Ultima generacion\n');
% disp(resultados_finales);
% fprintf('real de la ultima generacion\n');
% disp(real_total);
% disp('F(x):');
% disp(fx)

fprintf('maximo de la primera generacion\n');
disp(max(fx_primera_generacion))
fprintf('X1 primera generaicon: \n');
disp(max(real_total_primera(:,1)))
fprintf('X2 primera generaicon: \n');
disp(max(real_total_primera(:,2)))

fprintf('maximo de la ultima generacion\n');
disp(max(fx))
fprintf('X1 ultima generaicon: \n');
disp(max(real_total(:, 1)))
fprintf('X2 ultima generaicon: \n');
disp(max(real_total(:, 2)))

function valoresReales = binarioAReal(matrizBinaria,individuos,li,ls)
    Nb = size(matrizBinaria, 2);
    valoresReales = zeros(individuos,1);
    s = 2^Nb - 1;
    for j = 1:size(matrizBinaria, 1)
        suma_fila = 0;
        for k = 1:Nb
            bit = matrizBinaria(j, k);
            valor_bit = bit * 2^(Nb - k);
            suma_fila = suma_fila + valor_bit;
        end
        valoresReales(j) = li + (suma_fila / s) * (ls - li);
    end
end


