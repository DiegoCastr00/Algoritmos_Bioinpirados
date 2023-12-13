clc;
clear;

co=1;
corridas=[];
for w=1:co
    datos = [
        -100 100 3;
        -100 100 3;
    ];
    individuos = 50;
    iteraciones = 100;
    
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
        poblacion = round(rand(individuos, Nb));
        real = binarioAReal(poblacion,individuos,li,ls);
        if isempty(poblacion_total)
            poblacion_total = poblacion;
            real_total = real;
        else
            poblacion_total = [poblacion_total, poblacion];
            real_total = [real_total, real];
        end
    end
    
    real_total_primera = real_total;
    
    for a = 1:iteraciones
        x = real_total(:, 1);
        y = real_total(:, 2);
        
        f1 = 2*x + y - 8;
        f2 = 3 * x + 3 *y.^2 - 22;
        fx = f1.^2 + f2.^2;
     
         if a == 1
            x = real(:, 1);
            y = real(:, 1);
            f1 = 2*x + y - 8;
            f2 = 3 * x + 3 *y.^2 - 22;
            fx = f1.^2 + f2.^2;
         end
         
        normalizado = abs(fx/(min(fx)+1));
        E = sum(normalizado);
        
        p = normalizado / E;
        
        q = zeros(size(p));
        q(1) = p(1);
        
        for i = 2:length(p)
            q(i) = q(i - 1) + p(i);
        end 
        
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
    
            mitad1_padre1 = fila_padre1(1:floor(end/2));
            mitad2_padre1 = fila_padre1(floor(end/2)+1:end);
            
            mitad1_padre2 = fila_padre2(1:floor(end/2));
            mitad2_padre2 = fila_padre2(floor(end/2)+1:end);
            
            hijo1 = [mitad1_padre1, mitad2_padre2];
            hijo2 = [mitad1_padre2, mitad2_padre1];
    
            mut1 = hijo1;
            mut2 = hijo2;
    
            for k = 1:length(hijo1)
                if rand() < 0.1
                    mut1(k) = ~hijo1(k);
                    mut2(k) = ~hijo2(k);
                end
            end
    
            resultados_finales = [resultados_finales; mut1; mut2];
            while size(resultados_finales, 1) > size(poblacion_total, 1)
                if rand() < 0.5
                    resultados_finales(end,:) = [];
                else
                    resultados_finales(end-1,:) = [];
                end
            end
        end
        
        real_mutado = binarioAReal(resultados_finales,individuos,li,ls);
    
        poblacion_total = resultados_finales;
        num_cols = size(resultados_finales, 2);
        mitad_cols = num_cols / 2;
        
        resultados_finales_1 = resultados_finales(:, 1:mitad_cols);
        resultados_finales_2 = resultados_finales(:, (mitad_cols + 1):end);
        
        Nb = size(resultados_finales_1, 2);  
    
        real_resultados_finales_1 = binarioAReal(resultados_finales_1,individuos,li,ls);
        
        real_resultados_finales_2 = binarioAReal(resultados_finales_2,individuos,li,ls);
    
        real_total =[real_resultados_finales_1,real_resultados_finales_2];
    
        real_mutado_con = [real_mutado_con,real_mutado];
    
    end

    
    
    max_fx_primera_generacion = max(fx_primera_generacion);
    indice_maximo_fx_primera_generacion = find(fx_primera_generacion == max_fx_primera_generacion);
    x1_primera_generacion_fx = real_total_primera(indice_maximo_fx_primera_generacion, 1);
    x2_primera_generacion_fx = real_total_primera(indice_maximo_fx_primera_generacion, 2);
%     fprintf('Maximo de la primera generación: %f\n', max_fx_primera_generacion);
%     fprintf('X1 primera generación: %f\n', x1_primera_generacion_fx);
%     fprintf('X2 primera generación: %f\n', x2_primera_generacion_fx);

    max_fx_ultima_generacion = max(fx);
    indice_maximo_fx_ultima_generacion = find(fx == max_fx_ultima_generacion);
    x1_ultima_generacion_fx = real_total(indice_maximo_fx_ultima_generacion, 1);
    x2_ultima_generacion_fx = real_total(indice_maximo_fx_ultima_generacion, 2);
%     fprintf('\n\nMáximo de la última generación: %f\n', max_fx_ultima_generacion);
%     fprintf('X1 última generación: %f\n', x1_ultima_generacion_fx);
%     fprintf('X2 última generación: %f\n', x2_ultima_generacion_fx);

   corridas = [corridas;x1_ultima_generacion_fx,x2_ultima_generacion_fx,max_fx_ultima_generacion];
end
costo = iteraciones * individuos
disp(corridas)

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


