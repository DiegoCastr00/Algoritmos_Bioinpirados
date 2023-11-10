clc;
clear;
corridas = 30;
resultados = [];
for c = 1:corridas
    datos = [
        -3 3 3;
        -3 3 3;
    ];
    individuos =  10;
    iteraciones = 100;
    
    cr = 0.1; %cruza 
    xf=0.3;
    
    poblacion_total = [];
    real_total = [];
    con_final = [];
    for i = 1:size(datos, 1)
        fila = datos(i, :);
        li = fila(1);
        ls = fila(2);
        press = fila(3);
        Nb = ceil(log((ls - li) * 10^press) + 0.9);
        poblacion = (rand(individuos, Nb));
        poblacion_binario = round(poblacion);
        real = binarioAReal(poblacion_binario,individuos,li,ls);
        poblacion_total = [poblacion_total, poblacion];
        real_total = [real_total, real];
    end
    
    % disp('Primera generacion:');
    % disp(poblacion_total);
    % real_total
    
    
    
    final = [];
    finaldef = [];
    mejores = [];
    grafica = [];
    for g = 1:iteraciones
        sum_p = sum(poblacion_total, 2);
        for i = 1:size(poblacion_total, 1)
            numeros_aleatorios = randperm(individuos, 3);
    %         disp('    DIF1  DIF2 Mutacion');
    %         disp(numeros_aleatorios);
            
            indices_seleccionados = poblacion_total(numeros_aleatorios, :);
    %         disp('Filas correspondientes a los números aleatorios:');
    %         disp(indices_seleccionados);
            
            dif = indices_seleccionados(1,:) - indices_seleccionados(2,:);
    %         disp('Diferencia (DIF1 - DIF2):');
    %         disp(dif);
            
            peso = dif .* xf;
    %         disp('dif x factor');
    %         disp(xf);
            
            mutacion = indices_seleccionados(3,:) + peso;
    %         disp('xf + mutacion');
    %         disp(mutacion);
        
            %Manejador de cotas
            ruido = mutacion;
            mli = -0.3;
            mls = 0.3;
            
            for j = 1:numel(ruido)
                while ruido(j) >= mls || ruido(j) <= mli
                    if ruido(j) > mls 
                        ruido(j) = (2 * mls) - ruido(j);
                    elseif ruido(j) < mli
                        ruido(j) = (2 * mli) - ruido(j);
                    else
                        ruido(j) = ruido(j);
                    end
                end
        %         disp('Ruido por cotas');
        %         disp(ruido);
            end
        
    %         disp('Ruido por cotas');
    %         disp(ruido);
    %         disp('Padre ');
    %         disp(poblacion_total(i,:));
        
            %factor de cruza
            
            d = numel(ruido);
    %         fprintf('D: %d \n',d);
            
            j = randperm(d, 1);
    %         fprintf('J: %d \n\n',j);
        
            hijo = zeros(size(poblacion_total(i,:)));
    
    %%%% EXPONENCIAL        
%             for k = 1:numel(hijo)
%                 r = rand(1);
%                 if r >= cr
%                     hijo(k) = poblacion_total(i,k);
%                 elseif r <= cr 
%                     hijo(k) = ruido(k);
%                 elseif k == j
%                     hijo(k) = ruido(k);
%                 end
%             end
            
    %%%% BINOMIAL
            for k = 1:numel(hijo)
                r = rand(1);
                if k <= j && r > cr
                    hijo(k) = poblacion_total(i,k);
                elseif k <= j && r < cr 
                    hijo(k) = ruido(k);
                elseif k >= j
                    hijo(k) = ruido(k);
                end
            end
%             disp('Hijo');
%             disp(hijo);
            
            hijo_bi = round(hijo);
        
            mitad = numel(hijo) / 2; 
            parte1 = hijo(1:mitad);
            parte2 = hijo(mitad+1:end);
        
            hijo_realX = binarioAReal(parte1,1,li,ls);
            hijo_realY = binarioAReal(parte2,1,li,ls);
        
        
            real_hijo = [hijo_realX,hijo_realY];
        
            fo_hijo = fo(real_hijo);
    %         disp("FO Hijo")
    %         disp(fo_hijo)
    %         
    %         disp("X1, X2 padre")
    %         real_total(i,:)
            fo_padre = fo(real_total(i,:));
    %         disp("FO padre")
    %         disp(fo_padre)
        
            sum_hijo = sum(hijo, 2);
            sum_padre = sum_p(i,:);
        
            lideb = 0;
            lsdeb = 11;
            %DEB
            if lideb <= fo_hijo && fo_hijo <= lsdeb && lideb <= fo_padre && fo_padre <= lsdeb
                % Ambos sum_hijo y sum_padre están dentro de los límites
                if fo_hijo > fo_padre
                    final = hijo(1, :);
                    con_final = [con_final;real_hijo,fo_hijo];
                else
                    final = poblacion_total(i, :);
                    con_final = [con_final;real_total(i,:),fo_padre];
                end
            elseif (lideb <= fo_hijo && fo_hijo <= lsdeb) || (lideb <= fo_padre && fo_padre <= lsdeb)
                % Uno de los dos (ya sea sum_hijo o sum_padre) está dentro de los límites
                if lideb <= fo_hijo && fo_hijo <= lsdeb
                    final = hijo(1, :);
                    con_final = [con_final;real_hijo,fo_hijo];
                else
                    final = poblacion_total(i, :);
                    con_final = [con_final;real_total(i,:),fo_padre];
                end
            else
                % Ninguno de los dos está dentro de los límites
                distancia_hijo = min(abs(fo_hijo - lsdeb), abs(lideb - fo_hijo));
                distancia_padre = min(abs(fo_padre - lsdeb), abs(lideb - fo_padre));
                if distancia_hijo < distancia_padre
                    final = hijo(1, :);
                    con_final = [con_final;real_hijo,fo_hijo];
                else
                    final = poblacion_total(i, :);
                    con_final = [con_final;real_total(i,:),fo_padre];
                end
            end
        
            finaldef=[finaldef;final];
        end
        poblacion_total = finaldef;
        finaldef = [];
        real_total = [];
    
        real_total = con_final(:, 1:2);
        
    %     disp("nueva generacion")
    %     disp(poblacion_total)
    % 
    %     disp("nuevos reales")
    %     disp(real_total)
    
        con_final;
        num_iter = g;
        mejores=[mejores;max(con_final)];
        con_final = [];
    end
    
    mejores;
    max_mejores = max(mejores(:, 3));
    fila_maximo = find(mejores(:, 3) == max_mejores);
    fila_completa = mejores(fila_maximo, :);
    ultima_fila = fila_completa(end, :);
    resultados = [resultados;ultima_fila];
end

resultados;

resultados_ordenados = sortrows(resultados, -3)

costo = individuos * iteraciones;
fprintf("Costo: %d \n",costo);


function fx = fo(real_total)
    x = real_total(:, 1);
    y = real_total(:, 2);

    fx = 3 * (1 - x).^2 .* exp(-x.^2 - (y + 1).^2) + ...
         10 * (x / 5 - x.^3 - y.^5) .* exp(-x.^2 - y.^2) - ...
         1/3 * exp(-((x + 1).^2) - y.^2);
end


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