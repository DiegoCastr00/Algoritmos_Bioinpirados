clc;
clear;

corridas = 1;
iteraciones = 100000;

fuentes = 50;

equis= [
     0 1200;
     0 1200;
    -0.55 0.55;
    -0.55 0.55;
];

cero_gordo = 0.000001;
limiteMaximo = fuentes * size(equis, 1);
resultados = [];

for c = 1:corridas
    sigma = generarSigma(equis, fuentes);
    posicion = generarFuentes(equis,fuentes);
    limite = zeros(size(posicion, 1), 1);

    fx_fuente = funcionObjetivo(posicion);
    mayorePosi = [];
    mayoreFx = min(fx_fuente);

    for g = 1:iteraciones
        k = randi(fuentes, 1, fuentes);  %No unicas
        %k = randperm(fuentes); %Unicos
        v = nuevaSolucion(posicion,k,sigma);
        fx_v = funcionObjetivo(v);

        [g1_padre,g2_padre] = desigualdad(posicion);
        [h3_padre, h4_padre, h5_padre] = igualdad(posicion);

        h3_padre(h3_padre < cero_gordo) = 0;
        h4_padre(h4_padre < cero_gordo) = 0;
        h5_padre(h5_padre < cero_gordo) = 0;


        [g1_hijo,g2_hijo] = desigualdad(v);
        [h3_hijo, h4_hijo, h5_hijo] = igualdad(v);
        h3_hijo(h3_hijo < cero_gordo) = 0;
        h4_hijo(h4_hijo < cero_gordo) = 0;
        h5_hijo(h5_hijo < cero_gordo) = 0;

        SVR_padre = max(0,g1_padre) + max(0,g2_padre) + abs(h3_padre) + abs(h4_padre) + abs(h5_padre);
        SVR_hijo = max(0,g1_hijo) + max(0,g2_hijo) + abs(h3_hijo) + abs(h4_hijo) + abs(h5_hijo);

        for i = 1:size(posicion, 1)
%             SVR_hijo(i)
%             SVR_padre(i)
        
            switch true
                case SVR_hijo(i) < SVR_padre(i)
                    posicion(i, :) = v(i, :);
                    fx_fuente(i, :) = fx_v(i, :);
                    limite(i) = 0;
        
                case SVR_hijo(i) > SVR_padre(i)
                    limite(i) = limite(i) + 1;
        
                case SVR_hijo(i) == 0 && SVR_padre(i) == 0

        
                    if fx_v(i) < fx_fuente(i)
                        posicion(i, :) = v(i, :);
                        fx_fuente(i, :) = fx_v(i, :);
                        limite(i) = 0;
                    else
                        limite(i) = limite(i) + 1;
                    end
        
                case SVR_hijo(i) == SVR_padre(i)
                    posicion(i, :) = v(i, :);
                    fx_fuente(i, :) = fx_v(i, :);
                    limite(i) = 0;
        
                case SVR_hijo(i) == 0 && SVR_padre(i) ~= 0
                    posicion(i, :) = v(i, :);
                    fx_fuente(i, :) = fx_v(i, :);
                    limite(i) = 0;
            end
        end

        
        posicion;
        fx_fuente;
        limite;
        num_mejores = floor(size(posicion, 1)/2);
        
        [~, ordenado] = sort(fx_fuente);
        indices = ordenado(1:num_mejores);
        mejores = posicion(indices, :);
        mejores = [mejores, indices];
        mejores_fx = fx_fuente(indices, :);
        mejores_fx = [mejores_fx, indices];

        newPosicion = repmat(mejores, ceil(size(posicion, 1) / num_mejores), 1);
        newPosicion = newPosicion(1:size(posicion, 1), :);

        newfx = repmat(mejores_fx, ceil(size(fx_fuente, 1) / num_mejores), 1);
        newfx = newfx(1:size(posicion, 1), :);

        newPosicionAux = newPosicion(:, 1:end-1);
        newPosicionAux = newPosicionAux(1:size(posicion, 1), :);

        %k = randi(fuentes, 1, fuentes);  %No unicas
        k = randperm(fuentes); %Unicos
        sigma = generarSigma(equis, fuentes);

        kaux = k.';
        dif_poblacion = newPosicionAux - (posicion(kaux, :));
        v = newPosicionAux + sigma .* dif_poblacion;
        fx_newV = funcionObjetivo(v);


        [g1_padre,g2_padre] = desigualdad(newPosicion);
        [h3_padre, h4_padre, h5_padre] = igualdad(newPosicion);

        h3_padre(h3_padre < cero_gordo) = 0;
        h4_padre(h4_padre < cero_gordo) = 0;
        h5_padre(h5_padre < cero_gordo) = 0;


        [g1_hijo,g2_hijo] = desigualdad(v);
        [h3_hijo, h4_hijo, h5_hijo] = igualdad(v);
        h3_hijo(h3_hijo < cero_gordo) = 0;
        h4_hijo(h4_hijo < cero_gordo) = 0;
        h5_hijo(h5_hijo < cero_gordo) = 0;

        SVR_padre = max(0,g1_padre) + max(0,g2_padre) + abs(h3_padre) + abs(h4_padre) + abs(h5_padre);
        SVR_hijo = max(0,g1_hijo) + max(0,g2_hijo) + abs(h3_hijo) + abs(h4_hijo) + abs(h5_hijo);

        for i = 1:size(posicion, 1)
            indice_v = newfx(i,2);
            switch true
                case SVR_hijo(i) < SVR_padre(i)
                    posicion(indice_v,:) = v(i, :);
                    fx_fuente(indice_v,:) = fx_newV(i,:) ;
                    limite(indice_v,:) = 0;
        
                case SVR_hijo(i) > SVR_padre(i)
                    limite(indice_v,:) = limite(indice_v,:) + 1;
        
                case SVR_hijo(i) == 0 && SVR_padre(i) == 0
                    if fx_newV(i,:) < newfx(i,1)
                        posicion(indice_v,:) = v(i, :);
                        fx_fuente(indice_v,:) = fx_newV(i,:) ;
                        limite(indice_v,:) = 0;
                    else
                        limite(indice_v,:) = limite(indice_v,:) + 1;
                    end
        
                case SVR_hijo(i) == SVR_padre(i)
                    posicion(indice_v,:) = v(i, :);
                    fx_fuente(indice_v,:) = fx_newV(i,:) ;
                    limite(indice_v,:) = 0;
        
                case SVR_hijo(i) == 0 && SVR_padre(i) ~= 0
                    posicion(indice_v,:) = v(i, :);
                    fx_fuente(indice_v,:) = fx_newV(i,:) ;
                    limite(indice_v,:) = 0;
            end
        end

        for i = 1:size(posicion, 1)
            if limite(i) >= limiteMaximo
                resultados = [resultados; posicion(i,:), fx_fuente(i,:)];
                limite(i,:) = 0;
                posicion(i,:) =  generarFuentes(equis,1);
                fx_fuente(i,:) = funcionObjetivo(posicion(i,:));
            end
        end

    end
    final= [posicion, fx_fuente];
end

[g1,g2] = desigualdad(posicion);
[h3, h4, h5] = igualdad(posicion);
h3(h3 < cero_gordo) = 0;
h4(h4 < cero_gordo) = 0;
h5(h5 < cero_gordo) = 0;
SVR_final = max(0,g1) + max(0,g2) + abs(h3) + abs(h4) + abs(h5);
final = [final,SVR_final];


final_filtrado = final(final(:, 6) == 0, :);
[mejor_resultado, mejor_corrida] = min(final_filtrado(:, 5));
[peor_resultado, peor_corrida] = max(final_filtrado(:, 5));
fila_mejor = final_filtrado(mejor_corrida, :);
fila_peor = final_filtrado(peor_corrida, :);


[g1,g2] = desigualdad(resultados(:, 1:end-1));
[h3, h4, h5] = igualdad(resultados(:, 1:end-1));
h3(h3 < cero_gordo) = 0;
h4(h4 < cero_gordo) = 0;
h5(h5 < cero_gordo) = 0;
SVR_final = max(0,g1) + max(0,g2) + abs(h3) + abs(h4) + abs(h5);

resultados = [resultados,SVR_final]
final_filtrado = resultados(resultados(:, 6) == 0, :);
[mejor_resultado, mejor_corrida] = min(final_filtrado(:, 5));
[peor_resultado, peor_corrida] = max(final_filtrado(:, 5));
fila_mejor = final_filtrado(mejor_corrida, :)
fila_peor = final_filtrado(peor_corrida, :)


function posicion = generarFuentes(equis, fuentes)
    posicion = [];
    for i = 1:size(equis, 1)
        fila = equis(i, :);
        li = fila(1);
        ls = fila(2);
        poblacion = (ls - li) * rand(1, fuentes) + li;
        poblacion = transpose(poblacion);
        posicion = [posicion, poblacion];
    end
end


function sig = generarSigma(equis, fuentes)
    sig = [];
    for i = 1:size(equis, 1)
        sigma = (1 - (-1)) * rand(1, fuentes) - 1;
        sigma = transpose(sigma);
        sig = [sig,sigma];
    end
end

function v = nuevaSolucion(poblacion,k,sigma)
    kaux = k.';
    dif_poblacion = poblacion - poblacion(kaux, :);
    v = poblacion + sigma .* dif_poblacion;
end

function fx = funcionObjetivo(poblacion)
    x1 = poblacion(:, 1);
    x2 = poblacion(:, 2);
    fx = 3*x1 + 0.000001*x1.^3 + 2*x2 + (0.000002/3)*x2.^3;
end

function [g1,g2] = desigualdad(poblacion)
    x4 = poblacion(:, 4);
    x3 = poblacion(:, 3);

    g1 = -(x4) + (x3) - 0.55;
    g2 = -(x3) + (x4) - 0.55;
end

function [h3, h4, h5] = igualdad(poblacion)
    x1 = poblacion(:, 1);
    x2 = poblacion(:, 2);
    x3 = poblacion(:, 3);
    x4 = poblacion(:, 4);

    h3 = 1000*(sin(-(x3) - 0.25)) + 1000*(sin(-(x4) - 0.25)) + 894.8 - (x1);
    h4 = 1000*(sin((x3) - 0.25)) + 1000*(sin((x3) - (x4) - 0.25)) + 894.8 - (x2);
    h5 = 1000*(sin((x4) - 0.25)) + 1000*(sin((x4) - (x3) - 0.25)) + 1294.8;
end

