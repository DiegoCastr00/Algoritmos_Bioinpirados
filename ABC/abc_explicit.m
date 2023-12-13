clc;
clear;

corridas = 1;
iteraciones = 200;

fuentes = 5;
equis= [
     -5 5;
     -5 5;
];
cero = 0.001;
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

        for i = 1:size(posicion, 1)
            if fx_v(i) < fx_fuente(i)
                posicion(i, :) = v(i, :);
                fx_fuente(i, :) = fx_v(i, :);
                limite(i) = 0;
            else
                limite(i) = limite(i)+ 1;  
            end
        end

        limite;
        posicion;
        fx_fuente;
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

        k = randi(fuentes, 1, fuentes);
        sigma = generarSigma(equis, fuentes);

        kaux = k.';
        dif_poblacion = newPosicionAux - (posicion(kaux, :));
        v = newPosicionAux + sigma .* dif_poblacion;

        fx_newV = funcionObjetivo(v);

        for i = 1:size(posicion, 1)
            indice_v = newfx(i,2);
            if fx_newV(i) < newfx(i,1)
                limite(indice_v,:) = 0;
                posicion(indice_v,:) = v(i,:);
                fx_fuente(indice_v,:) = fx_newV(i,:);               
            else
                limite(indice_v,:) = limite(indice_v,:) + 1;
            end
        end
        posicion;
        fx_fuente;

        for i = 1:size(posicion, 1)
            if limite(i) >= limiteMaximo
                if fx_fuente(i) < mayoreFx
                    mayorePosi = posicion(i,:);
                    mayoreFx = fx_fuente(i,:);
                    limite(i,:) = 0;
                    posicion(i,:) = mayorePosi;
                    fx_fuente(i,:) = mayoreFx;
                end
            end
        end
    end

    final= [posicion, fx_fuente];
end
final
[mejor_resultado, mejor_corrida] = min(final(:, 3));
[peor_resultado, peor_corrida] = max(final(:, 3));
fila_mejor = final(mejor_corrida, :)
fila_peor = final(peor_corrida, :)


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
    cuadrado = poblacion .^ 2;
    fx = sum(cuadrado, 2);
end
