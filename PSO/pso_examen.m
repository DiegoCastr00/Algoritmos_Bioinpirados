clc;
clear;

corridas = 30;
iteraciones = 100;

li = -10;
ls = 10;
individuos = 10;
equis = 3;

w = 0.6;
c1 = 2;
c2 = 2.5;

mejores =[];
for c = 1:corridas
    velocidades = zeros(individuos, equis);
    poblacion = randi([li ls], individuos, equis);
    pbest = poblacion;
    fx = objetivo(poblacion);

    for g = 1:iteraciones
        minValue = min(fx);
        minIndex = find(fx == minValue, 1, 'first');
        gbest = poblacion(minIndex, :);

        for i = 1:size(poblacion, 1)
            rand1 = rand(1);
            rand2 = rand(1);
            parte1 = w * velocidades(i, :);
            parte2 = (c1 * rand1) * (pbest(i, :) - poblacion(i, :));
            parte3 = (c2 * rand2) * (gbest - poblacion(i, :));
            v = parte1 + parte2 + parte3;
            x = poblacion(i, :) + v;
            fxhijo = objetivo(x); 

            poblacion(i, :) = x;
            velocidades(i, :) = v;

            if fx(i, :) > fxhijo
               pbest(i, :) = x;
            end

            fx(i, :) = fxhijo;

        end
    end
    final = [poblacion, fx;];
    mejor = final(1, :);
    mejores = [mejores; mejor];
end
disp("Resultado de 30 corridas: ")
disp(mejores)
disp("Costo computacional: ")
costo = individuos * iteraciones;
disp(costo)
resultados_fo = final(:, end);

desviacion_estandar = std(resultados_fo)
[mejor_resultado, mejor_corrida] = min(resultados_fo);
[peor_resultado, peor_corrida] = max(resultados_fo);
corrida_promedio = mean(resultados_fo)

fila_mejor = final(mejor_corrida, :)
fila_peor = final(peor_corrida, :)


function fx = objetivo(poblacion)
    cuadrado = poblacion .^ 2;
    fx = sum(cuadrado, 2);
end


