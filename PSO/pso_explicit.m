clc;
clear;

corridas = 1;
iteraciones = 100;

li = -5;
ls = 5;
individuos = 4;
equis = 3;

w= 0.1;
c1 = 1.5;
c2 = 1.5;

velocidades = zeros(individuos,equis);

for c = 1:corridas
   poblacion = randi([li ls],individuos,equis);
    pbest = poblacion;
    fx = fo(poblacion);
    for g = 1:iteraciones
        minValue = min(fx);
        minIndex = find(fx == minValue, 1, 'first');
        gbest = poblacion(minIndex, :);

        disp("========  ENTRADA ===========")
        velocidades;
        poblacion;
        fx;
        gbest;
        disp("========  INICIO ===========")
        
        for i = 1:size(poblacion, 1)
            rand1 = rand(1)
            rand2 = rand(1)
            parte1 = w * velocidades(i,:)
            parte2 = (c1*rand1)*((pbest(i,:))-(poblacion(i,:)))
            parte3 = (c2*rand2)*((gbest)-(poblacion(i,:)))
            v = parte1 + parte2 + parte3
            x = poblacion(i,:) + v
            fxhijo = fo(x)
            poblacion(i,:) = x
            velocidades(i,:) = v
            
            if fx(i,:) > fxhijo
               pbest(i,:) = x
            end

            fx(i,:) = fxhijo
           
            disp("~~~~~~~~~~")

        end
    end
end

function fx = fo(poblacion)
    cuadrado = poblacion .^ 2;
    fx = sum(cuadrado, 2);
end

