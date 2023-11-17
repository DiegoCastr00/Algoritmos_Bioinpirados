clc;
clear;

corridas = 5;
iteraciones = 100;

li = -5;
ls = 5;
individuos = 6;
equis = 3;

w= 0.5;
c1 = 2;
c2 = 1.5;

for c = 1:corridas
    velocidades = zeros(individuos,equis);
    poblacion = randi([li ls],individuos,equis);
    pbest = poblacion;
    fx = fo(poblacion);

    for g = 1:iteraciones
        minValue = min(fx);
        minIndex = find(fx == minValue, 1, 'first');
        gbest = poblacion(minIndex, :);

        for i = 1:size(poblacion, 1)
            rand1 = rand(1);
            rand2 = rand(1);
            parte1 = w * velocidades(i,:);
            parte2 = (c1*rand1)*((pbest(i,:))-(poblacion(i,:)));
            parte3 = (c2*rand2)*((gbest)-(poblacion(i,:)));
            v = parte1 + parte2 + parte3;
            x = poblacion(i,:) + v;
            fxhijo = fo(x);

            poblacion(i,:) = x;
            velocidades(i,:) = v;
            
            if fx(i,:) > fxhijo
               pbest(i,:) = x;
            end
            
            fx(i,:) = fxhijo;
           

        end
    end
    final = [poblacion, fx;]
end

function fx = fo(poblacion)
    cuadrado = poblacion .^ 2;
    fx = sum(cuadrado, 2);
end

