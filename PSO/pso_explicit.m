clc;
clear;

corridas = 1;
iteraciones = 1;
li = -5;
ls = 5;
individuos = 4;
equis = 3;

velocidades = zeros(individuos,equis)


for c = 1:corridas
%    poblacion = randi([li ls],individuos,equis);
%    disp("poblacion: ")
%    disp(poblacion)
    poblacion = [-5, 1, 5;
                  2, -3, 1;
                  1, 1, 1;
                  -5 5 5];
    pbest = poblacion
    for g = 1:iteraciones
        fx = fo(poblacion);
        pobla = [poblacion, fx]
        minValue = min(fx);
        minIndex = find(fx == minValue);
        gbest = pobla(minIndex, :);
        
        disp("gbest: ")
        disp(gbest)
        for i = 1:size(poblacion, 1)
            
            
%           fo_padre = fo(fx(i,:))
        end
    end
end

function fx = fo(poblacion)
    cuadrado = poblacion .^ 2;
    fx = sum(cuadrado, 2);
end
