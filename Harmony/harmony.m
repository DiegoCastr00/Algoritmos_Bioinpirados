clc;
clear;

corridas = 1;
iteraciones = 100;

li = -10;
ls = 10;
armonias = 5;
equis = armonias;

raccept = 0.9;
rPA = 0.9;
bW = 0.1;

newX = 0;
for c = 1:corridas
    HM = randi([li ls], armonias, equis)
    fx = funcionObjetivo(HM)

    for g = 1:iteraciones
        xWorst = max(fx);
        xWorstIndex = find(fx == xWorst, 1, 'first');
        for i = 1:size(HM, 1)
            rand1 = rand(1);
            if rand1 < raccept 
                index = randi(numel(HM(i,:)));
                rand2 = rand(1);
                if rand2 < rPA
                    rand3 = -1 + 2 * rand;
                    newX(i,:) = HM(i,index) + bW * rand3;
                else
                    newX(i,:) = HM(i,index);
                end
            else 
                newX(i,:) =  randi([li ls]);
            end 

        end
        newX = transpose(newX);

        for j = 1:numel(newX)
            while newX(j) > ls || newX(j) < li
                if newX(j) > ls 
                    newX(j) = (2 * ls) - newX(j);
                elseif newX(j) < li
                    newX(j) = (2 * li) - newX(j);
                else
                    newX(j) = newX(j);
                end
            end
        end

        newXfo = funcionObjetivo(newX);
        if newXfo < xWorst
            HM(xWorstIndex, :) = newX;
            fx(xWorstIndex, :) = newXfo;
        else
            HM = HM;
            fx = fx;
        end
        newX = [];
    end
    HM
    fx
end

function fx = funcionObjetivo(poblacion)
    cuadrado = poblacion .^ 2;
    fx = sum(cuadrado, 2);
end
