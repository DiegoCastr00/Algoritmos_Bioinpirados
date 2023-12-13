clc;
clear;
rng('shuffle')
corridas = 1;
iteraciones = 100000;

equis= [
     0 1200;
     0 1200;
    -0.55 0.55;
    -0.55 0.55;
];
armonias = 50;

raccept = 0.9;
rPA = 0.3;
bW = 0.085;

cero_gordo =0.00001;
todos=[];
resultados = [];
for c = 1:corridas
    HM = [];
    for i = 1:size(equis, 1)
        fila = equis(i, :);
        li = fila(1);
        ls = fila(2);
        poblacion= (ls - li) * rand(1, armonias) + li;
        poblacion = transpose(poblacion);
        HM = [HM, poblacion];
    end
    HM;
    fx_padre = funcionObjetivo(HM);

    for g = 1:iteraciones
        Worst = max(fx_padre);
        WorstIndex = find(fx_padre == Worst, 1, 'first');

        [g1_padre,g2_padre] = desigualdad(HM(WorstIndex, :));
        [h3_padre, h4_padre, h5_padre] = igualdad(HM(WorstIndex, :));

%         padre_violaciones = [g1_padre,g2_padre,h3_padre, h4_padre, h5_padre];

        
        for i = 1:size(HM, 2)
            rand1 = rand(1);
            if rand1 < raccept 
                index = randi(numel(HM(i,:)));
                rand2 = rand(1);
                if rand2 < rPA
                    rand3 = -1 + 2 * rand;
                    HM_hijo(:,i) = HM(index,i) + (bW * rand3);
                else
                    HM_hijo(:,i) = HM(index,i);
                end
            else   
                HM_hijo(:,i) =  (equis(i,2) - equis(i,1)) * rand + equis(i,1);
            end 
        end
        %HM_hijo = transpose(HM_hijo)
        HM_hijo;
        for j = 1:numel(HM_hijo)
            fila = equis(j, :);
            li = fila(1);
            ls = fila(2);
            while HM_hijo(j) > ls || HM_hijo(j) < li
                if HM_hijo(j) > ls 
                    HM_hijo(j) = (2 * ls) - HM_hijo(j);
                elseif HM_hijo(j) < li
                    HM_hijo(j) = (2 * li) - HM_hijo(j);
                else
                    HM_hijo(j) = HM_hijo(j);
                end
            end
        end
        HM_hijo;
        fx_hijo = funcionObjetivo(HM_hijo);
        [g1_hijo,g2_hijo] = desigualdad(HM_hijo);
        [h3_hijo, h4_hijo, h5_hijo] = igualdad(HM_hijo);

        h3_padre(h3_padre < cero_gordo) = 0;
        h4_padre(h4_padre < cero_gordo) = 0;
        h5_padre(h5_padre < cero_gordo) = 0;


        h3_hijo(h3_hijo < cero_gordo) = 0;
        h4_hijo(h4_hijo < cero_gordo) = 0;
        h5_hijo(h5_hijo < cero_gordo) = 0;

%         hjo_violaciones = [g1_hijo,g2_hijo,h3_hijo, h4_hijo, h5_hijo];

        % Restricción 3

        SVR_padre = max(0,g1_padre) + max(0,g2_padre) + abs(h3_padre) + abs(h4_padre) + abs(h5_padre);
        SVR_hijo = max(0,g1_hijo) + max(0,g2_hijo) + abs(h3_hijo) + abs(h4_hijo) + abs(h5_hijo);

        % if SVR_hijo < SVR_padre
        %     HM(WorstIndex, :) = HM_hijo;
        %     fx_padre(WorstIndex, :) = fx_hijo;
        % else
        %     HM = HM;
        %     fx_padre = fx_padre;
        % end

        % Reglas de DEB
        % factibilidad del padre
        if g1_padre <= 0 && g2_padre <= 0 && h3_padre == 0 && h4_padre == 0 && h5_padre == 0
            if g1_hijo <= 0 && g2_hijo <= 0 && h3_hijo == 0 && h4_hijo == 0 && h5_hijo == 0
                % factibilidad del hijo
                % Regla 1: Entre dos soluciones factibles, se elige la de menor valor objetivo
                if fx_hijo < Worst
                    HM(WorstIndex, :) = HM_hijo;
                    fx_padre(WorstIndex, :) = fx_hijo;
%                     fprintf("regla 1, hijo fue menor fo\n")
                else
%                     fprintf("regla 1, padre fue menor fo\n");
                end

            elseif g1_hijo > 0 || g2_hijo > 0 || h3_hijo ~= 0 || h4_hijo ~= 0 || h5_hijo ~= 0
                % Hijo es infactible
                % Regla 2: Entre una solución factible y otra infactible, se prefiere la factible
                HM = HM;
                fx_padre = fx_padre;
%                 fprintf("regla 2 el padre es factible")

            else
                % Ambos son infactibles
                % Regla 3: Entre dos soluciones infactibles, se elige la que tenga menor suma de violaciones
                if SVR_hijo < SVR_padre
                    HM(WorstIndex, :) = HM_hijo;
                    fx_padre(WorstIndex, :) = fx_hijo;
%                     fprintf("regla 3, hijo fue menor SVR\n")
                else 
%                     fprintf("regla 3, padre tuvo menor SVR\n")
                end
            end

        elseif g1_hijo <= 0 && g2_hijo <= 0 && h3_hijo == 0 && h4_hijo == 0 && h5_hijo == 0
            % Padre es infactible, pero hijo es factible
            HM(WorstIndex, :) = HM_hijo;
            fx_padre(WorstIndex, :) = fx_hijo;
%             fprintf("regla 2, hijo es factible") % Reemplazar el peor del padre con el hijo
        else
            % Ambos son infactibles
            % Regla 3: Entre dos soluciones infactibles, se elige la que tenga menor suma de violaciones
            if SVR_hijo < SVR_padre
                HM(WorstIndex, :) = HM_hijo;
                fx_padre(WorstIndex, :) = fx_hijo;
%                 fprintf("regla 3, hijo tuvo menor SVR\n")
            else 
%                 fprintf("regla 3, padre tuvo menor SVR\n")
            end
        end
        HM_hijo = [];
        
        HM;
        fx_padre;
    end
    HM;
    fx_padre;
    [g1,g2] = desigualdad(HM);
    [h3, h4, h5] = igualdad(HM);

%     final_violaciones = [g1,g2,h3, h4, h5];
    h3(h3 < cero_gordo) = 0;
    h4(h4 < cero_gordo) = 0;
    h5(h5 < cero_gordo) = 0;

    SVR_final = max(0,g1) + max(0,g2) + abs(h3) + abs(h4) + abs(h5);

    best = min(SVR_final);
    bestIndex = find(SVR_final == best, 1, 'first');
    resultados = [resultados;HM(bestIndex, :), fx_padre(bestIndex), best];
    todos = [todos;HM, fx_padre,SVR_final]
    disp(c);
end

disp("Resultado de 30 corridas: "), 
disp("      x1        x2        x3        x4        fx        SVR")
resultados;
resultados = sortrows(resultados, size(resultados, 2));
disp(resultados)

disp("Costo computacional: ")
costo = armonias * iteraciones;
disp(costo)

desviacion_estandar = std(resultados(:, 5))
corrida_promedio = mean(resultados(:, 5))

[mejor_resultado, mejor_corrida] = min(resultados(:, 6));
[peor_resultado, peor_corrida] = max(resultados(:, 6));
fila_mejor = resultados(mejor_corrida, :)
fila_peor = resultados(peor_corrida, :)


function fx = funcionObjetivo(HM)
    x1 = HM(:, 1);
    x2 = HM(:, 2);
    fx = 3*x1 + 0.000001*x1.^3 + 2*x2 + (0.000002/3)*x2.^3;
end

function [g1,g2] = desigualdad(HM)
    x4 = HM(:, 4);
    x3 = HM(:, 3);

    g1 = -(x4) + (x3) - 0.55;
    g2 = -(x3) + (x4) - 0.55;
end

function [h3, h4, h5] = igualdad(HM)
    x1 = HM(:, 1);
    x2 = HM(:, 2);
    x3 = HM(:, 3);
    x4 = HM(:, 4);

    h3 = 1000*(sin(-(x3) - 0.25)) + 1000*(sin(-(x4) - 0.25)) + 894.8 - (x1);
    h4 = 1000*(sin((x3) - 0.25)) + 1000*(sin((x3) - (x4) - 0.25)) + 894.8 - (x2);
    h5 = 1000*(sin((x4) - 0.25)) + 1000*(sin((x4) - (x3) - 0.25)) + 1294.8;
end