clc;
clear;

tamX = 50;
tamY = 50;

%mundo = randi([0, 1], tamX, tamY);

mundo = zeros(tamX, tamY);
%Blinker
mundo(5,3) = 1;
mundo(5,4) = 1;
mundo(5,5) = 1;
%Glider
mundo(21,21) = 1;
mundo(22,22) = 1;
mundo(22,23) = 1;
mundo(21,23) = 1;
mundo(20,23) = 1;

figura = figure;
set(figura, 'Position', [100, 100, tamX * 8, tamY * 8]); 

umbral = 150;
generacion = 0;
num_celulas_vivas_anterior = 0;
contador_sin_cambios = 0;

while ishandle(figura)
    nuevo_mundo = toroidal(mundo);

    clf;
    imagesc(mundo);
    colormap([0,0,0; 1,1,1]);
    axis off;

    num_celulas_vivas = sum(nuevo_mundo(:));

    text(1, tamY + 5, ['Generacion: ' num2str(generacion)], 'Color', 'k', 'FontSize', tamX * 0.2);
    text(tamX / 2, tamY + 5, ['Celulas vivas: ' num2str(num_celulas_vivas)], 'Color', 'k', 'FontSize', tamX * 0.2);

    if num_celulas_vivas == num_celulas_vivas_anterior
        contador_sin_cambios = contador_sin_cambios + 1;
    else
        contador_sin_cambios = 0;
    end
    
    if contador_sin_cambios >= umbral
        break;
    end

    num_celulas_vivas_anterior = num_celulas_vivas;
    generacion = generacion + 1;

    pause(0.01);
    mundo = nuevo_mundo;
end

function nuevo_mundo = toroidal(mundo)
    [tamX, tamY] = size(mundo);
    nuevo_mundo = mundo;

    for y = 1:tamY
        for x = 1:tamX
            izquierda = mod(x - 2, tamX) + 1;
            derecha = mod(x, tamX) + 1;
            arriba = mod(y - 2, tamY) + 1;
            abajo = mod(y, tamY) + 1;

            vecinos = mundo(izquierda, arriba) + ...
                     mundo(x, arriba) + ...
                     mundo(derecha, arriba) + ...
                     mundo(izquierda, y) + ...
                     mundo(derecha, y) + ...
                     mundo(izquierda, abajo) + ...
                     mundo(x, abajo) + ...
                     mundo(derecha, abajo);

            switch mundo(x, y)
                case 0
                    if vecinos == 3
                        nuevo_mundo(x, y) = 1;
                    end
                case 1
                    if vecinos < 2 || vecinos > 3
                        nuevo_mundo(x, y) = 0;
                    end
            end
        end
    end
end


function nuevo_mundo = bordes(mundo)
    [tamX, tamY] = size(mundo);
    nuevo_mundo = mundo;

    for y = 1:tamY
        for x = 1:tamX
            izquierda = max(x - 1, 1);
            derecha = min(x + 1, tamX);
            arriba = max(y - 1, 1);
            abajo = min(y + 1, tamY);

            vecinos = mundo(izquierda, arriba) + ...
                     mundo(x, arriba) + ...
                     mundo(derecha, arriba) + ...
                     mundo(izquierda, y) + ...
                     mundo(derecha, y) + ...
                     mundo(izquierda, abajo) + ...
                     mundo(x, abajo) + ...
                     mundo(derecha, abajo);

            switch mundo(x, y)
                case 0
                    if vecinos == 3
                        nuevo_mundo(x, y) = 1;
                    end
                case 1
                    if vecinos < 2 || vecinos > 3
                        nuevo_mundo(x, y) = 0;
                    end
            end
        end
    end
end
