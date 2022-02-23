function [g,gfit] = PSO(Pa,lb,ub)

k = 0;
N = Pa(1);
d = Pa(2);
kmax = Pa(3);
c1 = Pa(4);
c2 = Pa(5);

%Inicialización de las partículas y la velocidad
for i = 1:N
  %  rng('multFibonacci');
    x(i,:)=rand(1,d).*(ub-lb)+lb; % Inicialización de las partículas
    v(i,:)=zeros(1,d); % Inicialización de la velocidad
end

%Evaluación de las partículas iniciales en la función objetivo
for i = 1:N
    xi=x(i,:); % Extracción de la partícula xi
       fx(i,:) = chi2(xi);
end

%Registro de la mejor partícula global y las mejores partículas locales
[gfit, ind]=min(fx); % Fitness de la mejor partícula global
g =x(ind,:); % Posición de la mejor partícula global
fp=fx; % Fitness de las mejores partículas locales
p=x; % Posición de las mejores partículas locales

%PROCESO ITERATIVO
while k < kmax
    k=k+1;

    %Cálculo de la nueva velocidad para cada partícula
    for i= 1:N
        xi = x(i,:);
        pi = p(i,:);
        v(i,:) = v(i,:)+c1*rand(1,d).*(pi-xi)+c2*rand(1,d).*(g-xi);
    end
    %Cálculo de la nueva posición de cada partícula
    x=x+v;
    % Verificar que las partículas no se salgan de los límites
    for i = 1:N 
        for j=1:d 
            if x(i,j) < lb(j)
                x(i,j) = lb(j);
            elseif x(i,j) > ub(j)
                x(i,j) = ub(j);
            end
        end
    end
    %Evaluación de las nuevas partículas en la función objetivo
    for i = 1:N
        xi = x(i,:);
        fx(i,:) = chi2(xi);
    end
    % Registro de la mejor partícula global y las mejores locales
    [gfitkmas1, ind] = min(fx);
    if gfitkmas1 < gfit
        gfit = gfitkmas1;
        g = x(ind,:);
    end
    %Actualiza el fitness de las mejores partículas locales
    for i = 1:N
        if fx(i,:) < fp(i,:)
            fp(i,:) = fx(i,:);
            p(i,:) = x(i,:);
        end
    end
    % Registro de las mejores soluciones encontradas en cada generación
%    Evolucion(k) = gfit;
end
%FIN DEL PROCESO ITERATIVO
