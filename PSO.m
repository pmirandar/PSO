clear
close

%Función objetivo
funObj=@(xi) 3*(1-xi(1))^2*exp(-(xi(1)^2)-(xi(2)+1)^2)-10*(xi(1)/5-xi(1)^3 - xi(2)^5)*exp(-xi(1)^2-xi(2)^2)-1/3*exp(-(xi(1)+1)^2 - xi(2)^2);

%Parámetros
N = 10;
d = 2; 
lb = [-3 -3];
ub = [3 3]; 
k = 0;
kmax = 100; 
c1 = 2; 
c2 = 2;

%Inicialización de las partículas y la velocidad
for i = 1:N
    x(i,:)=rand(1,d).*(ub-lb)+lb; % Inicialización de las partículas
    v(i,:)=zeros(1,d); % Inicialización de la velocidad
end

%Evaluación de las partículas iniciales en la función objetivo
for i = 1:N
    xi=x(i,:); % Extracción de la partícula xi
    fx(i,:) = funObj(xi);
end

%Registro de la mejor partícula global y las mejores partículas locales
[gfit, ind]=min(fx); % Fitness de la mejor partícula global
g =x(ind,:); % Posición de la mejor partícula global
fp=fx; % Fitness de las mejores partículas locales
p=x; % Posición de las mejores partículas locales

% Cálculo de la superficie del espacio de búsqueda para graficar
ejex=linspace(min(lb),max(ub),50); % Vector de soluciones para d=1
ejey=ejex; % Vector de soluciones para d=2
ejez=[]; % Matriz de fitness
for i = 1:length(ejex)
    for j = 1:length(ejey)
        ejez(i,j)=funObj([ejex(i) ejey(j)]);
    end
end
[ejey,ejex] = meshgrid(ejex,ejey); % Cálculo de la malla para superficie

%PROCESO ITERATIVO
while k < kmax
    k=k+1;
    %Dibuja la superficie de búsqueda
    figure(1);
    surf(ejex,ejey,ejez) 
    hold on
    %Dibuja las partículas
    plot3(x(:,1),x(:,2),fx,'o','MarkerFaceColor','m','MarkerSize',10)
    plot3(p(:,1),p(:,2),fp,'o','MarkerFaceColor','g','MarkerSize',10) %Dibuja las mejores partículas locales(verde)
    pause(0.3)
    hold off
    %Dibuja el contorno de la superficie del espacio de búsqueda
    figure(2) % Muestra figura 2
    contour(ejex,ejey,ejez,20) 
    hold on
    % Dibujar las partículas en el contorno
    plot(x(:,1),x(:,2),'o','MarkerFaceColor','m');
    plot(p(:,1),p(:,2),'o','MarkerFaceColor','g'); %Dibuja las mejores partículas locales(verde)
    pause(0.3)
    hold off

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
        fx(i,:) = funObj(xi);
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
    Evolucion(k) = gfit;
end
%FIN DEL PROCESO ITERATIVO
figure
plot(Evolucion) % Gráfica del proceso evolutivo del PSO
disp(['Mejor solución : ', num2str(g)]) % Mejor solución, mejor posición
disp(['Mejor fitness : ', num2str(gfit)]) % Mejor fitness, mejor costo

