clear
close

%ParÃ¡metros
N = 10000;
d = 6; 
k = 0;
kmax = 150; 
c1 = 0.9; 
c2 = 0.9;

Par = [N d kmax c1 c2];

global mu;
global md;
global ms;
global mc;
global mb;
global mt;

global vcb_th;
global vus_th;
global vub_th; 
global jarslkog_th;

global sigma_vcb;
global sigma_vus;
global sigma_vub;
global sigma_jars;

mu = 1.23;
md = 2.67;
ms = 53.16;
mc = 620;
mb = 2839;
mt = 168260;

%%/*Valores teóricos de los elementos de la matriz CKM*/
vcb_th = 0.04053;
vus_th = 0.22650;
vub_th = 0.00361;
jarslkog_th = 0.0000300;

%%/*donde la incertidumbre de cada uno es*/
sigma_vcb = 0.00061;
sigma_vus = 0.00048;
sigma_vub = 0.00009;
sigma_jars = 0.00009;

lb = [mc mu ms md 0 0];
ub = [mt mc mb ms 2*pi 2*pi];

Sols = 3; %% Numero de soluciones

%% Para probar soluciones
%% xi=[33389.8248 506.802639 53.16 15.9640073 0 6.28318531];
%% disp(['Comprobando el chi^2:', num2str(chi2(xi))]);

Data = fopen('run.dat','a');

for j = 1:Sols
   [g,gfit] = PSO(Par,lb,ub);
   disp(['Solution %d',j]);
   disp(j);
   disp(['Mejor soluciÃ³n : ', num2str(g)]) % Mejor soluciÃ³n, mejor posiciÃ³n
   disp(['Mejor fitness : ', num2str(gfit)]) % Mejor fitness, mejor costo
   
   if gfit < 4
       fprintf(Data,'%2.7f ',gfit);
       %%              Au  Eu  Ad  Ed  Phi1 Phi2
       fprintf(Data,'%7.7f %7.7f %7.7f %7.7f %7.7f %7.7f\n',g);
   end
end
fclose(Data);
