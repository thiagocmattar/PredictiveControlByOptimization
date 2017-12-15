clear all; clc; 

A1 = 60^2*pi;
A2 = A1; A3 = A1; A4 = A1;

a1 = 38; a2 = a1;
a3 = 17; a4 = a3;

k1 = 4300; k2 = k1;
gama1 = 0.43; gama2 = gama1;
g = 9800;

h1_ref = 156; h2_ref = 157;
h3_ref = 234; h4_ref = 282;
u1_ref = 8; u2_ref = u1_ref;

K11 = - (a1/A1) * sqrt(g/(2*h1_ref));
K12 = (a3/A1) * sqrt(g/(2*h3_ref));
K13 = gama1*k1/A1;

K21 = - (a2/A2) * sqrt(g/(2*h2_ref));
K22 = (a4/A2) * sqrt(g/(2*h4_ref));
K23 = gama2*k2/A2;

K31 = - (a3/A3) * sqrt(g/(2*h3_ref));
K32 = (1-gama2)*k2/A3;

K41 = - (a4/A4) * sqrt(g/(2*h4_ref));
K42 = (1-gama1)*k1/A4;

s = tf('s');
g11 = K13/(s-K11);
g12 = K12*K32/((s-K11)*(s-K31));
g21 = K22*K42/((s-K21)*(s-K41));
g22 = K23/(s-K21);
G = [g11 g12; g21 g22];

K = dcgain(G);
R = transpose(inv(K));
Lambda = R.*K;
N = det(K)/(dcgain(g11)*dcgain(g22));

G2 = [g21 g22; g11 g12];
K2 = dcgain(G2);
R2 = transpose(inv(K2));
Lambda2 = R2.*K2; 
N2 = det(K2)/(dcgain(g21)*dcgain(g12));

FCC1 = -g11/g12;
FCC2 = -g22/g21;



%%Controle Preditivo------------------------------------

%Trajetória de referência
t = 1:150; wt = 150/length(t); 
nt = 1500;
w = zeros(nt,2); w(1:150,1) = wt*t; w(150:400,1) = 150; w(401:end,1) = 170; 
w(1:150,2) = w(1:150,1); w(151:250,2) = 150; w(251:end,2) = 177; 
yr = zeros(nt,2); 

%Condições iniciais
u(1:nt,1:2) = 8; ur = u;
h0(1) = 0; h0(2) = 0; h0(3) = 0; h0(4) = 0;
yr(1,1) = h0(1); yr(1,2) = h0(2);

%Horizonte de predição
nu = 10;
nr = 20; 

%Restrições da função de otimização
lb(1:15,1:2) = 0;
ub(1:15,1:2) = 50;
A = [];
b = [];
Aeq = [];
beq = [];


for i=1:(nt-nr)
    
    %Otimização da ação de controle
    ur = u(i:end,:);
    w_eval = w(i:end,:);
    fun = @(ur)system_response(ur,w_eval,h0,nu,nr);
    
    x0 = randn(((nu)),2);
    %x = fminunc(fun,x0);
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub);
    
    %Resposta do sistema para a primeira ação de controle    
    [h1,h2,h3,h4] = calculah1h2(x,h0);
    
    %Novo vetor da ação de controle
     u(i,1) = x(1,1); u(i,2) = x(1,2);
     yr(i+1,1) = h1(2); yr(i+1,2) = h2(2);
     h0(1) = h1(2); h0(2) = h2(2); h0(3) = h3(2); h0(4) = h4(2);
    
     disp(i);
    
end

time = 1:length(yr);

figure
plot(time,w(:,1),'--r',time,w(:,2),'--b');
hold
plot(time,yr(:,1),'-m',time,yr(:,2),'-g',...
    'LineWidth',2);
title('Resposta do sistema');
legend('h1_r_e_f','h2_r_e_f','h1','h2');
xlim([0 (nt-nr-1)]);
xlabel('Tempo(s)');
ylabel('hij(mm)');

figure
plot(time,u(:,1),'-r',time,u(:,2),'-b',...
    'LineWidth',2);
title('Ação de controle');
legend('u_1','u_2','h1','h2');
xlim([0 (nt-nr-1)]);
xlabel('Tempo(s)');
ylabel('Tensão nas bombas (V)');



