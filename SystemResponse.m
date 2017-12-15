function [h1,h2,h3,h4] = calculah1h2(u,h0)

%Parâmetros do sistema
A1 = 60^2*pi;
A2 = A1; A3 = A1; A4 = A1;
a1 = 38; a2 = a1;
a3 = 17; a4 = a3;
k1 = 4300; k2 = k1;
gama1 = 0.43; gama2 = gama1;
g = 9800;
arg = [a3/A3 (1-gama2)*k2/A3 a4/A4 (1-gama1)*k1/A4 a1/A1 a3/A1 gama1*k1/A1 a2/A2 a4/A2 gama2*k2/A2 g];

%Horizonte de predição
n=length(u);

%Sinal de controle livre e forçado
v1 = u(:,1); v2 = u(:,2);

%Resposta do sistema
h1 = zeros(n,1); h2 = zeros(n,1); h3 = zeros(n,1); h4 = zeros(n,1);
h1(1) = h0(1); h2(1) = h0(2); h3(1) = h0(3); h4(1) = h0(4);
for k=1:(n-1)
    h3(k+1) = h3(k) - arg(1)*sqrt(2*arg(11)*h3(k)) + arg(2)*v2(k);
    h4(k+1) = h4(k) - arg(3)*sqrt(2*arg(11)*h4(k)) + arg(4)*v1(k);
    h1(k+1) = h1(k) - arg(5)*sqrt(2*arg(11)*h1(k)) + arg(6)*sqrt(2*arg(11)*h3(k)) + arg(7)*v1(k);
    h2(k+1) = h2(k) - arg(8)*sqrt(2*arg(11)*h2(k)) + arg(9)*sqrt(2*arg(11)*h4(k)) + arg(10)*v2(k);
end


%eq = sum(h1 - w(1:nr,1)).^2 + sum(h2 - w(1:nr,2)).^2;

%Variáveis de saída
% eq1 = (h1-w1(1:length(h1))).^2; eq2 = (h2-w2(1:length(h2))).^2;
% delta_u1 = zeros(nu,1); delta_u2 = zeros(nu,1);
% for k=1:(nu-1)
%     delta_u1(k+1) = u1(k+1)-u2(k);
%     delta_u2(k+1) = u2(k+1)-u2(k);
% end


