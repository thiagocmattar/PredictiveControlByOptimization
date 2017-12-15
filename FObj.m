function eq = J(u,w)

nu = 80;
nr = 120; 

[h1_hat,h2_hat] = system_response(u(:,1),u(:,2),nu,nr);

eq = sum(h1_hat - w(1:length(h1_hat),1)).^2 + sum(h2_hat - w(1:length(h2_hat),2)).^2;





