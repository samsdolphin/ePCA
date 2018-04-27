clear all
r = 10;
p = 1000;
gamma = 4;
n = p/gamma;
A = 25*(1+sqrt(gamma))*(1+sqrt(gamma));

a = rand(n,r);
v = rand(p,r);
for i=1:n
    norm_ai = sum(a(i,:));
    for j=1:r
        a(i,j) = a(i,j)/norm_ai*A;
    end
end
X = v*a'; % Xi = ai1*v1 + ... + air*vr

%Y = X + sqrt(diag(X))*kesi;
Y = poissrnd(X);
Y_bar = zeros(p,1);
for i=1:n
    Y_bar = Y_bar + Y(:,i);
end
Y_bar = Y_bar/n;

S = zeros(p,p);
for i=1:n
    S = S + (Y(:,i)-Y_bar)*(Y(:,i)-Y_bar)';
end
S = S/n;
Sw = diag(sqrt(Y_bar))\S*inv(diag(sqrt(Y_bar))) - eye(p);
eigenv = eig(Sw);

bp = (1+sqrt(gamma))*(1+sqrt(gamma));
bn = (1-sqrt(gamma))*(1-sqrt(gamma));
t = bn-1:0.2:bp-1;
histogram(real(eigenv),0:0.2:8,'Normalization','pdf');
hold on;
g = fplot(@(t) sqrt((bp-(t+1)).*((t+1)-bn))/(2*pi*gamma*(t+1))); % multiply gamma if gamma>1
g.Color = 'r';
g.LineWidth = 2;
g.XRange = [0 10];
hold off;