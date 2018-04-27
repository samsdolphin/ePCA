function purified_img = purifyImage(set)

[h,w,n] = size(set);
p = h*w;
gamma = p/n;
X = zeros(p,n);
Y = X;
for i = 1:n
    Y(:,i) = reshape(double(set(:,:,i)),[p,1]);
end

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
%Sd = S - diag(Y_bar);
Sw = diag(sqrt(Y_bar))\S*inv(diag(sqrt(Y_bar)));

[V,D] = eig(Sw);
D = real(D);
for i = 1:length(D)
    if D(i,i)>(1+sqrt(gamma))^2
        D(i,i) = (D(i,i)+1-gamma+sqrt((D(i,i)+1-gamma).^2-4.*D(i,i)))/2;
    else
        D(i,i) = 1;
    end
end
Sw = V*D*inv(V);
Sw = diag(sqrt(Y_bar))*Sw*diag(sqrt(Y_bar));

for i = 1:n
    X(:,i) = Sw*(diag(Y_bar)+Sw)\Y(:,i) + diag(Y_bar)*(diag(Y_bar)+Sw)\Y_bar;
end

purified_img = X;

end
