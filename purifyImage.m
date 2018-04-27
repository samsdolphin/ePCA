function purified_img = purifyImage(img,gamma)

[m,n,~] = size(img);
Y = img(:,:,1);
Y_bar = zeros(m,1);

for i=1:n
    Y_bar = Y_bar + Y(:,i);
end
Y_bar = Y_bar/n;

S = zeros(m,m);
for i=1:n
    S = S + (Y(:,i)-Y_bar)*(Y(:,i)-Y_bar)';
end
S = S/n;
Sw = diag(sqrt(Y_bar))*S*diag(sqrt(Y_bar));

[V,D] = eig(Sw);
for i = 1:length(D)
    if D(i,i)>(1+sqrt(gamma))^2
        D(i,i) = (D(i,i)+1-gamma+sqrt((D(i,i)+1-gamma).^2-4.*D(i,i)))/2;
    else
        D(i,i) = 1;
    end
end
Sw = V*D*inv(V);