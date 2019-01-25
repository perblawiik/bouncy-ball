h = 0.01;

k = 20; b = 10; g = 0;

BONDS = 3;
POINTS = 4;
m = [50; 1; 1; 1];
X = [1 5; 1 15; 1 25; 11 25];
I = [1 2; 2 3; 3 4];
V = [0 0; 0 0; 0 0; 5 0];
Vp = [0 0; 0 0; 0 0; 0 0];
F = [0; 0; 0];
Fp = [0; 0; 0];

for cycle=1:2000
    
    Vp = zeros(POINTS,2);
    Fp = zeros(BONDS,1);
    for n = 1:BONDS
        dif = X(I(n,1),:) - X(I(n,2),:);
       nDif = dif/norm(dif);
       Fp(n) = k * (dot(V(I(n,1),:)-V(I(n,2),:),nDif) - 1/b*F(n));
       Vp(I(n,1),:) = Vp(I(n,1),:) - 1/m(I(n,1)) * F(n)*nDif;
       Vp(I(n,2),:) = Vp(I(n,2),:) + 1/m(I(n,2)) * F(n)*nDif;
    end
    Vp(3,:) = Vp(3,:) - [0 g];
    Vp
    V = V + h*Vp
    F = F + h*Fp
    X = X + h*V
    
%     V(:,2) = (X(:,2)>0).*V(:,2);
%     X(:,2) = (X(:,2)>0).*X(:,2);
    
    clf;
    hold on;
    plot([-50 50], [0 0]); 
    plot(X(:,1),X(:,2),'or');
    xlim([-15 40]);
    ylim([-5 40]);
    
    pause(h);
end




