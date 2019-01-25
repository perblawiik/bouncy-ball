
% Stop the simulation by pressing CTRL + C in the command window



clear;

h = 0.01;   % this time is very important depending on later values, the smaller the better

k = 100;    % spring constant
b = 5;      % resistativitation constant
g = 9.8;    % gravitational constant


m = [1; 1; 1; 1];                   % masses [m]/ per particle
X = [10 10; 20 10; 15 15; 25 15];   % particle x, y Pos [Xx Xy] / per particle
I = [1 2; 2 3; 1 3; 3 4; 2 4];      % particle indices for spring bonds [i1 i2]/ per spring
V = [0 0; 0 0; 0 0; 0 0];           % starting velocity [Vx Vy]/ per particle
Vp = [0 0; 0 0; 0 0; 0 0];          % V'
Fk = [0; 0; 0; 0; 0];               % Fk spring starting force [F]/ per spring  (applied directionally later, depending on spring orientation)
Fkp = [0; 0; 0; 0; 0];              % Fk'

BONDS = size(I,1);
POINTS = size(X,1);

for cycle=1:5000    % how many steps will be run
    
    Vp = zeros(POINTS,2);       % set to zero so the components from each connected spring can be += and added separately
    for n = 1:BONDS     % Loop through the springs
       dif = X(I(n,1),:) - X(I(n,2),:);     % Gets vector from particle 1 to 2
       nDif = dif/norm(dif);                % normalises it, used to give the Fk and Fb direction
       dV = dot(V(I(n,1),:)-V(I(n,2),:),nDif);      % Gets deltaV, speed difference between the particles in the spring's direction
       Vp(I(n,1),:) = Vp(I(n,1),:) - 1/m(I(n,1)) * (b*dV + Fk(n))*nDif;     % apply spring influence to the connected particles
       Vp(I(n,2),:) = Vp(I(n,2),:) + 1/m(I(n,2)) * (b*dV + Fk(n))*nDif;     % first one's added, second is subtracted, in the springs direction since they will be either
                                                                            % both pulled towards eachother or drawn away from eachother 
       Fkp(n) = k * dV;     % the derivative for Fk...
    end
    Vp = Vp - [0 g];    % gravity is added for all points
    
    % approximating the new values using: X_n+1 = X_n + h*X'_n
    % they're not supressed for debugging purposes
    V  = V  + h*Vp
    Fk = Fk + h*Fkp
    X  = X  + h*V
    
    % Code that flips Y-ward velocity when the particle has Xy<0
    V(:,2) = (X(:,2)>0).*V(:,2)-(X(:,2)<0).*V(:,2);
    
    
    % draw each frame in "real timeish" (atm)
    clf;
    hold on;
    plot([-50 50], [0 0]); 
    plot(X(:,1),X(:,2),'or');
    for n = 1:BONDS 
        plot([X(I(n,1),1) X(I(n,2),1)],[X(I(n,1),2) X(I(n,2),2)],'b');
    end
    xlim([-15 40]);
    ylim([-5 40]);
    
    pause(h); % bad commmand since computations take time making each frame longer than h seconds
end




