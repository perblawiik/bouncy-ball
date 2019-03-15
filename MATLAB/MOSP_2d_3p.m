clear;

% CONSTANTS
    % time between each approximation
h = 0.01;
    % spring constant
k = 50;
    % resistativitation constant
b = 5;
    % gravitational constant
g = 9.8;
    % Floor bounciness multiplier (0-1 preferably :3)
fBounce = 0;

% PARTICLES AND SPRINGS
% ***DEFINING m,X,I IS ENOUGH FOR THE CODE TO RUN, THE REST IS ADAPTIVE TO THIS***
    % masses [m]/ per particle
m = [1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1;];
    % particle x, y Pos [Xx Xy] / per particle
X = [20 10; 24 11; 26 13; 27 17; 26 21; 24 23; 20 24; 16 23; 14 21; 13 17; 14 13; 16 11; 20 17];
    % particle indices for spring bonds [i1 i2]/ per spring
I = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 9 10; 10 11; 11 12; 12 1;
     1 13; 2 13; 3 13; 4 13; 5 13; 6 13; 7 13; 8 13; 9 13; 10 13; 11 13; 12 13;]; 
 
BONDS = size(I,1);
POINTS = size(X,1);

NORM = zeros(size(X));
for p = 1:12
    in1 = X(p,:)-X(mod(p-2,12)+1,:);
    in2 = X(p,:)-X(mod(p,12)+1,:);
    normDir = in1/norm(in1)+in2/norm(in2);
    NORM(p,:) = normDir/norm(normDir);
end

% DEFINING S.S.VARIABLES, STARTING VALUES
    % starting velocity [Vx Vy]/ per particle
V = zeros(POINTS,2);
% ***V(index,:) = [a b]; for initial velocities***
    % V'
Vp = zeros(POINTS,2);
    % Fk spring starting force [F]/ per spring  (applied directionally later, depending on spring orientation)
Fk = zeros(BONDS,2);               
    % Fk'
Fkp = zeros(BONDS,2);             

% SIMULATION
    % how many animation frames
CYCLES = 5000;
    % where the particle animation positions are stored
animation = zeros(POINTS,2,CYCLES); animation(:,:,1) = X;
for cycle=1:CYCLES -1    
    
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
    V  = V  + h*Vp;
    Fk = Fk + h*Fkp;
    X  = X  + h*V;
    
    % Pin = k/Vin, Pin = Fin/A
    % Pluft = k/Vut, Pluft = Fluft/A
    % Fin = k/Vin/A
    % 
    % Fin - Fut = k/A*(Vin -Vut) 
    
    
    
    % Code that flips Y-ward velocity when the particle has Xy<0
    V(:,2) = (X(:,2)>0).*V(:,2)-fBounce*(X(:,2)<0).*V(:,2);
    % Sets Xy values to 0 if they're below 0
    X(:,2) = (X(:,2)>0).*X(:,2);
    
    animation(:,:,cycle + 1) = X;
end

% ANIMATION
% Stop the animation by pressing CTRL + C in the command window
for i = 1:CYCLES
    clf;
    hold on;
    plot([-50 50], [0 0]); 
    plot(animation(:,1,i),animation(:,2,i),'ro')
    for p = 1:12
        in1 = animation(p,:,i)-animation(mod(p-2,12)+1,:,i);
        in2 = animation(p,:,i)-animation(mod(p,12)+1,:,i);
        normDir = in1/norm(in1)+in2/norm(in2);
        if(abs(norm(normDir))>.001)
            NORM(p,:) = sign(dot(NORM(p,:),normDir))*normDir/norm(normDir);
        end
    end
    for n = 1:BONDS 
        plot([animation(I(n,1),1,i) animation(I(n,2),1,i)],[animation(I(n,1),2,i) animation(I(n,2),2,i)],'b--');
    end
    for n = 1:12 
        plot([animation(n,1,i) animation(n,1,i)+10*NORM(n,1)],[animation(n,2,i) animation(n,2,i)+10*NORM(n,2)],'g');
    end
    xlim([-15 40]);
    ylim([-5 40]);
    
    % used to set frame time but drawing takes time so it's not accurate
    pause(h); 
end





