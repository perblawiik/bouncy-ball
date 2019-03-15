
clear;
close all;

% CONSTANTS
    % how many seconds to simulate
time = 30;
    % time between each approximation
h = 0.01;
    % spring constant
k = 1000;
    % resistativitation constant
b = 10;
    % gravitational constant
g = 9.8;
    % Floor bounciness multiplier (0-1 preferably :3)
fBounce = 0;
    % Floor friction (1 means full stopp, 0 means ice)
fFriction = 0;
    % Springs cant contract closer than this
minSpring = 0.8;


% PARTICLES AND SPRINGS
% ***DEFINING m,X,I IS ENOUGH FOR THE CODE TO RUN, THE REST IS ADAPTIVE TO THIS***
    % mass M total for system
M = 10;
    % particle x, y Pos [Xx Xy] / per particle
X = [0 20]+[5*cos((0:pi/5:2*pi-pi/5)') 5*sin((0:pi/5:2*pi-pi/5)')];
POINTS = size(X,1);

m = ones(POINTS)*M/POINTS;

%all to all
count = 1;
for i = 1:POINTS-1
    for j = i+1:POINTS
        I(count,:) = [i j];
        count = count + 1;
    end
end
BONDS = size(I,1);

% DEFINING S.S.VARIABLES, STARTING VALUES
    % starting velocity [Vx Vy]/ per particle
V = zeros(POINTS,2);
% ***V(index,:) = [a b]; for initial velocities***
    % V'
Vp = zeros(POINTS,2);
    % Fk spring starting force [F]/ per spring  (applied directionally later, depending on spring orientation)
Fk = zeros(BONDS,1);               
    % Fk'
Fkp = zeros(BONDS,1);             

% SIMULATION
    % how many animation frames
CYCLES = time/h;
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
    V  = V  + h * Vp;
    Fk = Fk + h*Fkp;
    X  = X  + h*V;
    
    % Code that flips Y-ward velocity when the particle has Xy<0
    V(:,2) = (X(:,2)>0).*V(:,2)-fBounce*(X(:,2)<0).*V(:,2);
    % Sets Xy values to 0 if they're below 0
    X(:,2) = (X(:,2)>0).*X(:,2);
    % Code that slows X-speed on floor collision
    V(:,1) = (X(:,2)>0).*V(:,1)+(1-fFriction)*(X(:,2)<=0).*V(:,1);
    
    animation(:,:,cycle + 1) = X;
end


scale = 20;

% ANIMATION
% Stop the animation by pressing CTRL + C in the command window
for i = 1:CYCLES
    clf;
    hold on;
%     for l = 1:size(linez,3)
%         plot([linez(1,1,l) linez(2,1,l)], [linez(1,2,l) linez(2,2,l)]); 
%     end
    line([-50 50], [0 0]); 


    % ONLY CONTOUR
    line([animation(:,1,i);animation(1,1,i)],[animation(:,2,i);animation(1,2,i)],'Color','black');

    xlim([-scale*0.5 scale*0.5]);
    ylim([-scale*0.2 scale*0.8]);
    
    % used to set frame time but drawing takes time so it's not accurate
    pause(0);   
end
