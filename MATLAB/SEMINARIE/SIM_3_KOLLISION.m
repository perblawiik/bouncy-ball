clear;
close all;

% CONSTANTS
    % time between each approximation
h = 0.01;
    % spring constant
k = 20;
    % resistativitation constant
b = 1;
    % gravitational constant
g = 9.8;
    % Floor bounciness multiplier (0-1 preferably)
fBounce = 0;

% PARTICLES AND SPRINGS
% ***DEFINING m,X,I IS ENOUGH FOR THE CODE TO RUN, THE REST IS ADAPTIVE TO THIS***
    % particle x, y Pos [Xx Xy] / per particle
r = 8;
seg = 10;
center = [20 30];
X = center + [r*cos(0:2*pi/seg:2*pi*(1-1/seg))' r*sin(0:2*pi/seg:2*pi*(1-1/seg))'];
X = cat(1, X, center);
POINTS = size(X,1);

    % masses [m]/ per particle
m_tot = 1;
m = ones(POINTS,1)*m_tot/POINTS;

    % particle indices for spring bonds [i1 i2]/ per spring
I = zeros(2*seg,2);
% Bonds along the circumference
for i = 1:seg-1
   I(i,:) = [i i+1];
end
I(seg,:) = [1 seg];
% Bonds to the center point
for i = 1:seg
   I(seg + i,:) = [i seg+1];
end
BONDS = size(I,1);

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
CYCLES = 2000;
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
       Fkp(n) = k * dV;
                                                                         % both pulled towards eachother or drawn away from eachother 
    end
    Vp = Vp - [0 g];    % gravity is added for all points
    
    % approximating the new values using: X_n+1 = X_n + h*X'_n
    V  = V  + h*Vp;
    Fk = Fk + h*Fkp;
    Xprevious = X;
    X  = X  + h*V;
    
    % COLLISION ------------------------------------------
    % The considered ledges. [x1 y1; x2 y2]
    linePoints(:, :, 1) = [10 5; 40 10];
    linePoints(:, :, 2) = [0 -40; 20 -15];
    linePoints(:, :, 3) = [-20 15; -5 -10];
    linePoints(:, :, 4) = [10 0; 10 5];
    linePoints(:, :, 5) = [20 -15; 10 0];
    linePoints(:, :, 6) = [-5 -10; -15 -20];
    linePoints(:, :, 7) = [-15 -20; -15 -30];
    linePoints(:, :, 8) = [-15 -30; -40 -30];
    linePoints(:, :, 9) = [-40 -50; 0 -40];
    
%     linePoints(:, :, 1) = [11 -20; 21 0];
%     linePoints(:, :, 2) = [21 0; 31 -20];
    
    
    for currentLedge = 1:size(linePoints, 3)
        for currentPoint = 1:POINTS
            L0 = Xprevious(currentPoint, : )';
            L1 = X(currentPoint,: )';
            P0 = linePoints(1, :, currentLedge)';
            P1 = linePoints(2, :, currentLedge)';

            dPx = P1(1)-P0(1);
            dPy = P1(2)-P0(2);

            Pnormal = [dPy -1*dPx];
            Pnormal = Pnormal/norm(Pnormal);

            L01 = L1 - L0;
            P01 = P1 - P0;

            if det([P01 -1*L01])~=0
                vt = [P01 -L01]\(L0-P0);
                v = vt(1);
                t = vt(2);

                if ( (0<=t)&&(t<=1)&&(0<=v)&&(v<=1) )
                    L = L0 + L01*t;
                    L = L';
                    X(currentPoint,: ) = L - 0.01*Pnormal;
                    V(currentPoint,: ) = (V(currentPoint,: )-Pnormal*dot(V(currentPoint,: ), Pnormal))*fBounce;
                end
            end

        end % End of the point
    end % End the ledge
    animation(:,:,cycle + 1) = X;
end

% ANIMATION
% Stop the animation by pressing CTRL + C in the command window
figure(1)

for i = 1:3:CYCLES
    clf;
    hold on;
    for ledge = 1:size(linePoints, 3)  
        plot([linePoints(1, 1, ledge) linePoints(2, 1, ledge)], [linePoints(1, 2, ledge) linePoints(2, 2, ledge)]);
    end
    
%     plot([10 50], [0 40]);
%     plot([-20 10], [10 -20]); 
%     
    for n = 1:BONDS 
        line([animation(I(n,1),1,i) animation(I(n,2),1,i)],[animation(I(n,1),2,i) animation(I(n,2),2,i)]);
    end
    xlim([-50 50]);
    ylim([-50 50]);
    
    % used to set frame time but drawing takes time so it's not accurate
    pause(0); 
end





