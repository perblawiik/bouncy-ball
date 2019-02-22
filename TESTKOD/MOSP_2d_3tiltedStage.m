clear;

% CONSTANTS
    % time between each approximation
h = 0.01;
    % spring constant
k = 200;
    % resistativitation constant
b = 5;
    % gravitational constant
g = 20;
    % Floor bounciness multiplier (0-1 preferably :3)
fBounce = 0.5; 

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
                                                                            % both pulled towards eachother or drawn away from eachother 
       Fkp(n) = k * dV;     % the derivative for Fk...
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
    plot(animation(:,1,i),animation(:,2,i),'ro')
    for n = 1:BONDS 
        plot([animation(I(n,1),1,i) animation(I(n,2),1,i)],[animation(I(n,1),2,i) animation(I(n,2),2,i)],'b--');
    end
    xlim([-50 50]);
    ylim([-50 50]);
    
    % used to set frame time but drawing takes time so it's not accurate
    pause(0); 
end





