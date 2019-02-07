clear;
close all

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
fBounce = 0.2; 

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
    
    % COLLISION ------------------------------------------
    % The considered ledges. [x1 y1; x2 y2]
    linePoints(:, :, 1) = [20 10; 50 40];
    linePoints(:, :, 2) = [10 -20; -20 10];
    
    % For each ledge
    for currentLedge = 1:2
        % CALCULATE THE LINE ------------------------------
        % y = k*x + m
        % Calculate "k"
        tilt= (linePoints(2, 2, currentLedge)-linePoints(1, 2, currentLedge))/(linePoints(2, 1, currentLedge)-linePoints(1, 1, currentLedge));
        
        % Calculate "m" (where the line crosses x = 0)
        cross = linePoints(1, 2, currentLedge)-tilt*linePoints(1, 1, currentLedge);
        
        % Set the line
        line = [-1*tilt 1 -1*cross];
        
        % Calculate the mid point of the ledge
        midPoint = 2*linePoints(:, 1, currentLedge)-linePoints(:, 2, currentLedge);
        
        normalRange = 0.5 *sqrt((linePoints(2, 2, currentLedge)-linePoints(1, 2, currentLedge))^2 + (linePoints(2, 1, currentLedge)-linePoints(1, 1, currentLedge))^2);
        dx = linePoints(1, 2, currentLedge)- linePoints(1, 1, currentLedge);
        dy = linePoints(2, 2, currentLedge)- linePoints(2, 1, currentLedge);
        normalV = [dy -dx];
        normalV= normalV/norm(normalV);
        %--------------------------------------------------
        % Each point's distance to the line
        distanceFromLine = (X(:,1)* line(1) + X(:,2)*line(2) + line(3))/(line(1).^2+line(2).^2 );
        
        % Distance from the normal
        distanceFromMid = sqrt(((X(:,1)-midPoint(1)).^2 + (X(:,2)-midPoint(2)).^2));
        
        % The indicies of the points below the line
        IndeciesBelowLine = zeros(1, size(X, 1));
        IndeciesBelowLine(distanceFromLine<0) = 1;
        
        % The indicies of the points close to the normal
        IndeciesCloseToNormal = zeros(1, size(X, 1));
        IndeciesCloseToNormal(abs(distanceFromMid) > normalRange) = 1;
        %--------------------------------------------------
        
        IndicesInContact = IndeciesCloseToNormal.*IndeciesBelowLine;
        TheIndices = find(IndicesInContact);
        
        if (length(TheIndices) > 0)
            % For each point below the ledge...
            for n = 1:length(TheIndices)
                % Stop velocity in normal direction
                V(TheIndices(n), :) = (V(TheIndices(n), :) - normalV*dot(V(TheIndices(n), :), normalV))*0.7; % normal * norm(V(IndeciesBelowLine(n), :) *fBounce);
            
                % Nudge all particles to the ledge
                X(TheIndices(n), :) =  X(TheIndices(n), :) - normalV * distanceFromLine(TheIndices(n));
            end
        end
    end
    animation(:,:,cycle + 1) = X;
end

% ANIMATION
% Stop the animation by pressing CTRL + C in the command window
for i = 1:CYCLES
    clf;
    hold on;
    plot([10 50], [0 40]);
    plot([-20 10], [10 -20]); 
    
    plot(animation(:,1,i),animation(:,2,i),'ro')
    for n = 1:BONDS 
        plot([animation(I(n,1),1,i) animation(I(n,2),1,i)],[animation(I(n,1),2,i) animation(I(n,2),2,i)],'b--');
    end
    xlim([-100 100]);
    ylim([-100 100]);
    
    % used to set frame time but drawing takes time so it's not accurate
    pause(h); 
end





