clear;

% CONSTANTS
    % how many seconds to simulate
time = 30;
    % time between each approximation
h = 0.01;
    % spring constant
k = 100;
    % resistativitation constant
b = 5;
    % gravitational constant
g = 9.8;
    % Floor bounciness multiplier (0-1 preferably :3)
fBounce = 0;
    % Floor friction (1 means full stopp, 0 means ice)
fFriction = 0;
    % Pressure constant
L = 10000;

% PARTICLES AND SPRINGS
% ***DEFINING m,X,I IS ENOUGH FOR THE CODE TO RUN, THE REST IS ADAPTIVE TO THIS***
    % mass M total for system
M = 10;
    % particle x, y Pos [Xx Xy] / per particle
X = [20 10; 24 11; 26 13; 27 17; 26 21; 24 23; 20 24; 16 23; 14 21; 13 17; 14 13; 16 11];
    % particle indices for spring bonds [i1 i2]/ per spring
I = [1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8; 8 9; 9 10; 10 11; 11 12; 12 1]; 

X = [0 40]+[10*cos((0:pi/20:2*pi-pi/20)') 10*sin((0:pi/20:2*pi-pi/20)')];
POINTS = size(X,1);

m = ones(POINTS)*M/POINTS;

% Bonds in a circle, 1-2, 2-3, 3-4, 4-1
for i = 1:POINTS-1
    I(i,:) = [i i+1];
end
I(POINTS,:) = [POINTS 1];
BONDS = size(I,1);


% all to all
% count = 1;
% for i = 1:POINTS
%     for j = i+1:POINTS
%         I(count,:) =[i j];
%         count = count +1;
%     end
% end
% BONDS = size(I,1);

% calculating normals
NORM = zeros(size(X));
for p = 1:POINTS
    in1 = X(p,:)-X(mod(p-2,POINTS)+1,:);
    in2 = X(p,:)-X(mod(p,POINTS)+1,:);
    normDir = in1/norm(in1)+in2/norm(in2);
    NORM(p,:) = normDir/norm(normDir);
end
NORMcop = NORM; % copy for animation

% calculate starting volume
V0 = 0;
mid = mean(X(:,:));
for p = 1:POINTS-1
    v1 = [X(p,:)-mid 0];
    v2 = [X(p+1,:)-mid 0];
    V0 = V0 + norm(cross(v1,v2))/2;
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
    
    
    % calculating Area->Volume and resulting pressure
    Area = 0;
    mid = mean(X(:,:));
    for p = 1:POINTS-1
        v1 = [X(p,:)-mid 0];
        v2 = [X(p+1,:)-mid 0];
        Area = Area + norm(cross(v1,v2))/2;
    end
    Vratio = Area/V0;
    Fpressure = L*(1 - Vratio)/POINTS;
    
    
    % calculate normals based on adjacent points and last normal
    for p = 1:POINTS
        in1 = X(p,:)-X(mod(p-2,POINTS)+1,:);
        in2 = X(p,:)-X(mod(p,POINTS)+1,:);
        normDir = in1/norm(in1)+in2/norm(in2);
        if(abs(norm(normDir))>.001)
            NORM(p,:) = sign(dot(NORM(p,:),normDir))*normDir/norm(normDir);
        end
    end
    
    % apply pressure
    for p = 1:POINTS
        v1 = [X(p,:)-mid];
        Vp(p,:) = Vp(p,:) + 1/m(p)*Fpressure*v1/norm(v1);
    end
    
    Vp = Vp - [0 g];    % gravity is added for all points
    
    % approximating the new values using: X_n+1 = X_n + h*X'_n
    V  = V  + h*Vp;
    Fk = Fk + h*Fkp;
    X  = X  + h*V;
    
%     % The considered ledges. [x1 y1; x2 y2]
%     linePoints(:, :, 1) = [10 0; 50 40]; 
%     linePoints(:, :, 2) = [-20 10; 10 -20]; 
%     linez = linePoints;
%     
%     % For each ledge
%     for currentLedge = 1:2
%         % CALCULATE THE LINE ------------------------------
%         % y = k*x + m
%         % Calculate "k"
%         tilt= (linePoints(2, 2, currentLedge)-linePoints(1, 2, currentLedge))/(linePoints(2, 1, currentLedge)-linePoints(1, 1, currentLedge));
%         
%         % Calculate "m" (where the line crosses x = 0)
%         crossover = linePoints(1, 2, currentLedge)-tilt*linePoints(1, 1, currentLedge);
%         
%         % Set the line
%         theLine = [-1*tilt 1 -1*crossover];
%         
%         % Calculate the mid point of the ledge
%         midPointY = linePoints(2, 2, currentLedge)-linePoints(1, 2, currentLedge);
%         
%         % Calculate the normal-line of the line
%         normalLine = [tilt 1 midPointY-crossover];
%         % -------------------------------------------------
%         % Normal _vector_
%         normalV = [theLine(1) theLine(2)];
%         normalV = normalV/norm(normalV);
%         %--------------------------------------------------
%         % Each points' distance to the line
%         distanceFromLine = (X(:,1)* theLine(1) + X(:,2)*theLine(2) + theLine(3))/(theLine(1).^2+theLine(2).^2 );
%         
%         % The indicies of the points below the line
%         IndeciesBelowLine = find(distanceFromLine<0);
%         
%         % Distance from the normal
%         distanceFromNormal = (X(:,1)* normalLine(1) + X(:,2)*normalLine(2) + normalLine(3))/(normalLine(1).^2+normalLine(2).^2 );
%         
%         % The indicies of the points below the line
%         IndeciesCloseToNormal = find(abs(distanceFromNormal)<sqrt((linePoints(2, 2, currentLedge)-linePoints(1, 2, currentLedge))^2 + (linePoints(2, 1, currentLedge)-linePoints(1, 1, currentLedge))^2) );
%         %--------------------------------------------------
%         % Initialize the amount of points in contact with ledge
%         countIndiciesInContact= 0;
%         
%         if ((length(IndeciesBelowLine)*length(IndeciesCloseToNormal))~= 0)
%             % Create a temporaty list of all indicies that are in contact with
%             % the ledge
%             tempIndiciesInContact = zeros(length(IndeciesBelowLine));
% 
%             % For each point below the ledge..
%             for belowLine = 1:length(IndeciesBelowLine)
%                 % For each point far away from the normal..
%                 for closeToNormal = 1:length(IndeciesCloseToNormal)
%                     % Chech if the point below the ledge is also far away from
%                     % the normal
%                     if IndeciesBelowLine(belowLine) == IndeciesCloseToNormal(closeToNormal)
%                         % If it is, add it to the temporary list
%                         tempIndiciesInContact( countIndiciesInContact + 1) = IndeciesBelowLine(belowLine);
%                         % One additional point in contact with ledge
%                         countIndiciesInContact= countIndiciesInContact + 1;
%                         % No need to continue search
%                         break;
%                     end
%                 end  
%             end
%             % 
%             IndiciesInContact = zeros(1, countIndiciesInContact);
%             for inContact = 1:countIndiciesInContact
%                 IndiciesInContact(inContact) = tempIndiciesInContact(inContact);
%             end
%             IndiciesInContact= IndiciesInContact';
%         end
%         % Check if there are any points below the line
%         if (countIndiciesInContact ~= 0)
%             % For each point below the ledge...
%             for n = 1:length(IndiciesInContact)
%                 % flip velocities of those below
%                 V(IndiciesInContact(n), :) = (V(IndiciesInContact(n), :) - normalV*dot(V(IndiciesInContact(n), :), normalV))*1.0; % normal * norm(V(IndeciesBelowLine(n), :) *fBounce);
%             
%                 % Nudge all particles below the line to the line
%                 X(IndiciesInContact(n), :) =  X(IndiciesInContact(n), :) - normalV * distanceFromLine(IndiciesInContact(n));
%             end
%         end
%     end
    
    % Code that flips Y-ward velocity when the particle has Xy<0
    V(:,2) = (X(:,2)>0).*V(:,2)-fBounce*(X(:,2)<0).*V(:,2);
    % Sets Xy values to 0 if they're below 0
    X(:,2) = (X(:,2)>0).*X(:,2);
    % Code that slows X-speed on floor collision
    V(:,1) = (X(:,2)>0).*V(:,1)+(1-fFriction)*(X(:,2)<=0).*V(:,1);
    
    animation(:,:,cycle + 1) = X;
end


scale = 100;

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






%%
NORM = NORMcop;
% ANIMATION
% Stop the animation by pressing CTRL + C in the command window
for i = 1:CYCLES
    % CALCULATIONS FOR EXTRA STUFF
    for p = 1:POINTS
        in1 = animation(p,:,i)-animation(mod(p-2,POINTS)+1,:,i);
        in2 = animation(p,:,i)-animation(mod(p,POINTS)+1,:,i);
        normDir = in1/norm(in1)+in2/norm(in2);
        if(abs(norm(normDir))>.001)
            NORM(p,:) = sign(dot(NORM(p,:),normDir))*normDir/norm(normDir);
        end
    end
    Area = 0;
    mid = mean(animation(:,:,i));
    for p = 1:POINTS-1
        v1 = [animation(p,:,i)-mid 0];
        v2 = [animation(p+1,:,i)-mid 0];
        Area = Area + norm(cross(v1,v2))/2;
    end

%     Vratio = Area/V0;
%     Fpressure = L*(1 - Vratio)/POINTS;
    
    
    clf;
    hold on;
    plot([-50 50], [0 0]); 

    % ONLY CONTOUR
    for n = 1:BONDS 
        plot([animation(I(n,1),1,i) animation(I(n,2),1,i)],[animation(I(n,1),2,i) animation(I(n,2),2,i)],'b-');
    end

    %EXTRA STUFF
    plot(animation(:,1,i),animation(:,2,i),'ro')
    plot(mid(1),mid(2),'og')
    for n = 1:POINTS 
        plot([animation(n,1,i) animation(n,1,i)+10*NORM(n,1)],[animation(n,2,i) animation(n,2,i)+10*NORM(n,2)],'g');
    end
    
    scale = 50;
    xlim([-scale*0.5 scale*0.5]);
    ylim([-scale*0.2 scale*0.8]);
    
    % used to set frame time but drawing takes time so it's not accurate
    pause(0);   
end




