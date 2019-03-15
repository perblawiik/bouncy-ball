
clear;
close all;

% CONSTANTS
    % how many seconds to simulate
time = 10;
    % time between each approximation
h = 0.01;
    % spring constant
k = 0.2;
    % resistativitation constant
b = 0.01;
    % gravitational constant
g = 9.8;
    % Floor bounciness multiplier (0-1 preferably :3)
fBounce = 0;
    % Floor friction (1 means full stopp, 0 means ice)
fFriction = 0;
    % Pressure constant
L = 1000;


% PARTICLES AND SPRINGS
% ***DEFINING m,X,I IS ENOUGH FOR THE CODE TO RUN, THE REST IS ADAPTIVE TO THIS***
    % mass M total for system
M = 10;
    % particle x, y Pos [Xx Xy] / per particle
X = [0 100]+[10*cos((0:pi/20:2*pi-pi/20)') 10*sin((0:pi/20:2*pi-pi/20)')];
POINTS = size(X,1);

m = ones(POINTS)*M/POINTS;

    % particle indices for spring bonds [i1 i2]/ per spring
% Bonds in a circle, 1-2, 2-3, 3-4, 4-1
count = 1;
for i = 1:POINTS
    for j = i+1:POINTS
        I(count,:) = [i j];
        count = count + 1;
    end
end
I(POINTS,:) = [POINTS 1];
BONDS = size(I,1);

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
    
    
    % calculating Area->Volume and resulting pressure
    circum = 0;
    Area = 0;
    mid = mean(X(:,:));
    for p = 1:POINTS-1
        v1 = [X(p,:)-mid 0];
        v2 = [X(p+1,:)-mid 0];
        Area = Area + norm(cross(v1,v2))/2;
        circum = circum + norm(X(p,:) - X(p+1,:));
    end
    
    Vratio = Area/V0;
    Fpressure = L*(1 - Vratio)*(circum/POINTS);
    
    
    % calculate normals based on adjacent points and last normal
    for p = 1:POINTS
        in1 = X(p,:)-X(mod(p-2,POINTS)+1,:);
        in2 = X(p,:)-X(mod(p,POINTS)+1,:);
        normDir = in1/norm(in1)+in2/norm(in2);
        if(abs(norm(normDir))>.001)
            NORM(p,:) = sign(dot(NORM(p,:),normDir))*normDir/norm(normDir);
        end
    end
    
    % apply pressure force
    for p = 1:POINTS
        v1 = [X(p,:)-mid];
        Vp(p,:) = Vp(p,:) + 1/m(p)*Fpressure*v1/norm(v1);
        %Vp(p,:) = Vp(p,:) + 1/m(p)*Fpressure*NORM(p,:);
    end
    
    Vp = Vp - [0 g];    % gravity is added for all points
    
    % approximating the new values using: X_n+1 = X_n + h*X'_n
    V  = V  + h*Vp;
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
    circum = 0;
    Area = 0;
    mid = mean(animation(:,:,i));
    for p = 1:POINTS-1
        v1 = [animation(p,:,i)-mid 0];
        v2 = [animation(p+1,:,i)-mid 0];
        Area = Area + norm(cross(v1,v2))/2;
        circum = circum + norm(animation(p,:,i) - animation(p+1,:,i));
    end
    Vratio = Area/V0;
    Fpressure = L*(1 - Vratio)*(circum/POINTS);

    
    
    
    clf;
    hold on;
    plot([-50 50], [0 0]); 

    % ONLY CONTOUR
    plot([animation(1:POINTS,1,i)' animation(1,1,i)],[animation(1:POINTS,2,i)' animation(1,2,i)],'b-');

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




