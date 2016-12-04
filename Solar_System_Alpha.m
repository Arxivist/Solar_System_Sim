clear all;
close all;

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
General Info

0 - a (Semimajor axis, Au)
1 - e (eccentricity)
2 - i (incluniation, Deg)
3 - omega (longitude of ascending node, Deg)
4 - w (argument of periapsis)
5 - M0 (true anomoly at epoch, longitude Jan 2000, deg )
6 - Mass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PARAMATERS

%% Initialize figure
 model = figure('Name','Orbits');
  hold on;
  grid on;
  %set(model, 'Position', get(0, 'Screensize'));
  axis equal;
  
  MSun = 1.989e30;
  
  %% Number of planets
  N = 4;

  %%  initial step distance
  DELTA = 10000;

  %% Gravitational Constant
  GAMMA = 6.67408e-11;

%mercuy inclination

  planets = zeros (N,7);
    %% Initialize Mercury
    planets(1,1:7) = [0.38709893*1.496e+11,0.20563069,7.00487,48.33167,77.45645,252.25084,3.285e23];

    %% Initialize Venus
    planets(2,1:7) = [0.72333199*1.496e+11,0.00677323,3.39471,76.68069,131.53298,181.97973,4.87e24];

    %% Initialize Earth
    planets(3,1:7) = [1.00000011*1.496e+11,0.01671022,0.00005,-11.26064,102.94719,100.46435, 5.972e24];

    %% Initialize Mars
    planets(4,1:7) = [1.52366231*1.496e+11,0.09341233,1.85061,49.57854,336.04084,355.45332,6.39e23];

    %% initialize the mass of the sun
    planets(N+1,7) = 1.989e30;
  %initial positions and velocity
  for i = 1:N
    [Pos0(i,1),Pos0(i,2),Pos0(i,3),Ux,Uy,Uz] = Plot_Orbit(planets(i,1),planets(i,2),planets(i,3),planets(i,4),planets(i,5),planets(i,6),i);
    % position comes out of PLot_Orbit
    
    %  Velocity is harder. first need to calculate precise orbital speed: v = sqrt(G*M*(2/r-1/a))
    % first calculate r
    r = sqrt(Pos0(i,1)^2+Pos0(i,2)^2+Pos0(i,3)^2);
    % then calculate v
    v = sqrt(GAMMA*planets(N+1,7)*(2/r-1/planets(i,1)));
    % then turn v into a vector
    Vel0(i,1) = v*Ux;Vel0(i,2) = v*Uy;Vel0(i,3) = v*Uz;
    %plot3([Pos0(i,1) Ux],[Pos0(i,2) Uy],[Pos0(i,3) Uz],8,'k');
  endfor
  %initial position and velocity of the Sun
  Pos0(N+1,3) = 0;
  Vel0(N+1,3) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Orbital Calculations

% see what would happen to the solar system if the sun changed mass
%planets(N+1,7) = 3e30;
  
%plot initial points
%hPlot = plot3(Pos0(1,1),Pos0(1,2),Pos0(1,3),'r*');

for i = 1:N
  hPlot(i) = plot3(Pos0(i,1),Pos0(i,2),Pos0(i,3),'ko');
  %hScatter(i) = scatter3(Pos0(i,1),Pos0(i,2),Pos0(i,3),8,'k');
  %hScatter(i) = scatter3('XData',Pos0(i,1),'YData',Pos0(i,2),'ZData',Pos0(i,3),'MarkerSize',8,'Color','k');
endfor


for LOOP_NUM = 1:1000
  %% determine Gravity
  G = zeros(N+1,3);
  A = zeros(N+1,3);
  for i = 1:N+1
    for j = 1:N+1
      if (j != i)          
         for x = 1:3
          G(i,x) = G(i,x) + GAMMA*((planets(i,7)*planets(j,7)*(Pos0(j,x)-Pos0(i,x)))/((pdist(Pos0(j,1:3),Pos0(i,1:3)))^3));
          %% determine planetary acceleration from gravitational force (F = MA, A =F/M)
          A(i,x) = G(i,x)/planets(i,7);
         endfor 
      endif
    endfor
  endfor
  


  %% sympletic Euler methods
  for i = 1:N
    for x = 1:3;
      Vel0(i,x) = Vel0(i,x) + (DELTA.*A(i,x));
      Pos0(i,x) = Pos0(i,x) + (DELTA.*Vel0(i,x));
    endfor
  endfor

  
  %% Draw the points


  for i = 1:N
    %hScatter(i) = scatter3(Pos0(i,1),Pos0(i,2),Pos0(i,3));
    set (hPlot(i),'xdata',Pos0(i,1),'ydata',Pos0(i,2),'zdata',Pos0(i,3));
  endfor
  drawnow
  
endfor
   
