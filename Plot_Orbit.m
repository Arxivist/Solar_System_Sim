function [x0,y0,z0,Vel_ux, Vel_uy, Vel_uz] = Plot_Orbit (a,e,i,omega,w,M0,num_color)

%{
clear all;
close all;
clc;

model = figure('Name','Orbit');
hold on;
grid on;
set(model, 'Position', get(0, 'Screensize'));
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Orbital Elements

0 - a (Semimajor axis, Au)
1 - e (eccentricity)
2 - i (incluniation, Deg)
3 - omega (longitude of ascending node, Deg)
4 - w (argument of periapsis, Deg)
5 - M0 (true anomoly at epoch, longitude Jan 2000, deg )
6 - Mass

h = angular momentum (Km/s)
r = radius from centre


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}



%% Define Paramaters
Mercury(1:7) = [0.38709893,0.20563069,7.00487,48.33167,77.45645,252.25084,3.285e23];
  % test cases
  %{
  model = figure('Name','Orbit');
  hold on;
  grid on;
  set(model, 'Position', get(0, 'Screensize'));
  axis([-2,2,-2,2,-2,2])
  a = 2;
  e = 0.5;
  i = 5;
  omega = 0;
  w = 0;
  num_color = 3;
  M0 = 135;
  %}
  
  
  %convert to radians
  i = i*(pi/180);
  w = w*(pi/180);
  omega = omega*(pi/180);
  M0 = M0*(pi/180);
  
  % convert orb_color to string
  if (num_color == 1)
    orb_color = 'm';
  elseif (num_color == 2)
    orb_color = 'g';
  elseif (num_color == 3)
    orb_color = 'b';
  elseif (num_color == 4)
    orb_color = 'r';
  else
    orb_color = 'k';
  endif
  
  %Centered on 0
  Cx = 0;
  Cy = 0;
  Cz = 0;
  
  t = linspace(0,2*pi,200);
  
  %finding semiminor axis (b) from semimajor axis (a) and eccentricity (e)
  % e = sqrt(1-(b^2/a^2)); 
  b = a*sqrt(1-e^2);
  
  % find linear eccentricity (c)
  % distance between two focal paramaters
  c = a*e;
  
  
  % find semi Latus Rectum - line going through one of the foci
  
  
  
  %find x,y,z components of semi-major and semi-minor vectors (U for semi major, and V for semi minor)
  
  Vz = a*sin(i);
  Vxy = a*cos(i);
  Vy = Vxy * sin(w);
  Vx = Vxy * cos(w);
  %plot semi-major axis
  plot3([0 Vx], [0 Vy],[0 Vz], 'Color', orb_color);
  % turn into unit vector
  Vr = sqrt(Vx^2 + Vy^2+Vz^2);
  Vz = Vz/Vr; Vx = Vx/Vr; Vy = Vy/Vr;
  
  %semi-minor axis lines up with ascending node, no z component
  Uz = 0;
  Uy = b*sin(w-(pi/2));
  Ux = b*cos(w-(pi/2));
  
  %plot semi-minor axis
  plot3([0 Ux], [0 Uy],[0 Uz], 'Color', orb_color,'linestyle','-.');
  % turn into unit vector
  Ur = sqrt(Ux^2 + Uy^2+Uz^2);
  Uz = Uz/Ur; Uy = Uy/Ur; Ux = Ux/Ur;
  
  
  
  % Calculate entire orbit
  x = Cx + a.*cos(t) * Vx + b.*sin(t)*Ux;
  y = Cy + a.*cos(t) * Vy + b.*sin(t)*Uy;
  z = Cz + a.*cos(t) * Vz + b.*sin(t)*Uz;
  
   
  % Calculate starting position of orbital body
  x0 = Cx + a.*cos(M0) * Vx + b.*sin(M0)*Ux;
  y0 = Cy + a.*cos(M0) * Vy + b.*sin(M0)*Uy;
  z0 = Cz + a.*cos(M0) * Vz + b.*sin(M0)*Uz;
  %scatter3(x0,y0,z0,8,'k');
  

  % calculate starting velocity unit vector for the orbital body     
  Vel_x = -a.*sin(M0) * Vx + b.*cos(M0)*Ux;
  Vel_y = -a.*sin(M0) * Vy + b.*cos(M0)*Uy;
  Vel_z = -a.*sin(M0) * Vz + b.*cos(M0)*Uz;
  Vel = sqrt(Vel_x^2+Vel_y^2+Vel_z^2);
  Vel_ux = Vel_x/Vel; Vel_uy = Vel_y/Vel; Vel_uz = Vel_z/Vel; 
  
  %figure;
  plot3(x,y,z,'Color',orb_color);
  xlabel('X Axis');
  ylabel('Y Axis');
  zlabel('Z Axis');
  
  hold on;
 
  %{

% 2d ellipse
t = linspace(0,2*pi,100);
  % ellipse formula
  % [(x-x0)/a]^2 + [(y-y0)/b]^2 = 1
  % x(t) = a*sin(t)
  % y(t) = b*cos(t)
  %%Plot ellipse. Easiest to plot  generalized, then transform it.
  x = a.*sin(t);
  y = b.*cos(t);
  x0 = (x-Cx).*cos(w) + (y-Cy).*sin(w);
  y0 = -1*(x-Cx).*sin(w) + (y-Cy).*cos(w);
  plot(x0,y0, 'color','r');
  axis([-1,1,-1,1])
%}
%{ 
Handy reference equations
  [x, y, z] = ellipsoid(0,0,0,2.0,4.0,0,20);
  surf(x, y,z)
  colormap copper
  axis equal
  
  majorAxis = 2;
  minorAxis = 1;
  centerX = 10;
  centerY = 15;
  orientation = -45;
  theta = linspace(0,2*pi,150);
  orientation=orientation*pi/180;
  xx = (majorAxis/2) * sin(theta) + centerX;
  yy = (minorAxis/2) * cos(theta) + centerY;

  xx2 = (xx-centerX)*cos(orientation) - (yy-centerY)*sin(orientation) + centerX;
  yy2 = (xx-centerX)*sin(orientation) + (yy-centerY)*cos(orientation) + centerY;
  plot(xx2,yy2)
  axis equal
  grid
  
  
  ellipse in 3 space
  x(t) = Cx + a cos(t) Ux + b sin(t) Vx
  y(t) = Cy + a cos(t) Uy + b sin(t) Vy
  z(t) = Cz + a cos(t) Uz + b sin(t) Vz
  %}
