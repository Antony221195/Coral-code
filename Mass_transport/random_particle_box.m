function [xpos,ypos,zpos] = random_particle_box(num_particle,xrow,yrow,H)

% Define the box limits in um
xmin = xrow(2); xmax = xrow(end-1);
ymin = yrow(2); ymax = yrow(end-1);
zmin = 0; zmax = H;

% Generate random coordinates
xpos = xmin + (xmax - xmin) * rand(num_particle, 1);
ypos = ymin + (ymax - ymin) * rand(num_particle, 1);
zpos = zmin + (zmax - zmin) * rand(num_particle, 1);

end