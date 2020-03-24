function [obj_func] = objective_function(x,y)

% clear 
% close all

% theta = rad2deg(mod(atan2(y,x),2*pi)); % Between 0 and 2pi
theta = atan2(y,x); % Between -pi and +pi
mag = sqrt(x^2 + y^2);

% obj_func = ((x+y)/(x^2 + y^2)) * atan2(y,x);
obj_func = mag*theta;

% create figure if one doesn't exist
g = groot;
if isempty(g.Children)
    figure;
end

hold on
plot([0 x/mag],[0 y/mag],'LineWidth',2)
xlim([-1 1])
ylim([-1 1])
axis equal
grid on

end
