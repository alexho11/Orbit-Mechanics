function [ yp ] = yprime( ~,y )
GM=3.986005e14;
R=sqrt((y(1)^2)+(y(2)^2)+(y(3)^2));
acc=(-(GM/(R^3)).*y(1:3));
yp=[y(4);y(5);y(6);acc(1);acc(2);acc(3)];
end

