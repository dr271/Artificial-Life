function [z] = RippleLandscape(x)
	z=(4*cos((x(1)^2+x(2)^2)/4)/(x(1)^2+x(2)^2+1))*2;
end