function [z] = CamelLandscape(x)
    term1 = (4-2.1*x(1)^2+(x(1)^4)/3) * x(1)^2;
    term2 = x(1)*x(2);
    term3 = (-4+4*x(2)^2) * x(2)^2;

    z = term1 + term2 + term3;
end