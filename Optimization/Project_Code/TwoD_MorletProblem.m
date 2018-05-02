function val = TwoD_MorletProblem(x) 
% Returns the value of a two-dimensional Morlet wavelet, defined in HW#3
val = exp(-1*(x(1)).^2).*cos(pi*(x(1)))-pi/2.*exp(-2*(x(2)).^2).*sin(2*pi*x(2)) ;
end