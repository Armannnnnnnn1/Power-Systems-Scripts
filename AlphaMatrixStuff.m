clc; clear all;
format shortg;
%IMPORTANT: This can be used for current as well because they follow the same formulas.
%To vectorially add, just do sum(RectangularVariable) then abs(...) rad2deg(angle(...))

valFind = input('1: Find Vabc\n2: Find V123\n')


angleRad = deg2rad(120); angleRad2 = deg2rad(240);
a = exp(1i*angleRad); a2 = exp(1i*angleRad2);

alphaMatrix = [1 1 1; 1 a2 a; 1 a a2];

Vtn = zeros(3,1); Vabc = zeros(3,1); Vtn = zeros(3,1); VtnPOL = zeros(3,2); VabcPOL = zeros(3,2);
if valFind == 1
    for i=1:3
        Vtn(i) = input("vtn(" + i + '):\n');
    end
    Vabc = alphaMatrix * Vtn; % Rectangular form
    
	for i = 1:1:length(Vabc) %Polar form
    	VabcPOL(i,1) = abs(Vabc(i));
    	VabcPOL(i,2) = rad2deg(angle(Vabc(i)));
	end
disp(VabcPOL)
end

if valFind == 2
    for i=1:3
        Vabc(i) = input("vabc(" + i + '):\n');
    end
    Vtn = 1/3 * alphaMatrix * Vabc % Rectangular form
    
	for i = 1:1:length(Vtn) % Polar form
    	VtnPOL(i,1) = abs(Vtn(i));
    	VtnPOL(i,2) = rad2deg(angle(Vtn(i)));
	end
disp(VtnPOL)
end


