function [x, y] = Gen2DPoints(N,a,b,c,d)
%functie de generare a unui numar N de 
%puncte aleatoriu si uniform distribuite
x=a+(b-a)*rand(1,N); %genereaza N numere uniform distribuite in intervalul [a,b]
y=c+(d-c)*rand(1,N); %genereaza N numere unform distribuite in intervalul [c,d]
endfunction
