function I=Trapez(x,y)
%functie ce intoarce valorea integralei calculata cu metoda directa a trapezului
I=0;%initializare la 0
nr=length(x);%obtinem dimensiunea vectorilor pentru parcurgere
for i=1:nr-1
I=I+((x(i+1)-x(i))/2)*(y(i+1)+y(i));%aplicam formula directa a trapezelor prezentata in curs
endfor
I=abs(I);%intoarcem integrala in valore absoluta pentru ca este o arie
endfunction