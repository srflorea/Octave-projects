function I2=MonteCarlo(x,y,tol)
%functia calculeaza aria unei suprafete printr-o metoda de integrare de tip Monte Carlo.
x=x(:); %liniarizam cei doi vectori
y=y(:);
nr=length(x); %obtinem lungimile acestora
[minx,imin]=min(x); %obtinem valorea minima si maxima din vectorul coordonatelor x ,precum si pozitiile lor
[maxx,imax]=max(x); 
miny=min(y); %obtinem valorea minima si maxima din vectorul coordonatelor y 
maxy=max(y);
v1(1:imax-imin+1,1:2)=[x(imin:imax) y(imin:imax)]; %impartim valorile din vectorul x in 2 vectoriastfel: valorile intre minim si maxim intr-un vector,iar celelalte in cel de-al doilea
v2(1:nr-imax,1:2)=[x(imax+1:nr) y(imax+1:nr)];
v2(nr-imax+1:nr-imax+imin-1,1:2)=[x(1:imin-1) y(1:imin-1)];
N=100;% pornim iteratia cu un N intial
nrint=0; %numarul punctelor interioare curbei este intializata la 0
%aplicam o functie de interpolare definita in matlab pentru calcularea polinomului de interpolare:aceasta este 'interp1'.
pp1=interp1(v1(:,1),v1(:,2),'pp');
pp2=interp1(v2(:,1),v2(:,2),'pp');

[a,b]=Gen2DPoints(N,minx,maxx,miny,maxy);%generam N puncte aleatoriu definite pe [minx,maxx]x[miny,maxy]
int1=ppval(pp1,a(:));%aflam valorile polinomului de interpolare in punctele generate
int2=ppval(pp2,a(:));
   DA=((b(:)>=int1(:)) & (b(:)<=int2(:))); %obtinem o matrice de dimensiunile urmatoare:nr de linii=nr de linii a lui 'a', nr de coloane=2. Va avea 0 acolo unde punctul [a,b] este in afara curbei si 1 pentru punctele [a,b] ce sunt in interiorul lui a
   nrint=sum(DA);%se inumara punctele din interiorul curbei
I=(maxy-miny)*(maxx-minx)*(nrint/N); %se calculeaza valoarea integralei dupa formula precizata in cerinta pentru 100 de puncte

N=2*N; %se dubleaza numarul de puncte ce vor fi generate dupa care se obtine o noua valoare a integralei pentru 200 de puncte la fel ca mai sus
nrint=0;
[a,b]=Gen2DPoints(N,minx,maxx,miny,maxy);
int1=ppval(pp1,a(:));
int2=ppval(pp2,a(:));
   DA=((b(:)>=int1(:)) & (b(:)<=int2(:)));
   nrint=sum(DA);
I2=(maxy-miny)*(maxx-minx)*(nrint/N);

N=2*N;
while (abs(I2-I))>tol & N<=2000000 % cat timp diferenta dintre 2 valori ale integralei este mai mare decat valoarea tol sau numaruld e puncte generate este mai mic sau egal cu 2000000 repeta algoritmul se mai sus
I=I2; 
nrint=0;
[a,b]=Gen2DPoints(N,minx,maxx,miny,maxy);
int1=ppval(pp1,a(:));
int2=ppval(pp2,a(:));
   DA=((b(:)>=int1(:)) & (b(:)<=int2(:)));
   nrint=sum(DA);
I2=(maxy-miny)*(maxx-minx)*(nrint/N);
N=2*N;
endwhile
endfunction