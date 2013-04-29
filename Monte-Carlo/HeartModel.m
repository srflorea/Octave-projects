function HeartModel(file_in)
%functie ce afiseaza ariile sectiunilor iniimii calculate cu metoda a), pe o singura linie, calculate cu metoda b),pe o singura linie si valorile volumului inimii, calculat cu cele doua metode
for i=1:9 % pentru cele 9 fisiere de intrare calculam datele necesare
  file=strcat(file_in,"/","heart",num2str(i),".dat"); %concatenam calea primita ca parametru cu numele fiecarui fisier
  a=load(file); %incarcam fisierele
  x=a(:,1); % obtinem valorile coordonatei x
  y=a(:,2); % obtinem valorile coordonatei y
  z(i)=a(1,3); % obtinem valorile coordonatei z
  I=Trapez(x,y); % apelam functia Trapez
  arii1(i)=I; % depunem valorile obtinute cu funtia Trapez in vectorul arii1
  aria=MonteCarlo(x,y,0.001); % apelam functia MonteCarlo
  arii2(i)=aria; % depunem valorile obtinute cu functia MonteCarlo in vectorul arii2
endfor
z(10:18)=z(9:-1:1); % pentru calcularea volumului construim un z convenabil: primele 9 valori sunt valorile obtinute mai sus, iar urmatoarele 9 sunt aceleasi in sens invers
arii1(10:18)=0;% obtinem 2 vectori convenabili arii1 si arii22 pentru apelarea functiei Volum. 'arii1': primele 9 elemente sunt valorile ariilor calculate cu metoda a) de integrare,urmatoarele 9 sunt 0 
arii22(1:9)=0; % 'arii22': primele 9 elemente sunt 0, iar urmatoarele 9 sunt valorile ariilor calculate cu metoda b)
arii22(10:18)=arii2(1:9);
[V1,V2]=Volum(arii1,arii22,z); % se apeleaza functia ce calculeaza cele doua volume
m(1,1:9)=arii1(1:9); %se construieste o matrice cu ariile calculate cu cele 2 metode pe cate o linie
m(2,1:9)=arii2(1:9);
disp(m); % se afiseaza matricea construita
n(1)=V1; % se construieste un vector cu cele doua volume
n(2)=V2;
disp(n); %se afiseaza vectorul construit
endfunction