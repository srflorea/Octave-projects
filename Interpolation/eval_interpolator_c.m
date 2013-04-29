function N = eval_interpolator_c(tip,eps,kmax)
%functia eval_interpolator_c, determina cat de repede converge un polinom 
%si calculeaza eroarea pe care o obtine pentru un anumit suport de interpolare 

% Semnificatia parametrilor: ===============================================
%	tip  - tipul de interpolare dorit; Se poate da atat indicele cat si |
%	       numele 							    |
%       eps  - toleranta minima 					    |	
%       kmax - k-ul maxim  dupa care ciclul se opreste daca metoda nu       |
%              converge ( numarul maxim de pasi la care se refera cerinta ) |
% Rezultatul functiei : ====================================================
%	N    - numarul de noduri pentru care interpolantul converge 	    |
%       Functia va afisa de asemenea si erorile obtinute la fiecare pas.    |
% ==========================================================================

a=3;
N=1000; 
k=1;
E=0; 
E0=2*eps;
 
 if (exist('kmax')==0)        %in cazul in care parametrul optional nu exista se 
   kmax=5;  		      %initializeaza cu 5
 endif  
 
  while (abs(E-E0)>eps )
    k=k+1;
    Nk=2^(k);                
    xk=linspace(-pi,pi,Nk);   %Generam suportul de interpolare  
    y(1:Nk)=exp(a*cos(xk(1:Nk)))/(2*pi*besseli(0,a));
                              %Calculam functia in punctele din xk   
    x=linspace(-pi,pi,N+1);   %Generam punctele in care se v-a analiza functia 
    sumaE=0;                      %initializam eroarea cu zero 
    
   for l=1:N  
    val=x(l);   
    %val=pi/2;
   

    switch tip 
 
    case {"1","lagrange" }    %Daca s-a ales Lagrange se v-a calcula 
                              %polinomul de interpolare in punctul dat         
    n=Nk; 
    rez=0;
    
    %Aplicam algoritmul de calculare al polinomului Lagrange  
    
     for i=1:n 
      prod=y(i); 
          for j=1:n
            if i~=j 
             prod=prod*(val-xk(j))/(xk(i)-xk(j));
            endif
           endfor 
      rez=rez+prod; 	     %Se actualizeaza rezultatul 
     endfor 
   
     rez;
     f(val); 

  

     %In cazul 2 se va calcula polinomul de interpolare Newton 
     case { "2", "newton" }  
     n=Nk; 
     c=y;
  
      %Calculam diferentele divizate
      for i=1:n-1 
       c(i+1:n)=(c(i+1:n)-c(i))./(xk(i+1:n)-xk(i));
      endfor 
    
      %Aplicam formula lui Newton  
      rez=c(1);
      p=1;  
      for i=2:n
        p=p*(val-xk(i-1)); 
        rez=rez+p*c(i);   
      endfor 
       
     rez; 
     f(val);   
  

    %In cazul 3 se va face interpolarea cu functii spline liniare ax+b
    case {"3", "linear spline"} 
      
    n=Nk;

    %Se cauta splineul liniar caruia ii apartine punctul curent (val)
    for i=1:n-1
     if ( xk(i)>= val ||  val<=xk(i+1)) 
      ind=i;
      break; 
     endif 
    endfor 
    
    %Se calculeaza coeficientii     
    a=(y(ind+1)-y(ind))/(xk(ind+1)-xk(ind));
    b=(xk(ind+1)*y(ind)-xk(ind)*y(ind+1))/(xk(ind+1)-xk(ind));
    
    %Se calculeaza interpolantul dupa formula f(x)=a*x+b, functie care 
    %reprezinta ecuatia unei drepte 
    rez=a*val+b;
    f(val);

    



    %In cazul 4 se va face interpolarea cu functii spline cubice naturale
    case {"4", "natural"}
        
     n=Nk;
     h=abs(xk(1)-xk(2));
     A=eye(n,n);      
     b=zeros(n,1); 

     %Vom calcula coeficientii splineurilor doar o singura data 
     %pentru fiecare spatiu de interpolare generat

     if (l==1)
      for i=2:n-1 
       %Completam matricea si vectorul coloana b;
       A(i,i-1)=h; 
       A(i,i)=2*h; 
       A(i,i+1)=h; 
       b(i,1)=3*(y(i+1)-y(i))/h - 3*(y(i)- y(i-1))/h;
      endfor
     
     %Se rezolva sistemul de forma A*x=b; ci reprezinta solutia sistemului 
     ci=A\b;
     ci=ci'; 
    
     %Se calculeaza ceilalti coeficienti 
     ai(1:n)=y(1:n);
     di(1:n-1)=(ci(2:n)-ci(1:n-1))/3*h;
     bi(1:n-1)=(ai(2:n)-ai(1:n-1))/h - (2*ci(1:n-1)+ci(2:n))*h/3;    
     endif     

     %Se cauta splineul caruia ii apartine punctul 
     for i=1:n-1 
      if (val>=xk(i) && val<=xk(i+1)) 
        ind=i;  
        break;
      endif 
     endfor 
    
    %Se calculeaza valoarea interpolantului in punctul corespunzator 
    bin=val-xk(ind);
    rez=ai(ind)+bin*bi(ind)+bin^2*ci(ind)+bin^3*di(ind);
    rez;
    f(val); 




    %In cazul 5 se va face interpolarea cu splineuri cubice tensionate 
  
    %Algoritmul este destul de asemanator cu cel de la splineuri naturale 
    %cu exceptia faptului ca se pun conditiile derivatelor de ordin doi 
    %in capete 
 
    case {"5", "cubic spline"}
      
     n=Nk;
     h=abs(xk(1)-xk(2));
     A=zeros(n,n);      
     b=zeros(n,1); 

     if (l==1) 
     %Vom calcula coeficientii splineurilor doar o singura data 
     %pentru fiecare spatiu de interpolare generat


     %Punem condiitille pentru derivatele splineurilor din capete 
     A(1,1)=2*h; A(1,2)=h; A(n,n-1)=h; A(n,n)=h;
     f1d=(y(2)-y(1))/(xk(2)-xk(1));
     fNd=(y(n)-y(n-1))/(xk(n)-xk(n-1));

     b(1,1)=3*(y(2)-y(1))/h-3*f1d;
     b(n,1)=3*fNd - 3*(y(n)-y(n-1))/h;
    
     %Completam restul matricei si vectorul coloana b
      for i=2:n-1 
       A(i,i-1)=h; 
       A(i,i)=2*h; 
       A(i,i+1)=h; 
       b(i,1)=3*(y(i+1)-y(i))/h - 3*(y(i)- y(i-1))/h;
      endfor
     
     %Se rezolva sistemul A*x=b; ci reprezinta solutia sistemului 
     ci=A\b;
     ci=ci'; 
     %Apoi se calculeza restul coeficientilor dupa formulele date 
     ai(1:n)=y(1:n);
     di(1:n-1)=(ci(2:n)-ci(1:n-1))/3*h;
     bi(1:n-1)=(ai(2:n)-ai(1:n-1))/h - (2*ci(1:n-1)+ci(2:n))*h/3;    
     endif     

     %Se cauta splineul caruia ii apartine punctul pentru care dorim 
     %sa calculam interpolantul 
     for i=1:n-1 
      if (val>=xk(i) && val<=xk(i+1))
        ind=i;
        break;
      endif
     endfor

    %Se calculeaza interpolantul in punctul respectiv pentru splineul 
    %corespunzator 
    bin=val-xk(ind);
    rez=ai(ind)+bin*bi(ind)+bin^2*ci(ind)+bin^3*di(ind);
    f(val);




    %Polinomul de aproximare trigonometric Fourrier 
    case {"6", "fourrier" } 
       
    rez=0;
    n=(Nk-2)/2;

 
     for j=1:n    
       bj=0; aj=0;
       for i=1:Nk 
        bj=bj+y(i)*sin(xk(i)*j);    %Se calculeaza coeficientii bj si aj 
        aj=aj+y(i)*cos(xk(i)*j);
       endfor 
      %Se inmultesc coeficientii cu termenii corespunzatori din baza canonica 
      %si se actualizeaza rezultatul. 
      rez=rez+(1/(Nk+1))*bj*sin(val*j)+(1/(Nk+1))*aj*cos(val*j); 
     endfor   
  

    a0=0;
    for i=1:2*n+2
      a0=a0+y(i);    
    %Se calculeaza a0 separat, pentru ca formula este diferita
    %fata de cea a celorlalti coeficienti.
    endfor 
    
    rez=((1/sqrt(2))*(1/(Nk+1)))*a0 + rez;  
     %Se actualizeaza din nou rezultatul => rezultat final  
    rez; 
    f(val);

    otherwise 
    printf("Tipul ales nu exista\n");  
    endswitch  

 sumaE=sumaE+(abs(f(val)-rez))^2; %suma din formula erorii este actualizata
 
 endfor
   k                              %Se va afisa pe ecran pasul curent  
   E0=E;
   E=(((2*pi)/(N+1))*sumaE)^(1/2) %Se aplica formula finala pentru calculul
				  %erorii si va fi afisata pe ecran 
   N=Nk; 
    if (k>kmax && E>E0)             
     %In cazul in care s-a depasit un numarul de iteratii prestabilit
     %ciclul se opreste si functia returneaza inf => metoda nu converge  
        N=inf;
        break;
    endif 

 endwhile 

endfunction 

%Se defineste o functie Matlab corespunzatore lui f(x) dat in cerinta 
%si care reprezinta modul de evolutie al petelor solare 
function ans=f(x)
a=3;
ans=exp(a*cos(x))/(2*pi*besseli(0,a));
endfunction
