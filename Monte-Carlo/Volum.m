function [V1,V2]=Volum(arii1,arii2,z)
%functie ce calculeaza volumul functiei din fisierele de intrare
%primeste ca parametrii 3 vectori:ariile calculate cu metoda trapezelor, ariile calculate cu metodata de integrare Monte Carlo si coordonatele pe axa z
V1=Trapez(z,arii1); %se calculeaza volumul prin integrarea functiei reprezentata prin punctele z in punctele ce reprezinta ariile calculate cu prima metoda
V2=MonteCarlo(z,arii2,0.001); %se calculeaza volumul prin integrarea functiei reprezentata prin punctele z in punctele ce reprezinta ariile calculate cu a doua metoda
endfunction