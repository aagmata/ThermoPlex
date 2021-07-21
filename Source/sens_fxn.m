function y = sens_fxn(T,MgCl,a)
Seq1 = 'ACTAGCATACTGGACATCTG';
Seq2 = seqrcomplement(Seq1);
c0t = 5e-7./a;

Xt = @(G,T,a,c0t) 0.5.*((1+a+(1/c0t)*exp(G*1000/(1.987.*(273.15+T))))...
- sqrt((1+a+(1/c0t)*exp(G*1000/(1.987.*(273.15+T)))).^2-(4.*a)));
dG = ThermoDhyb(Seq1,Seq2,T,MgCl);
y = Xt(dG,T,a,c0t);