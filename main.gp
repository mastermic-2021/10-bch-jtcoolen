/* x is a generator of Fq obtained via ffgen(q) */
[q,t,m] = readvec("input.txt");
f = ffgen(q,'a);
encodefq(i,x)=subst(Pol(digits(i,x.p)),'x,x);
decodefq(z,x)=fromdigits(Vec(z.pol),x.p);
int2fqx(i,x)=Polrev([encodefq(z,x)|z<-digits(i,x.p^x.f)]);
fqx2int(p,x)=fromdigits([decodefq(z,x)|z<-Vecrev(p)],x.p^x.f);

\\print(q, "\n", t, "\n", m);

\\print("\n", 139*x*x^70);

y=int2fqx(m, f);
\\ distance prescrite : d = 2*35+1 = 71

\\printf("\n");
lm=digits(m,128);
\\print(length(lm));
\\a=f;
syndrome(y, a, c) = sum(i=0, 2 * t - 1, subst(y, 'x, a^(i+c)) * 'x^i);
job() = {
for(k=0, 1000,

  a = ffprimroot(f);
  for(j=0, 238,
  \\printf(ffprimroot(f));
  \\a = ffprimroot(f);
  ratio=bestapprPade(Mod(syndrome(y, a, j), x^(2*t)));
  R = numerator(ratio);
  E = denominator(ratio);
  ll = List();
  for(i = 0, 126, pol=subst(E, 'x, a^(-i)); if (pol == 0, val=subst(R/deriv(E) * (x^(j-1)), 'x, a^(-i)); l=fqx2int(val, f); listput(ll, l)));
crush=Strchr(Vec(ll));
\\print(a, " ", j);
if(crush == "Je t'aime Adele, epouse-moi !", print(crush); return);
);
);
\\syndrome(y, a) = sum(i=0, 2 * t - 1, subst(y, 'x, a^(i)) * 'x^i); \ \\alpha racine de l'unité tq c(x) in C ssi pour tout 0<=j<= 2t-1, c(alpha^(a+j)) = 0
\\print("\n");
\\printf(syndrome(y, a));
};


job();
\\locator(a) = prod(i=0, t, (1-a^i * x));

\\print("\n", lift(Mod(syndrome(y, a), x^(2*t))));
\\printf(locator(y, a));
\\ratio=bestapprPade(Mod(syndrome(y, a), x^(2*t))); \\ S = R/E
\\print("\n", ratio);
\\R = numerator(ratio);
\\E = denominator(ratio);

\\for(i = 0, 127,pol=subst(E, 'x, a^(-i)); if (pol == 0, val=subst(R/deriv(E), 'x, a^(-i)); l=fqx2int(val, f); print(Strchr(l))));

