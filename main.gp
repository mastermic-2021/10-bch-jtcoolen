/* x is a generator of Fq obtained via ffgen(q) */
[q,t,m] = readvec("input.txt");
f = ffgen(q,'a);
encodefq(i,x)=subst(Pol(digits(i,x.p)),'x,x);
decodefq(z,x)=fromdigits(Vec(z.pol),x.p);
int2fqx(i,x)=Polrev([encodefq(z,x)|z<-digits(i,x.p^x.f)]);
fqx2int(p,x)=fromdigits([decodefq(z,x)|z<-Vecrev(p)],x.p^x.f);

/* Basile entend par "stéganographie" qu'il a ajouté à certaines positions
 * d'un code BCH correct un entier (code ASCII) correspondant aux lettres
 * du message "secret", ainsi "stéganographié" dans le code.
 *
 * Pour retrouver le message, il suffit d'identifier les positions des
 * erreurs et corriger ces erreurs. Comme le code est 35-correcteur,
 * le message est formé d'au plus 35 caractères.
 *
 * On va donc appliquer la méthode de décodage d'un code BCH vue en cours,
 * en déterminant le polynôme syndrome pour déterminer l'approximation
 * rationnelle R/E de S modulo X^(2t) avec Padé.
 * Muni de E le poly localisateur d'erreurs, on détermine les positions contenant des erreurs :
 * les indices i tels que E(alpha^(-i)) = 0.
 * La valeur de l'erreur e_i en position i étant donnée par :
 * - R(alpha^(-i)) / E'(alpha^(-i)) * (alpha^(-i * (a - 1))).
 */

y = int2fqx(m, f); \\ message reçu
syndrome(y, a, c) = sum(i=0, 2 * t - 1, subst(y, 'x, a^(i+c)) * 'x^i); \\ polynôme syndrome

/* On va rechercher le bon paramètre c satisfaisant la propriété (*):
 * alpha élément primitif (générateur de Fq) tq c(x) \in C ssi pour tout 0 <= j <= 2t-1, c(alpha^(a+j)) = 0. (*)
 * On pourrait directement calculer le polynôme générateur du code donné par le ppcm des polynômes minimaux de alpha^a,...,alpha^(a+d-2)
 * avec d=2t+1 la distance prescrite (reformulation de (*)).
 */
search_msg() = {
  \\ recherche du paramètre alpha et c :
  for(k = 0, q,
    a = ffprimroot(f); \\ on teste une nouvelle racine primitive
    for(j = 0, q - 1, \\ ordre de (F_q)*
      rational_approx = bestapprPade(Mod(syndrome(y, a, j), x^(2 * t)));
      R = numerator(rational_approx);
      E = denominator(rational_approx); \\ polynôme localisateur des erreurs
      msg = List();
      for(i = 0, q - 2, \\ taille du message : q - 1
        pol = subst(E, 'x, a^(-i));
	if (pol == 0,
	  val = subst(R/deriv(E) * (x^(j - 1)), 'x, a^(-i));
	  chr = fqx2int(val, f); \\ code ASCII de l'erreur
	  listput(msg, chr) \\ on accumule les valeurs des erreurs dans une liste pour obtenir un message potentiel
	);
      );
      candidate_msg = Strchr(Vec(msg)); \\ on convertit la liste de codes ASCII en chaîne ASCII
      if(#candidate_msg > 10, print(candidate_msg); return); \\ On s'arrête dès que le nombre d'erreur est significatif
    );
  );
};

search_msg();
