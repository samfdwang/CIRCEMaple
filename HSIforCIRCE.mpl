
restart;

with(PDEtools);
with(plots);
declare((f, g, h, p, u, v, F, w)(seq(x[i], i = 1 .. 3)));

alias(u=u(seq(x[i],i=1..3)),v=v(seq(x[i],i=1..3)),w=w(seq(x[i],i=1..3)))
;
OrgIto:=diff(u,x[2]$2)+diff(u,x[2]$3,x[1])+6*diff(u,x[1])*diff(u,x[2])+3*diff(u,x[2],x[1])*u+3*diff(u,x[2]$2)*diff(v,x[1])
;
GBD := proc(M, p, ff, gg) local LM, n, i, j, Mul_M, k, dd, D_seq, cont, q; dd := 0; LM := nops(M); n := [seq(0, i = 1 .. LM)]; Mul_M := mul(i + 1, i in M); for i to Mul_M do for j to LM do if j = 1 then n[j] := iquo(i - 1, mul(M[k] + 1, k = j + 1 .. LM), 'r'); else n[j] := iquo(i - 1 - add(n[q]*mul(M[k] + 1, k = q + 1 .. LM), q = 1 .. j - 1), mul(M[k] + 1, k = j + 1 .. LM), 'r'); n; end if; end do; D_seq := seq(n[j], j = 1 .. LM); cont := mul(binomial(M[j], n[j]), j = 1 .. LM)*(-1)^(add(M[j] - n[j], j = 1 .. LM) mod p); dd := dd + cont*D[seq(j $ D_seq[j], j = 1 .. LM)](ff)(seq(x[i], i = 1 .. LM))*D[seq(j $ (M[j] - D_seq[j]), j = 1 .. LM)](gg)(seq(x[i], i = 1 .. LM)); convert(dd, diff); end do; end proc;
OrgBiForm:=simplify(GBD([1,3,0],2,f,f)+c[1]*GBD([1,0,1],2,f,f)+c[2]*GBD([0,2,0],2,f,f)+c[3]*GBD([0,1,1],2,f,f)+c[4]*GBD([1,1,0],2,f,f)+c[5]*GBD([0,0,2],2,f,f)+c[6]*GBD([2,0,0],2,f,f))
;
#OrgBiForm:=simplify(GBD([1,3,0],2,f,f)+GBD([1,0,1],2,f,f)+c[1]*GBD([0,2,0],2,f,f))
;
xi[1]:=sum(a[i]*x[i],i=1..3)+a[4]
;
xi[2]:=sum(b[i]*x[i],i=1..3)+b[4]
;
xi[3]:=exp(sum(p[i]*x[i],i=1..3)+p[4])
;
let_f:=sum(xi[i]^2,i=1..2)+e+xi[3]*0
;
algeq:=collect(simplify(subs(f(seq(x[i],i=1..3))=let_f,OrgBiForm)),[exp,seq(x[i],i=1..3)],'distributed'):
sys:=map(expand,[coeffs(algeq,[exp(p[1]*x[1] + p[2]*x[2] + p[3]*x[3] + p[4]),seq(x[i],i=1..3)])]):
#sys
;

#use RealDomain in sol:=solve(sys,{seq(a[i],i=1..4),seq(b[i],i=1..4),seq(c[i],i=1..6),e}) end use
;
sol:=solve(sys,{seq(a[i],i=1..4),seq(b[i],i=1..4),seq(c[i],i=1..5),e})
;
save sol, "myresult.mpl"
;
