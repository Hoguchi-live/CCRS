/*

This file contains the formulas for radical N-isogenies for N = 2, ..., 13, 15.
More precisely, let E be an elliptic curve in a general form such that P = (0,0) is
an N-torsion point. We give formulas for an isomorpic form Eiso of EA = E/<P> such
that Eiso has the same form as E and (0,0) on Eiso extends E -> EA cyclically.

For N = 2 we use the Montgomery form of an elliptic curve.
For N = 3 we use a Weierstrass form that is isomorphic to a Hessian form.
For N > 3 we use the Tate normal form of an elliptic curve. From N = 6 onwards, we use
the optimized equations by Sutherland that can be found here:
http://math.mit.edu/~drew/X1_optcurves.html

Remark that instead of using the assertion IsIsomorphic(Eiso,EA) we use the assertion
jInvariant(EA) eq jInvariant(Eiso), or - starting from N = 4 - a self-written function
to compare the j-invariant of an elliptic curve with the j-invariant of a curve from
just the b and c parameters of the Tate normal form. The reason for the self-written
function is that it avoids inversion. The reason for avoiding IsIsomorphic is that it
is too hard to compute for larger values of N.

*/


clear;

function same_j_invars(EA, newb, newc)

  a_1 := 1-newc; a_2 := -newb; a_3 := -newb; a_4 := 0; a_6 := 0;
  b_2 := a_1^2 + 4*a_2;
  b_4 := a_1*a_3 + 2*a_4;
  b_6 := a_3^2 + 4*a_6;
  b_8 := a_1^2*a_6 + 4*a_2*a_6 - a_1*a_3*a_4 + a_2*a_3^2 - a_4^2;
  c_4 := b_2^2 - 24*b_4;
  c_6 := -b_2^3 + 36*b_2*b_4 - 216*b_6;
  d := -b_2^2*b_8 - 8*b_4^3 - 27*b_6^2 + 9*b_2*b_4*b_6;
  return jInvariant(EA)*d eq c_4^3;

end function;



N := 2;
QA<A> := FunctionField(Rationals());
QAB<B> := FunctionField(QA);
E := EllipticCurve([0,A,0,B,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QAz<z> := PolynomialRing(QAB);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QAB | z^N - t>;
newA := -3*w + A;
newB := -2*A*w + 8*B;
Eiso := EllipticCurve([0,newA,0,newB,0]);
assert jInvariant(EA) eq jInvariant(Eiso);
"Formulas verified for N equal", N;



N := 3;
QA<A> := FunctionField(Rationals());
QAB<B> := FunctionField(QA);
E := EllipticCurve([A,0,B,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QAz<z> := PolynomialRing(QAB);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QAB | z^N - t>;
newA := -6*w + A;
newB := 3*A*w^2 - A^2*w + 9*B;
Eiso := EllipticCurve([newA,0,newB,0,0]);
assert jInvariant(EA) eq jInvariant(Eiso);
"Formulas verified for N equal", N;



N := 4;
QA<A> := FunctionField(Rationals());
b := -A; c := 0;
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QAz<z> := PolynomialRing(QA);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QA | z^N - t>;
numA := w*(4*w^2+1);
denA := (2*w+1)^4;
newA := numA/denA;
newb := -newA; newc := 0;
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 5;
QA<A> := FunctionField(Rationals());
b := A; c := A;
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QAz<z> := PolynomialRing(QA);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QA | z^N - t>;
numA := w * (w^4 + 3*w^3 + 4*w^2 + 2*w + 1);
denA := (w^4 - 2*w^3 + 4*w^2 - 3*w + 1);
newA := numA/denA;
newb := newA; newc := newA;
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 6;
QA<A> := FunctionField(Rationals());
r := A; s := 1;
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QAz<z> := PolynomialRing(QA);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QA | z^N - t>;
numA := (-3*A + 2)*w^4 + 3*A^2*w^2 + 2*A*w - 3*A^3 + 4*A^2;
denA := w^4 + 2*A*w^2 + 3*A*w + A^2;
newA := numA/denA;
newr := newA; news := 1;
newb := newr*news*(newr-1); newc := news*(newr-1);
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 7;
QA<A> := FunctionField(Rationals());
r := A; s := A;
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QAz<z> := PolynomialRing(QA);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QA | z^N - t>;
numA := w^6 + A*w^5 + 2*A^3*w^2 - A^3*w + A^4;
denA := -w^6 + A*w^4 + A^3*w^2 - 2*A^3*w + A^4;
newA := numA/denA;
newr := newA; news := newA;
newb := newr*news*(newr-1); newc := news*(newr-1);
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 8;
QA<A> := FunctionField(Rationals());
r := 1/(2-A); s := A;
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QAz<z> := PolynomialRing(QA);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QA | z^N - t>;
numA := -2*A*(A - 2)*w^2 - A*(A - 2);
denA := (A - 2)^2*w^4 - A*(A-2)*w^2 - A*(A-2)*w + A;
newA := numA/denA;
newr := 1/(2-newA); news := newA;
newb := newr*news*(newr-1); newc := news*(newr-1);
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 9;
QA<A> := FunctionField(Rationals());
r := A^2-A+1; s := A;
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QAz<z> := PolynomialRing(QA);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QA | z^N - t>;
numA := ((A*(A^2 - A + 1))*((w^3 + A*(A^2 - A + 1))*w^2 + (A*(A^2 - A + 1))^2));
denA := (((w^3 - A*(A^2 - A + 1)*(A-1))*w^3 - A^3*(A^2 - A + 1)^2)*w + (A*(A^2 - A + 1))^3);
newA := numA/denA;
newr := newA^2-newA+1; news := newA;
newb := newr*news*(newr-1); newc := news*(newr-1);
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 10;
QA<A> := FunctionField(Rationals());
r := -A^2/(A^2-3*A+1); s := A;
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QAz<z> := PolynomialRing(QA);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QA | z^N - t>;
numA := (A^6 - 5*A^5 + 5*A^4 + 5*A^3 - 5*A^2 + A)*w^6 + (2*A^6 - 14*A^5 + 34*A^4 -
    34*A^3 + 14*A^2 - 2*A)*w^5 + (-2*A^6 + 7*A^5 - 3*A^4 - 6*A^3 + 5*A^2 -
    A)*w^4 + (-4*A^6 + 16*A^5 - 17*A^4 + 7*A^3 - A^2)*w^3 + (6*A^6 - 7*A^5 +
    2*A^4)*w^2 + (2*A^7 + 3*A^6 - 4*A^5 + A^4)*w;
denA := (A^6 - 2*A^5 - 15*A^4 + 50*A^3 - 45*A^2 + 16*A - 2)*w^6 + (3*A^6 - 18*A^5 +
    32*A^4 - 12*A^3 - 8*A^2 + 6*A - 1)*w^5 + (2*A^5 - 5*A^4 - 2*A^3 + 4*A^2 -
    A)*w^4 + (-4*A^6 + 18*A^5 - 24*A^4 + 12*A^3 - 2*A^2)*w^3 + (8*A^6 - 6*A^5 +
    A^4)*w^2 + (6*A^7 - 5*A^6 + A^5)*w + 2*A^7 - 3*A^6 + A^5;
newA := numA/denA;
newr := -newA^2/(newA^2-3*newA+1); news := newA;
newb := newr*news*(newr-1); newc := news*(newr-1);
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 11;
QA<A> := FunctionField(Rationals());
QAB<B> := PolynomialRing(QA);
modpol := B^2 + (A^2 + 1)*B + A;
QB<B> := ext<QA | modpol>;
r := A*B+1; s := -A+1;
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QBz<z> := PolynomialRing(QB);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QB | z^N - t>;
numA := (3*A + 2*B)*w^7 + (2*A - B^2 + 1)*w^6 + ((-7*B^2 - 5*B)*A + (5*B^2 + 3*B))*w^5 + ((-B^2 - B)*A - 2*B^3 - 3*B^2 - B)*w^4 + ((3*B^4 + B^3 + B^2)*A - 4*B^4 - B^3 + 2*B^2 + 2*B)*w^3 + ((-B^5 -
    6*B^4 - 6*B^3)*A - B^5 - B^4 - 2*B^3 - 3*B^2)*w^2 + (-B^5*A + (2*B^6 + 10*B^5 + 14*B^4 + 9*B^3 + 2*B^2))*w + (-4*B^5 - 6*B^4 - 2*B^3)*A - 2*B^6 - 2*B^5;
denA := 3*w^8 + (3*A - 3*B - 5)*w^7 + (2*B*A + (B^2 + 7*B + 4))*w^6 + ((-3*B^2 - 2*B)*A - 2*B^2 - 3*B)*w^5 + ((-5*B^2 - 2*B)*A - 2*B^3 + B^2)*w^4 + ((B^4 - 5*B^3 - B^2)*A - 5*B^4 - 2*B^3 - B^2 +
    B)*w^3 + ((-5*B^4 - 3*B^3 - B^2)*A + (2*B^5 + 14*B^4 + 15*B^3 + 6*B^2))*w^2 + ((-7*B^5 - 11*B^4 - 5*B^3)*A - 6*B^6 - 15*B^5 - 20*B^4 - 15*B^3 - 5*B^2)*w + (7*B^5 + 9*B^4 + 3*B^3)*A +
    2*B^7 + 13*B^6 + 21*B^5 + 18*B^4 + 9*B^3 + 2*B^2;
newA := numA/denA;
numB := w^8 + (4*A + B)*w^7 + (3*B^2 - 2*B + 2)*w^6 + ((-5*B^2 - 4*B)*A + (3*B^3 + 2*B^2))*w^5 + ((B^3 + 12*B^2)*A - B^3 - 10*B^2 + 2*B)*w^4 + ((-B^4 + 3*B^3 - 2*B^2)*A + (9*B^4 + 13*B^3 +
    10*B^2))*w^3 + ((6*B^4 - 3*B^3 - B^2)*A + (5*B^5 + B^4 + 2*B^3 - 2*B^2))*w^2 + ((8*B^5 + 12*B^4 + 9*B^3)*A - 2*B^6 - 15*B^5 - 16*B^4 - 4*B^3 + 4*B^2)*w + (-4*B^5 + 6*B^4 + 2*B^3)*A -
    2*B^7 - 6*B^6 + 6*B^5 + 8*B^4 + 6*B^3;
denB := -w^8 - 6*A*w^7 + (2*B*A - 6*B^2 - 4*B + 2)*w^6 + ((3*B^2 + 7*B)*A + (9*B^2 + 5*B))*w^5 + ((-10*B^2 - B)*A - 3*B^3 - 6*B)*w^4 + ((-6*B^4 - 19*B^3 - 4*B^2)*A - 9*B^4 - 13*B^3 - 15*B^2 -
    2*B)*w^3 + ((3*B^5 + 2*B^4 + 7*B^3 + 3*B^2)*A - 7*B^5 - 9*B^4 - 9*B^3 - 2*B^2)*w^2 + ((-10*B^5 + 3*B^4)*A + (3*B^6 + 32*B^5 + 32*B^4 + 16*B^3))*w + (2*B^5 - 10*B^4 - 4*B^3)*A + 4*B^7 +
    14*B^6 + 8*B^5 + 4*B^4 - 2*B^3;
newB := numB/denB;
newr := newA*newB+1; news := -newA+1;
newb := newr*news*(newr-1); newc := news*(newr-1);
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 12;
QA<A> := FunctionField(Rationals());
r := (2*A^2-2*A+1)/A; s := (3*A^2-3*A+1)/A^2;
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QAz<z> := PolynomialRing(QA);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QA | z^N - t>;
numA := A^5*w^7 + (A^6 - A^5)*w^6 + (4*A^6 - 6*A^5 + 4*A^4 - A^3)*w^5 + (4*A^7 - 10*A^6 + 10*A^5 - 5*A^4 + A^3)*w^4 + (12*A^7 - 42*A^6 + 64*A^5 - 55*A^4 + 28*A^3 - 8*A^2 + A)*w^3 - 24*A^9 + 84*A^8
    - 140*A^7 + 140*A^6 - 90*A^5 + 37*A^4 - 9*A^3 + A^2;
denA := A^5*w^7 + (A^5 - A^4)*w^6 + (4*A^6 - 6*A^5 + 4*A^4 - A^3)*w^5 + (8*A^7 - 20*A^6 + 20*A^5 - 10*A^4 + 2*A^3)*w^4 + (24*A^7 - 84*A^6 + 128*A^5 - 110*A^4 + 56*A^3 - 16*A^2 + 2*A)*w^3 + (24*A^7
    - 84*A^6 + 140*A^5 - 140*A^4 + 90*A^3 - 37*A^2 + 9*A - 1)*w^2 + (-24*A^8 + 84*A^7 - 140*A^6 + 140*A^5 - 90*A^4 + 37*A^3 - 9*A^2 + A)*w - 48*A^9 + 192*A^8 - 364*A^7 + 420*A^6 - 320*A^5 +
    164*A^4 - 55*A^3 + 11*A^2 - A;
newA := numA/denA;
newr := (2*newA^2-2*newA+1)/newA; news := (3*newA^2-3*newA+1)/newA^2;
newb := newr*news*(newr-1); newc := news*(newr-1);
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 13;
QA<A> := FunctionField(Rationals());
QAB<B> := PolynomialRing(QA);
modpol := B^2 + (A^3 + A^2 + 1)*B - A^2 - A;
QB<B> := ext<QA | modpol>;
r := -A*B + 1; s := 1 - A*B/(B + 1);
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QBz<z> := PolynomialRing(QB);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QB | z^N - t>;
numA := ((A^6 - A^4 + 2*A^3 + A^2 - 2*A + 1)*B + (A^9 + A^8 - A^7 + 2*A^6 + 4*A^5 - 2*A^4 + 3*A^2 - A))*w^9 + ((-A^5 - A^4 + A^3 - A^2 - 2*A + 1)*B - A^8 - 2*A^7 - A^5 - 5*A^4 - A^3 + A^2 -
    2*A)*w^8 + (A^2*B + (A^5 + A^4 + A^2 + A))*w^7 + (A*B + (A^4 + A^3 + A))*w^6 + (-A*B - A^4 - A^3 - A^2 - A)*w^5 + (A^3 + A^2)*w^4 + ((A^4 - A^3 + A^2)*B - A^3 + A^2)*w^3 + ((-A^5 + A^3
    - A^2)*B + (A^4 - A^2))*w^2 + (A^9 - A^7 + 2*A^6 + A^5 - 2*A^4 + A^3)*B - A^8 + A^6 - A^5 - A^4 + A^3;
denA := ((-A^6 + A^4 - 2*A^3 - A^2 + 2*A - 1)*B - A^9 - A^8 + A^7 - 2*A^6 - 4*A^5 + 2*A^4 - 3*A^2 + A)*w^9 + ((A^3 + 1)*B + (A^6 + A^5 + 2*A^3 + 2*A^2 + 1))*w^8 + ((A^3 - A^2 - A + 1)*B + (A^6 -
    2*A^4 + A^3 + A^2 - 2*A))*w^7 + (A*B + (A^4 + A^3 + A))*w^5 - A^3*w^4 + ((A^3 - A^2)*B - A^2)*w^3 + ((A^6 + A^5 - A^4 + A^2)*B - A^5 - A^4 + A^3 + A^2)*w^2 + ((A^7 + A^4)*B - A^6)*w +
    (-A^9 + A^7 - 2*A^6 - A^5 + 2*A^4 - A^3)*B + A^8 - A^6 + A^5 + A^4 - A^3;
newA := numA/denA;
numB := ((-A^7 - 3*A^6 - A^5 - A^4 - 6*A^3 - 3*A^2 + A - 2)*B - A^10 - 4*A^9 - 4*A^8 - 3*A^7 - 11*A^6 - 13*A^5 - 4*A^4 - 7*A^3 - 8*A^2 - A - 1)*w^10 + ((2*A^6 + 2*A^5 - 3*A^4 + 3*A^3 + 4*A^2 - 3*A
    + 1)*B + (2*A^9 + 4*A^8 - A^7 + 2*A^6 + 11*A^5 - 2*A^3 + 6*A^2 - A + 1))*w^9 + ((-A^7 - 2*A^6 + A^5 - A^4 - 2*A^3 + A^2 + 2)*B - A^10 - 3*A^9 - A^8 - A^7 - 6*A^6 - 2*A^5 + A^4 + 3*A^2 +
    2)*w^8 + ((A^5 - 2*A^4 - 3*A)*B + (A^8 - A^7 - 2*A^6 + A^5 - 4*A^4 - 5*A^3 - A^2 - 3*A))*w^7 + ((-A^5 - 2*A^4 - 2*A^3 + 2*A^2 + A)*B - A^8 - 3*A^7 - 4*A^6 - A^5 + A^4 + 3*A^2 + A)*w^6 +
    ((-3*A^3 + 2*A^2)*B - 3*A^5 - A^4 - A^3 + 2*A^2)*w^5 + ((-A^6 + A^4 + 2*A^2)*B + (2*A^5 + 3*A^4 + 2*A^3 + 2*A^2))*w^4 + ((-A^8 - A^7 + A^6 - A^5 + A^4 - A^2)*B + (A^7 + A^6 - A^5 - A^3
    - A^2))*w^3 + ((2*A^10 + A^9 - A^8 + 2*A^7 + 3*A^6 - A^5 - A^4 + A^3)*B - 2*A^9 - A^8 + A^7 - 2*A^5 + A^3)*w^2 + ((A^10 - A^9 - 3*A^8 + 3*A^6 - 2*A^5 - 3*A^4 + A^3 + A^2)*B - A^9 + A^8
    + 3*A^7 + A^6 - 4*A^5 - A^4 + 2*A^3 + A^2)*w + (-A^10 + 2*A^8 - 2*A^7 - 2*A^6 + 3*A^5 - A^3)*B + A^9 - 2*A^7 + A^6 + 2*A^5 - A^4 - A^3;
denB := ((4*A^6 + 2*A^5 - 2*A^4 + 9*A^3 + 6*A^2 - 4*A + 5)*B + (4*A^9 + 6*A^8 + 11*A^6 + 21*A^5 + 2*A^4 + 8*A^3 + 16*A^2 + 3))*w^10 + ((2*A^6 + 2*A^5 + 2*A^3 + 4*A^2)*B + (2*A^9 + 4*A^8 + 2*A^7 +
    4*A^6 + 10*A^5 + 6*A^4 + 2*A^3 + 4*A^2 + 2*A))*w^9 + ((A^7 - 2*A^5 + A^4)*B + (A^10 + A^9 - 2*A^8 + 2*A^6 - 2*A^5 - A^4 + 2*A))*w^8 + ((2*A^6 + A^5 + 2*A^4 + A^3 - A^2 + 2*A)*B + (2*A^9
    + 3*A^8 + 3*A^7 + 5*A^6 + 3*A^5 + 4*A^4 + 5*A^3 + 2*A))*w^7 + ((A^5 + 4*A^4 + 2*A^3 + 2*A^2)*B + (A^8 + 5*A^7 + 6*A^6 + 5*A^5 + 9*A^4 + 4*A^3 + 2*A^2))*w^6 + ((-A^4 + 3*A^3 - 4*A^2 -
    3*A)*B - A^7 + 3*A^5 - 2*A^4 - 3*A^3 - 7*A^2 - 3*A)*w^5 + ((-2*A^5 - A^3 - A^2 - A)*B - 4*A^5 - A^4 - 3*A^3 - 2*A^2 - A)*w^4 + ((-A^8 + 2*A^7 + 2*A^6 - A^5 - A^4 - A^3)*B + (A^7 - 2*A^6
    - 2*A^5 - 2*A^4 - A^3))*w^3 + ((-3*A^10 - A^9 + A^8 - 8*A^7 - 2*A^6 + 5*A^5 - 3*A^4 - 2*A^3 + A^2)*B + (3*A^9 + A^8 - A^7 + 5*A^6 + A^5 - 4*A^4 - A^3 + A^2))*w^2 + ((-A^10 + A^9 + A^8 +
    A^3)*B + (A^9 - A^8 - A^7 - A^6 + A^5 + A^4 + A^3))*w + (3*A^10 - 2*A^9 - 6*A^8 + 8*A^7 + 4*A^6 - 11*A^5 + 2*A^4 + 3*A^3)*B - 3*A^9 + 2*A^8 + 6*A^7 - 5*A^6 - 6*A^5 + 5*A^4 + 3*A^3;
newB := numB/denB;
newr := -newA*newB + 1; news := 1 - newA*newB/(newB + 1);
newb := newr*news*(newr-1); newc := news*(newr-1);
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 14;
QA<A> := FunctionField(Rationals());
QAB<B> := PolynomialRing(QA);
modpol := B^2 + (A^2 + A)*B + A;
QB<B> := ext<QA | modpol>;
r := 1 - (A + B)/((B + 1)*(A + B + 1));
s := (1 - A)/(B + 1);
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QBz<z> := PolynomialRing(QB);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QB | z^N - t>;
newA := ((-A^15 - 12*A^14 - 59*A^13 - 164*A^12 - 297*A^11 - 320*A^10 - 18*A^9 + 342*A^8 + 24*A^7 - 204*A^6 + 720*A^5 + 692*A^4 - 634*A^3 - 384*A^2 + 265*A - 14)/(A^8 - 2*A^7 - 10*A^6 + 30*A^5 -
    16*A^4 - 22*A^3 + 26*A^2 - 6*A - 1)*B + (-A^17 - 13*A^16 - 72*A^15 - 225*A^14 - 437*A^13 - 520*A^12 - 275*A^11 + 208*A^10 + 607*A^9 + 757*A^8 + 412*A^7 - 302*A^6 - 157*A^5 + 389*A^4 -
    153*A^3 - 357*A^2 + 76*A - 1)/(A^8 - 2*A^7 - 10*A^6 + 30*A^5 - 16*A^4 - 22*A^3 + 26*A^2 - 6*A - 1))*w^12 + ((2*A^13 + 21*A^12 + 96*A^11 + 271*A^10 + 520*A^9 + 550*A^8 + 67*A^7 - 88*A^6
    + 700*A^5 + 382*A^4 - 1025*A^3 - 529*A^2 + 328*A - 15)/(A^7 - A^6 - 11*A^5 + 19*A^4 + 3*A^3 - 19*A^2 + 7*A + 1)*B + (2*A^15 + 24*A^14 + 123*A^13 + 363*A^12 + 695*A^11 + 893*A^10 +
    745*A^9 + 323*A^8 - 254*A^7 - 729*A^6 - 306*A^5 + 297*A^4 - 309*A^3 - 482*A^2 + 88*A - 1)/(A^7 - A^6 - 11*A^5 + 19*A^4 + 3*A^3 - 19*A^2 + 7*A + 1))*w^11 + ((-A^12 - 11*A^11 - 40*A^10 -
    64*A^9 - 73*A^8 - 118*A^7 - 31*A^6 + 276*A^5 + 178*A^4 - 231*A^3 - 53*A^2 + 148*A - 12)/(A^7 - A^6 - 11*A^5 + 19*A^4 + 3*A^3 - 19*A^2 + 7*A + 1)*B + (-A^14 - 11*A^13 - 49*A^12 -
    116*A^11 - 158*A^10 - 101*A^9 + 57*A^8 + 152*A^7 + 83*A^6 + 115*A^5 + 161*A^4 - 92*A^3 - 124*A^2 + 53*A - 1)/(A^7 - A^6 - 11*A^5 + 19*A^4 + 3*A^3 - 19*A^2 + 7*A + 1))*w^10 + ((-A^12 -
    7*A^11 - 12*A^10 + 6*A^9 + 39*A^8 + 105*A^7 + 155*A^6 - 101*A^5 - 248*A^4 + 186*A^3 + 109*A^2 - 221*A + 22)/(A^7 - A^6 - 11*A^5 + 19*A^4 + 3*A^3 - 19*A^2 + 7*A + 1)*B + (-A^14 - 8*A^13
    - 20*A^12 - 5*A^11 + 67*A^10 + 154*A^9 + 144*A^8 - A^7 - 65*A^6 - 87*A^5 - 230*A^4 + 2*A^3 + 167*A^2 - 87*A + 2)/(A^7 - A^6 - 11*A^5 + 19*A^4 + 3*A^3 - 19*A^2 + 7*A + 1))*w^9 + ((A^8 +
    15*A^7 + 49*A^6 + 50*A^5 + 11*A^4 + 19*A^3 + 11*A^2 - 20*A + 8)/(A^6 - 11*A^4 + 8*A^3 + 11*A^2 - 8*A - 1)*B + (2*A^10 + 17*A^9 + 54*A^8 + 85*A^7 + 78*A^6 + 33*A^5 - 33*A^4 - 37*A^3 -
    6*A^2 - 18*A + 1)/(A^6 - 11*A^4 + 8*A^3 + 11*A^2 - 8*A - 1))*w^7 + ((-A^7 - 6*A^6 - 4*A^5 + 15*A^4 + 6*A^3 - 15*A^2 + 7*A + 6)/(A^5 - A^4 - 10*A^3 + 18*A^2 - 7*A - 1)*B + (-A^9 - 6*A^8
    - 9*A^7 + 2*A^6 + 14*A^5 + 20*A^4 + 10*A^3 - 17*A^2 - 6*A + 1)/(A^5 - A^4 - 10*A^3 + 18*A^2 - 7*A - 1))*w^5 + ((A^5 + 6*A^4 + 11*A^3 + 4*A^2 + 6)/(A^4 - 10*A^2 + 8*A + 1)*B + (A^7 +
    7*A^6 + 16*A^5 + 15*A^4 + 6*A^3 - 3*A^2 - 7*A + 1)/(A^4 - 10*A^2 + 8*A + 1))*w^4 + ((-4*A^3 - 4*A^2 + 8*A + 4)/(A^4 - 10*A^2 + 8*A + 1)*B + (-A^6 - 4*A^5 + 7*A^3 + A + 1)/(A^4 - 10*A^2
    + 8*A + 1))*w^3 + ((-A^3 + A^2 + A - 5)/(A^4 - 10*A^2 + 8*A + 1)*B + (-A^4 - 4*A^3 - 2*A^2 + 4*A - 1)/(A^4 - 10*A^2 + 8*A + 1))*w^2 + (A - 1)/(A^3 + A^2 - 9*A - 1)*B + (A^2 + 3*A)/(A^3
    + A^2 - 9*A - 1);
newB := ((-5*A^17 - 80*A^16 - 564*A^15 - 2321*A^14 - 6218*A^13 - 11253*A^12 - 13617*A^11 - 11417*A^10 - 8960*A^9 - 3429*A^8 + 15382*A^7 + 23069*A^6 - 4850*A^5 - 15487*A^4 + 10023*A^3 + 10477*A^2 -
    2583*A + 73)/(A^10 - 21*A^8 + 16*A^7 + 122*A^6 - 176*A^5 - 58*A^4 + 176*A^3 - 43*A^2 - 16*A - 1)*B + (-5*A^19 - 86*A^18 - 650*A^17 - 2874*A^16 - 8345*A^15 - 16824*A^14 - 23872*A^13 -
    22724*A^12 - 10755*A^11 + 6456*A^10 + 17387*A^9 + 14998*A^8 + 11477*A^7 + 15348*A^6 + 4742*A^5 - 11676*A^4 - 2108*A^3 + 5730*A^2 - 543*A + 4)/(A^10 - 21*A^8 + 16*A^7 + 122*A^6 - 176*A^5
    - 58*A^4 + 176*A^3 - 43*A^2 - 16*A - 1))*w^13 + ((-A^17 - 9*A^16 - 5*A^15 + 186*A^14 + 904*A^13 + 2398*A^12 + 4490*A^11 + 4425*A^10 - 1726*A^9 - 5635*A^8 + 3195*A^7 + 4632*A^6 -
    9724*A^5 - 6656*A^4 + 6344*A^3 + 2389*A^2 - 1429*A + 62)/(A^10 - 21*A^8 + 16*A^7 + 122*A^6 - 176*A^5 - 58*A^4 + 176*A^3 - 43*A^2 - 16*A - 1)*B + (-A^19 - 10*A^18 - 18*A^17 + 173*A^16 +
    1173*A^15 + 3504*A^14 + 6190*A^13 + 6420*A^12 + 2588*A^11 - 2556*A^10 - 6577*A^9 - 8791*A^8 - 2775*A^7 + 5904*A^6 + 380*A^5 - 4758*A^4 + 1719*A^3 + 2158*A^2 - 375*A + 4)/(A^10 - 21*A^8
    + 16*A^7 + 122*A^6 - 176*A^5 - 58*A^4 + 176*A^3 - 43*A^2 - 16*A - 1))*w^12 + ((2*A^16 + 24*A^15 + 130*A^14 + 378*A^13 + 510*A^12 - 224*A^11 - 1738*A^10 - 1292*A^9 + 1602*A^8 - 804*A^7 -
    7014*A^6 - 418*A^5 + 8414*A^4 + 284*A^3 - 3794*A^2 + 900*A - 32)/(A^10 - 21*A^8 + 16*A^7 + 122*A^6 - 176*A^5 - 58*A^4 + 176*A^3 - 43*A^2 - 16*A - 1)*B + (2*A^18 + 28*A^17 + 160*A^16 +
    476*A^15 + 736*A^14 + 282*A^13 - 1190*A^12 - 2826*A^11 - 3372*A^10 - 2196*A^9 + 1086*A^8 + 3100*A^7 - 1236*A^6 - 2898*A^5 + 3530*A^4 + 2418*A^3 - 1890*A^2 + 208*A - 2)/(A^10 - 21*A^8 +
    16*A^7 + 122*A^6 - 176*A^5 - 58*A^4 + 176*A^3 - 43*A^2 - 16*A - 1))*w^11 + ((-5*A^15 - 61*A^14 - 288*A^13 - 704*A^12 - 1006*A^11 - 526*A^10 + 1545*A^9 + 2735*A^8 - 1223*A^7 - 2287*A^6 +
    5132*A^5 + 3272*A^4 - 5150*A^3 - 1998*A^2 + 1123*A - 47)/(A^10 - 21*A^8 + 16*A^7 + 122*A^6 - 176*A^5 - 58*A^4 + 176*A^3 - 43*A^2 - 16*A - 1)*B + (-5*A^17 - 66*A^16 - 354*A^15 - 999*A^14
    - 1564*A^13 - 1048*A^12 + 922*A^11 + 2886*A^10 + 3214*A^9 + 2541*A^8 - 272*A^7 - 4061*A^6 - 662*A^5 + 2888*A^4 - 1440*A^3 - 1754*A^2 + 289*A - 3)/(A^10 - 21*A^8 + 16*A^7 + 122*A^6 -
    176*A^5 - 58*A^4 + 176*A^3 - 43*A^2 - 16*A - 1))*w^10 + ((A^15 + 16*A^14 + 106*A^13 + 342*A^12 + 592*A^11 + 694*A^10 + 617*A^9 - 532*A^8 - 1993*A^7 + 26*A^6 + 1850*A^5 - 1212*A^4 -
    1136*A^3 + 1048*A^2 - 165*A + 2)/(A^10 - 21*A^8 + 16*A^7 + 122*A^6 - 176*A^5 - 58*A^4 + 176*A^3 - 43*A^2 - 16*A - 1)*B + (A^17 + 17*A^16 + 119*A^15 + 446*A^14 + 974*A^13 + 1234*A^12 +
    652*A^11 - 634*A^10 - 1368*A^9 - 977*A^8 - 1039*A^7 - 492*A^6 + 1508*A^5 + 334*A^4 - 948*A^3 + 456*A^2 - 27*A)/(A^10 - 21*A^8 + 16*A^7 + 122*A^6 - 176*A^5 - 58*A^4 + 176*A^3 - 43*A^2 -
    16*A - 1))*w^9 + ((-A^13 - 23*A^12 - 153*A^11 - 487*A^10 - 909*A^9 - 1336*A^8 - 1521*A^7 + 219*A^6 + 2678*A^5 + 530*A^4 - 2386*A^3 - 408*A^2 + 628*A - 31)/(A^9 + A^8 - 20*A^7 - 4*A^6 +
    118*A^5 - 58*A^4 - 116*A^3 + 60*A^2 + 17*A + 1)*B + (-A^15 - 21*A^14 - 167*A^13 - 675*A^12 - 1547*A^11 - 2069*A^10 - 1389*A^9 + 157*A^8 + 806*A^7 + 390*A^6 + 1140*A^5 + 1036*A^4 -
    1070*A^3 - 736*A^2 + 180*A - 2)/(A^9 + A^8 - 20*A^7 - 4*A^6 + 118*A^5 - 58*A^4 - 116*A^3 + 60*A^2 + 17*A + 1))*w^8 + ((2*A^12 + 25*A^11 + 147*A^10 + 419*A^9 + 538*A^8 + 303*A^7 +
    371*A^6 + 395*A^5 - 533*A^4 - 588*A^3 + 254*A^2 + 86*A - 11)/(A^9 + A^8 - 20*A^7 - 4*A^6 + 118*A^5 - 58*A^4 - 116*A^3 + 60*A^2 + 17*A + 1)*B + (2*A^14 + 29*A^13 + 171*A^12 + 530*A^11 +
    935*A^10 + 952*A^9 + 438*A^8 - 315*A^7 - 555*A^6 - 34*A^5 - 108*A^4 - 435*A^3 + 14*A^2 + 41*A - 1)/(A^9 + A^8 - 20*A^7 - 4*A^6 + 118*A^5 - 58*A^4 - 116*A^3 + 60*A^2 + 17*A + 1))*w^7 +
    ((-13*A^8 - 126*A^7 - 357*A^6 - 245*A^5 + 84*A^4 - 420*A^3 - 525*A^2 + 215*A - 21)/(A^7 + A^6 - 19*A^5 - 3*A^4 + 99*A^3 - 61*A^2 - 17*A - 1)*B + (-A^11 - 19*A^10 - 126*A^9 - 386*A^8 -
    608*A^7 - 581*A^6 - 315*A^5 + 226*A^4 + 171*A^3 - 294*A^2 + 79*A - 2)/(A^7 + A^6 - 19*A^5 - 3*A^4 + 99*A^3 - 61*A^2 - 17*A - 1))*w^6 + ((-2*A^8 + 6*A^7 + 48*A^6 + 89*A^5 + 73*A^4 +
    76*A^3 + 166*A^2 + 133*A - 13)/(A^7 + A^6 - 19*A^5 - 3*A^4 + 99*A^3 - 61*A^2 - 17*A - 1)*B + (3*A^9 + 32*A^8 + 125*A^7 + 242*A^6 + 246*A^5 + 107*A^4 - 31*A^3 - 12*A^2 + 57*A - 1)/(A^7 +
    A^6 - 19*A^5 - 3*A^4 + 99*A^3 - 61*A^2 - 17*A - 1))*w^5 + ((A^8 + 5*A^7 + 25*A^6 + A^5 - 117*A^4 + 39*A^3 + 211*A^2 - 77*A + 8)/(A^7 + A^6 - 19*A^5 - 3*A^4 + 99*A^3 - 61*A^2 - 17*A -
    1)*B + (A^10 + 9*A^9 + 22*A^8 - 5*A^7 - 51*A^6 + 31*A^5 + 41*A^4 - 11*A^3 + 146*A^2 - 24*A + 1)/(A^7 + A^6 - 19*A^5 - 3*A^4 + 99*A^3 - 61*A^2 - 17*A - 1))*w^4 + ((-A^7 - 15*A^6 - 43*A^5
    - 7*A^4 + 81*A^3 - 5*A^2 - 85*A + 11)/(A^7 + A^6 - 19*A^5 - 3*A^4 + 99*A^3 - 61*A^2 - 17*A - 1)*B + (-2*A^9 - 17*A^8 - 44*A^7 - 26*A^6 + 14*A^5 - 8*A^4 + 24*A^3 + 34*A^2 - 40*A +
    1)/(A^7 + A^6 - 19*A^5 - 3*A^4 + 99*A^3 - 61*A^2 - 17*A - 1))*w^3 + ((A^6 + 26*A^5 + 49*A^4 - 24*A^3 - 97*A^2 + 14*A - 1)/(A^7 + A^6 - 19*A^5 - 3*A^4 + 99*A^3 - 61*A^2 - 17*A - 1)*B +
    (3*A^8 + 25*A^7 + 59*A^6 + 7*A^5 - 67*A^4 - 21*A^3 - 43*A^2 + 5*A)/(A^7 + A^6 - 19*A^5 - 3*A^4 + 99*A^3 - 61*A^2 - 17*A - 1))*w^2 + ((-9*A^4 - 16*A^3 - 38*A^2 - 40*A + 7)/(A^6 + 2*A^5 -
    17*A^4 - 20*A^3 + 79*A^2 + 18*A + 1)*B + (-5*A^6 - 37*A^5 - 73*A^4 - 18*A^3 + 21*A^2 - 17*A + 1)/(A^6 + 2*A^5 - 17*A^4 - 20*A^3 + 79*A^2 + 18*A + 1))*w + (2*A^7 + 5*A^6 + 16*A^5 +
    16*A^4 - 14*A^3 + 13*A^2 - 4*A - 2)/(A^9 + 3*A^8 - 17*A^7 - 42*A^6 + 91*A^5 + 154*A^4 - 119*A^3 - 114*A^2 - 20*A - 1)*B + (6*A^8 + 43*A^7 + 67*A^6 - 93*A^5 - 164*A^4 + 73*A^3 + 59*A^2 +
    9*A)/(A^9 + 3*A^8 - 17*A^7 - 42*A^6 + 91*A^5 + 154*A^4 - 119*A^3 - 114*A^2 - 20*A - 1);
newr := 1 - (newA + newB)/((newB + 1)*(newA + newB + 1));
news := (1 - newA)/(newB + 1);
newb := newr*news*(newr-1); newc := news*(newr-1);
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;



N := 15;
QA<A> := FunctionField(Rationals());
QAB<B> := PolynomialRing(QA);
modpol := B^2 + (A^2 + A + 1)*B + A^2;
QB<B> := ext<QA | modpol>;
r := 1 + (A*B + B^2) / (A^2*(A + B + 1));
s := 1 + B/(A*(A + 1));
b := r*s*(r-1); c := s*(r-1);
E := EllipticCurve([1-c,-b,-b,0,0]);
P := E ! [0,0];
assert N*P eq E ! 0;
QBz<z> := PolynomialRing(QB);
EA, phi := IsogenyFromKernel(E, &*{(z-(n*P)[1]) : n in [1..N-1]});
t := TatePairing(P,-P,N);
Qw<w> := ext<QB | z^N - t>;
numA := ((A^12 - 2*A^11 - 3*A^10 + 5*A^9 - 7*A^8 - 23*A^7 - 11*A^6 - 11*A^5 - 27*A^4 - 23*A^3 - 8*A^2 - A)*B - 6*A^11 - 9*A^10 - 11*A^9 - 36*A^8 - 50*A^7 - 39*A^6 - 45*A^5 - 51*A^4 - 31*A^3 - 9*A^2
    - A)*w^12 + ((-A^12 + 2*A^11 + 3*A^10 - 9*A^9 + A^8 + 19*A^7 + 2*A^6 - A^5 + 21*A^4 + 22*A^3 + 8*A^2 + A)*B - A^12 + 3*A^11 + 2*A^10 - 4*A^9 + 18*A^8 + 32*A^7 + 17*A^6 + 27*A^5 + 44*A^4
    + 30*A^3 + 9*A^2 + A)*w^11 + ((-2*A^12 - A^11 + 11*A^10 - A^9 - 29*A^8 + 5*A^7 + 31*A^6 - 17*A^5 - 22*A^4 + 28*A^3 + 38*A^2 + 15*A + 2)*B - 2*A^12 + 7*A^10 - 15*A^9 - 29*A^8 + 18*A^7 +
    26*A^6 - 12*A^5 + 19*A^4 + 68*A^3 + 53*A^2 + 17*A + 2)*w^10 + ((2*A^12 - 3*A^11 - 6*A^10 + 4*A^9 - 5*A^8 - 22*A^7 - 22*A^6 - 14*A^5 - 14*A^4 - 25*A^3 - 22*A^2 - 8*A - 1)*B - 10*A^11 -
    14*A^10 - 14*A^9 - 36*A^8 - 53*A^7 - 51*A^6 - 44*A^5 - 46*A^4 - 48*A^3 - 30*A^2 - 9*A - 1)*w^9 + ((-4*A^12 + 4*A^11 + 16*A^10 - 28*A^9 - 39*A^8 + 22*A^7 + 10*A^6 - 41*A^5 - 7*A^4 +
    62*A^3 + 66*A^2 + 27*A + 4)*B - 6*A^12 - 9*A^10 - 61*A^9 - 43*A^8 + 19*A^7 + 3*A^6 - A^5 + 78*A^4 + 132*A^3 + 93*A^2 + 31*A + 4)*w^8 + ((-4*A^11 + 6*A^10 + 14*A^9 - 2*A^8 + 24*A^7 +
    40*A^6 + 15*A^5 + 20*A^4 + 14*A^3 - 4*A^2 - 5*A - 1)*B + (2*A^12 + 4*A^11 + 28*A^10 + 38*A^9 + 38*A^8 + 74*A^7 + 67*A^6 + 34*A^5 + 30*A^4 + 9*A^3 - 9*A^2 - 6*A - 1))*w^7 + ((2*A^10 -
    3*A^9 - 6*A^8 + A^7 - 2*A^6 - 12*A^5 - 8*A^4 + 3*A^3 - 2*A^2 - 2*A)*B - 2*A^11 - 2*A^10 - 9*A^9 - 14*A^8 - 6*A^7 - 14*A^6 - 20*A^5 - 7*A^4 + A^3 - 4*A^2 - 2*A)*w^6 + ((3*A^8 + 5*A^7 +
    5*A^6 + 4*A^5 + 7*A^4 + 3*A^3 - 3*A^2 - 3*A - 2)*B + (2*A^11 + 4*A^10 + 6*A^9 + 14*A^8 + 14*A^7 + 12*A^6 + 11*A^5 + 7*A^4 - A^3 - 6*A^2 - 5*A - 2))*w^5 + ((2*A^9 + 4*A^8 + 8*A^7 - A^6 -
    5*A^5 - A^4 - 7*A^3 - 13*A^2 - 16*A - 2)*B + (2*A^11 + 6*A^10 + 7*A^9 + 9*A^8 + 5*A^7 - 6*A^6 - 22*A^5 - 14*A^4 - 26*A^3 - 30*A^2 - 17*A - 3))*w^4 + ((-2*A^9 - 5*A^8 - 3*A^7 - 11*A^6 -
    2*A^5 - 11*A^4 - 9*A^3 + A^2 - 2*A - 3)*B - 2*A^11 - 7*A^10 - 4*A^9 - 18*A^8 - 14*A^7 - 12*A^6 - 17*A^5 - 12*A^4 - 6*A^3 - 2*A^2 - A - 3)*w^3 + ((-3*A^6 + 7*A^5 + 5*A^4 + 10*A^2 + 2*A +
    9)*B - A^8 + 5*A^6 + A^5 + 17*A^4 + A^3 + 11*A^2 + 11*A + 5)*w^2 + ((2*A^8 + 3*A^7 - 5*A^6 + 10*A^5 - 9*A^4 + 3*A^3 + 5*A^2 - 7*A + 3)*B + (2*A^10 + 5*A^9 - 2*A^8 + 13*A^7 + 3*A^6 -
    3*A^5 + 14*A^4 - 4*A^3 + 3*A^2 - A))*w + (A^8 + A^7 + A^6 + 3*A^5 + A^4 + 2*A^3 + A^2)*B + A^10 + 2*A^9 + 2*A^8 + 5*A^7 + 4*A^6 + 4*A^5 + 4*A^4 + 2*A^3 + A^2;
denA := ((A^12 - 2*A^11 - 3*A^10 + 5*A^9 - 7*A^8 - 23*A^7 - 11*A^6 - 11*A^5 - 27*A^4 - 23*A^3 - 8*A^2 - A)*B - 6*A^11 - 9*A^10 - 11*A^9 - 36*A^8 - 50*A^7 - 39*A^6 - 45*A^5 - 51*A^4 - 31*A^3 - 9*A^2
    - A)*w^12 + ((-A^10 - 12*A^7 - 20*A^6 - 15*A^5 - 22*A^4 - 33*A^3 - 24*A^2 - 8*A - 1)*B - A^11 - 5*A^10 - 10*A^9 - 22*A^8 - 45*A^7 - 56*A^6 - 55*A^5 - 62*A^4 - 58*A^3 - 32*A^2 - 9*A -
    1)*w^11 + ((3*A^12 - 3*A^11 - 10*A^10 + 14*A^9 + 14*A^8 - 14*A^7 + 18*A^6 + 47*A^5 + 23*A^4 + 17*A^3 + 23*A^2 + 12*A + 2)*B + (2*A^12 - 5*A^11 - A^10 + 31*A^9 + 28*A^8 + 24*A^7 + 75*A^6
    + 85*A^5 + 50*A^4 + 42*A^3 + 35*A^2 + 14*A + 2))*w^10 + ((-A^12 - 5*A^11 + 11*A^10 + 8*A^9 - 40*A^8 + 42*A^6 - 24*A^5 - 26*A^4 + 27*A^3 + 25*A^2 + 5*A)*B - 2*A^12 - 6*A^11 + 9*A^10 -
    11*A^9 - 44*A^8 + 19*A^7 + 30*A^6 - 30*A^5 + 6*A^4 + 52*A^3 + 30*A^2 + 5*A)*w^9 + ((2*A^12 + 4*A^11 - 14*A^10 + 12*A^9 + 77*A^8 + 21*A^7 - 3*A^6 + 103*A^5 + 78*A^4 - 22*A^3 - 22*A^2 +
    2*A + 2)*B + (7*A^12 + 18*A^11 + 24*A^10 + 91*A^9 + 149*A^8 + 82*A^7 + 96*A^6 + 161*A^5 + 56*A^4 - 42*A^3 - 20*A^2 + 4*A + 2))*w^8 + ((-2*A^11 + 4*A^10 - 5*A^9 - 18*A^8 + 15*A^7 - 3*A^6
    - 42*A^5 - 7*A^4 - A^3 - 35*A^2 - 28*A - 6)*B - 2*A^12 - 5*A^11 - A^10 - 18*A^9 - 23*A^8 - 43*A^6 - 68*A^5 - 30*A^4 - 42*A^3 - 63*A^2 - 34*A - 6)*w^7 + ((-2*A^10 + 4*A^9 - 5*A^8 -
    17*A^7 + 5*A^6 + 13*A^5 - 3*A^4 + A^3 + 26*A^2 + 23*A + 6)*B - 3*A^11 - 4*A^10 - 4*A^9 - 14*A^8 - 12*A^7 + 16*A^6 + 25*A^5 + 15*A^4 + 33*A^3 + 49*A^2 + 29*A + 6)*w^6 + ((-A^9 + 4*A^8 +
    10*A^7 + 2*A^6 + A^5 + 8*A^4 + 11*A^3 - 4*A - 2)*B + (A^11 + 2*A^10 + 5*A^9 + 12*A^8 + 17*A^7 + 10*A^6 + 10*A^5 + 17*A^4 + 9*A^3 - 4*A^2 - 6*A - 2))*w^5 + ((A^9 - 4*A^8 - A^7 + 5*A^6 +
    10*A^5 - A^4 + 12*A^3 + 16*A^2 + 10*A + 6)*B + (A^11 - A^10 - 3*A^8 + 9*A^7 + 14*A^6 + 19*A^5 + 23*A^4 + 25*A^3 + 29*A^2 + 15*A + 5))*w^4 + ((5*A^8 + 3*A^7 + 7*A^6 + 2*A^5 + 8*A^4 +
    2*A^3 + A^2 - 2*A - 3)*B + (5*A^10 + 6*A^9 + 13*A^8 + 15*A^7 + 16*A^6 + 13*A^5 + 9*A^4 + 8*A^3 - 3*A^2 - 3*A - 2))*w^3 + ((-2*A^8 + 2*A^7 - 2*A^6 + 6*A^5 - 9*A^4 + 12*A^3 - 2*A^2 + 3*A
    + 2)*B - 2*A^10 - 6*A^8 + 5*A^7 - 5*A^6 + 4*A^5 + 4*A^4 + 3*A^3 + 9*A^2 + 3)*w^2 + ((4*A^8 + 2*A^7 + 2*A^6 + 7*A^5 - 2*A^4 - 2*A^3 + 4*A^2 - 7*A - 3)*B + (4*A^10 + 6*A^9 + 4*A^8 +
    13*A^7 + 8*A^6 + 4*A^4 - 3*A^3 - 6*A^2 - 5*A - 5))*w + (-2*A^8 - A^7 - A^6 - 7*A^5 + A^4 - 5*A^3 - 4*A^2 + A - 2)*B - 2*A^10 - 3*A^9 - 2*A^8 - 10*A^7 - 5*A^6 - 6*A^5 - 10*A^4 - 2*A^3 -
    4*A^2 - A;
newA := numA/denA;
numB := ((-A^11 + 2*A^10 + 2*A^9 - 4*A^8 + 6*A^7 + 12*A^6 + 2*A^5 + 5*A^4 + 11*A^3 + 6*A^2 + A)*B + (5*A^10 + 5*A^9 + 5*A^8 + 20*A^7 + 21*A^6 + 12*A^5 + 17*A^4 + 17*A^3 + 7*A^2 + A))*w^9 + ((A^12 -
    A^11 - 3*A^10 + 10*A^9 - 2*A^8 - 16*A^7 + 8*A^6 - A^5 - 22*A^4 - 9*A^3 + 2*A^2 + A)*B + (2*A^12 - A^11 + 8*A^9 - 15*A^8 - 21*A^7 - 2*A^6 - 22*A^5 - 30*A^4 - 7*A^3 + 3*A^2 + A))*w^8 +
    ((-4*A^10 + 3*A^9 + 6*A^8 - 16*A^7 - 14*A^6 - A^5 - 13*A^4 - 21*A^3 - 15*A^2 - 6*A - 1)*B - A^11 - 6*A^10 - A^9 - 11*A^8 - 36*A^7 - 29*A^6 - 25*A^5 - 39*A^4 - 37*A^3 - 21*A^2 - 7*A -
    1)*w^7 + ((3*A^9 - 2*A^8 - 12*A^7 + 2*A^6 + 2*A^5 - 20*A^4 - 12*A^3 + 3*A^2 + 2*A)*B - A^10 - A^9 - 12*A^8 - 20*A^7 - 7*A^6 - 17*A^5 - 30*A^4 - 9*A^3 + 5*A^2 + 2*A)*w^6 + ((A^9 - 3*A^8
    - 5*A^7 + A^5 - 5*A^4 - 3*A^3 + 5*A^2 + 3*A)*B - 2*A^10 - 2*A^9 - 7*A^8 - 8*A^7 - A^6 - 2*A^5 - 5*A^4 + 2*A^3 + 8*A^2 + 3*A)*w^5 + ((-A^9 + 2*A^8 + 2*A^7 + 2*A^6 + 3*A^5 + 2*A^4 + 4*A^3
    + 2*A^2 - A - 1)*B - A^11 + 2*A^9 + 2*A^8 + 7*A^7 + 6*A^6 + 6*A^5 + 6*A^4 + 5*A^3 + A^2 - 2*A - 1)*w^4 + ((2*A^8 + 2*A^6 + 3*A^5 + 2*A^4 + 3*A^3 + A^2 + A)*B + (2*A^10 + 2*A^9 + 5*A^8 +
    7*A^7 + 6*A^6 + 9*A^5 + 7*A^4 + 4*A^3 + 3*A^2 + A))*w^3 + ((-A^7 + 3*A^6 + 4*A^4 + 4*A^3 + 4*A + 1)*B - A^9 + 2*A^8 + A^6 + 7*A^5 + A^4 + 4*A^3 + 3*A^2 + 2*A + 1)*w^2 + ((2*A^6 - 4*A^5
    + 2*A^4 - 4*A^3 - 4*A^2 + 2*A - 4)*B + (2*A^8 - 2*A^7 + 2*A^6 + 2*A^5 - 6*A^4 + 4*A^3 - 2*A^2))*w + (A^6 - A^5 + A^4 - A^2 + A - 1)*B + A^8 + A^6 + 2*A^5 - A^4 + 2*A^3;
denB := ((A^12 - 6*A^10 + 4*A^9 + 12*A^8 - 14*A^7 - 13*A^6 + 12*A^5 - 3*A^4 - 21*A^3 - 12*A^2 - 2*A)*B + (A^12 - A^11 - 5*A^10 + 7*A^9 + 3*A^8 - 25*A^7 - 14*A^6 - A^5 - 26*A^4 - 33*A^3 - 14*A^2 -
    2*A))*w^9 + ((-A^12 + 6*A^10 - 5*A^9 - 7*A^8 + 11*A^7 - 2*A^6 - 14*A^5 - 5*A^4 - 4*A^3 - 4*A^2 - A)*B - A^12 + 2*A^11 + 5*A^10 - 9*A^9 - 5*A^8 + 4*A^7 - 18*A^6 - 22*A^5 - 10*A^4 - 8*A^3
    - 5*A^2 - A)*w^8 + ((-A^11 + 2*A^10 + A^9 - 9*A^8 - 2*A^7 + 4*A^6 - 6*A^5 - 6*A^4 + 3*A^3 + 8*A^2 + 5*A + 1)*B - 2*A^11 - 6*A^9 - 14*A^8 - 3*A^7 - A^6 - 7*A^5 + A^4 + 12*A^3 + 13*A^2 +
    6*A + 1)*w^7 + ((-A^9 + 9*A^7 - 7*A^5 + 13*A^4 + 10*A^3 - 7*A^2 - 6*A - 1)*B + (A^10 + A^9 + 6*A^8 + 14*A^7 + 2*A^6 + 3*A^5 + 18*A^4 + 2*A^3 - 13*A^2 - 7*A - 1))*w^6 + ((-2*A^9 - A^8 +
    3*A^7 + 4*A^6 + 2*A^5 + 5*A^4 + 9*A^3 + 6*A^2 + 3*A + 1)*B - A^11 - A^9 + 3*A^8 + 10*A^7 + 11*A^6 + 12*A^5 + 16*A^4 + 16*A^3 + 9*A^2 + 4*A + 1)*w^5 + ((A^9 + 4*A^7 + 3*A^6 - A^5 + 4*A^4
    + 4*A^3 - A^2 - 2*A)*B + (A^11 + 2*A^10 + 2*A^9 + 7*A^8 + 7*A^7 + 6*A^6 + 4*A^5 + 6*A^4 + 3*A^3 - 3*A^2 - 2*A))*w^4 + ((-2*A^8 + 2*A^7 - 2*A^6 - A^5 + A^2 + 1)*B - 2*A^10 - A^9 - 3*A^8
    - 3*A^7 - A^6 - 4*A^5 + A^4 + A^3 + A^2 + A + 1)*w^3 + ((A^8 - A^6 + 2*A^5 - 6*A^4 + A^3 - 3*A^2 - 3*A - 1)*B + (A^10 + A^9 - 2*A^8 + 3*A^7 - 5*A^6 - 4*A^5 - 4*A^4 - 6*A^3 - 4*A^2 - 4*A
    - 1))*w^2 + ((-2*A^7 - A^6 - 5*A^4 + 2*A^3 - A^2 + 2)*B - 2*A^9 - 3*A^8 - A^7 - 7*A^6 - 2*A^5 - A^4 - 3*A^3 + 2*A^2 + A + 1)*w + (-A^6 - A^4 - 2*A^3 + A^2 - 2*A)*B - A^8 - A^7 - A^6 -
    4*A^5 - A^4 - 3*A^3 - 2*A^2 - A - 1;
newB := numB/denB;
newr := 1 + (newA*newB + newB^2) / (newA^2*(newA + newB + 1));
news := 1 + newB/(newA*(newA + 1));
newb := newr*news*(newr-1); newc := news*(newr-1);
assert same_j_invars(EA, newb, newc);
"Formulas verified for N equal", N;
