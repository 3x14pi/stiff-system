#ifndef FUNZIONI2_H
#define FUNZIONI2_H

#include <cmath>
#include <cstdlib>

typedef double Real;

// Dimensione del sistema e parametro stiff
const int N = 9;
const Real stiff = 1e7;

// Variabili globali per Eulero Implicito:
extern Real g_dt;           // g_dt: il passo temporale corrente
extern Real g_y_old[N];     // g_y_old: vettore che contiene la soluzione al passo precedente

// Dichiarazioni delle funzioni

// LU decomposizione
void lu(Real *A, int P[], int n);

// risist: risolve il sistema lineare
void risist(Real *A, int P[], Real x[], Real b[], int n);

// Calcola la norma 2
Real norm_2(Real *v, int n);

// Funzione newton_sist
void newton_sist(void (*f)(Real*, Real*), void (*Jf)(Real*, Real*),
                 int n, Real x[], int *nit, Real toll, int nitmax);

// f_res: calcola il vettore F(Y) = Y - g_dt * f(Y) - g_y_old
void f_res(Real* F, Real* Y);

// J: calcola la matrice jacobiana del sistema F(Y) rispetto a Y
void J_res(Real* J, Real* Y);

// integrate: integra il sistema da t = 0 a T usando il metodo di Eulero
// e restituisce la soluzione finale in y_out.
void integrazione(Real dt, Real T, Real y_out[N]);

#endif
