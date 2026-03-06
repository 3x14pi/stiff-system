#include "funzioni2.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
using namespace std;
// Definizione delle variabili globali
Real g_dt; // Rappresenta il passo temporale corrente h
Real g_y_old[N]; //Contiene il vettore della soluzione al passo temporale precedente
                 //permette di accedere al passo temporale corrente in ogni chiamata alle funzioni che calcolano f
                 //e jacobiano.

// Funzione di fattorizzazione LU
void lu(Real *A, int P[], int n)
{
    Real aux;
    int ia;
    for (int k = 0; k < n; k++)
    {
        P[k] = k;
    }
    for (int k = 0; k < n - 1; k++)
    {
        // Cerco il pivot sulla colonna
        Real pivot = fabs(*(A + k*n + k));
        int ipiv = k;
        for (int kk = k + 1; kk < n; kk++)
        {
            aux = fabs(*(A + kk*n + k));
            if (aux > pivot)
            {
                ipiv = kk;
                pivot = aux;
            }
        }
        // Scambio le righe
        if (ipiv != k)
        {
            for (int kk = 0; kk < n; kk++)
            {
                aux = *(A + k*n + kk);
                *(A + k*n + kk) = *(A + ipiv*n + kk);
                *(A + ipiv*n + kk) = aux;
            }
            ia = P[k];
            P[k] = P[ipiv];
            P[ipiv] = ia;
        }
        for (int i = k + 1; i < n; i++)
        {
            Real m = *(A + i*n + k) / (*(A + k*n + k));
            *(A + i*n + k) = m;
            for (int j = k + 1; j < n; j++)
            {
                *(A + i*n + j) -= m * (*(A + k*n + j));
            }
        }
    }
}

// Funzione che risolve il sistema lineare
void risist(Real *A, int P[], Real x[], Real b[], int n)
{
    Real sum = 0.0;
    for (int k = 0; k < n; k++)
    {
        x[k] = b[P[k]];
    }
    for (int k = 1; k < n; k++)
    {
        sum = 0.0;
        for (int i = 0; i < k; i++)
        {
            sum += *(A + k*n + i) * x[i];
        }
        x[k] = x[k] - sum;
    }
    x[n-1] = x[n-1] / (*(A + (n-1)*n + n-1));
    for (int i = n-2; i >= 0; i--)
    {
        sum = 0.0;
        for (int j = i+1; j < n; j++)
        {
            sum += *(A + i*n + j) * x[j];
        }
        x[i] = (x[i] - sum) / (*(A + i*n + i));
    }
}

//Funzione norma 2
Real norm_2(Real *v, int n) {
    Real sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

//Funzione di Newton
void newton_sist(void (*f)(Real*, Real*), void (*Jf)(Real*, Real*),
                 int n, Real x[], int *nit, Real toll, int nitmax)
{
    Real *J = new Real[n*n];
    Real *F = new Real[n];
    Real *delta = new Real[n];
    int *P = new int[n];
    Jf(J, x);
    lu(J, P, n);
    // Inizio ciclo di Newton
    Real delta1 = 1e200;
    Real delta2 = 0.0;
    *nit = 0;
    do {
        *nit += 1;
        f(F, x);
        for (int i = 0; i < n; i++) {
            F[i] = -F[i];
        }
        risist(J, P, delta, F, n);
        delta2 = norm_2(delta, n);
        if (delta2 > delta1) {
            delete [] J;
            delete [] F;
            delete [] delta;
            delete [] P;
            cout << "delta2 > delta1: Newton non converge" << endl;
            return;
        } else {
            delta1 = delta2;
        }
        for (int i = 0; i < n; i++) {
            x[i] += delta[i];
        }
    } while ((delta1 > toll) && (*nit < nitmax));
    delete [] J;
    delete [] F;
    delete [] delta;
    delete [] P;
    if (*nit == nitmax) {
        cout << "Newton non converge in " << nitmax << " iterazioni" << endl;
        return;
    }
}


// Funzione f: funzione di integrazione
void f_res(Real* F, Real* Y) {
    // Primo componente (i = 0)
    F[0] = Y[0] - g_dt * (stiff * (Y[1] - 2 * Y[0]) - exp(Y[0])) - g_y_old[0];
    // Componenti centrali (i = 1, …, N-2)
    for (int i = 1; i < N - 1; i++) {
        F[i] = Y[i] - g_dt * (stiff * (Y[i+1] - 2 * Y[i] + Y[i-1]) - exp(Y[i])) - g_y_old[i];
    }
    // Ultima componente (i = N-1)
    F[N-1] = Y[N-1] - g_dt * (stiff * (Y[N-2] - 2 * Y[N-1]) - exp(Y[N-1])) - g_y_old[N-1];
}


// Funzione J: calcola la matrice jacobiana del sistema F(Y) rispetto a Y

void J_res(Real* J, Real* Y) {
    // Inizializza tutti gli elementi a 0
    for (int i = 0; i < N * N; i++)
        J[i] = 0.0;
    
    // Riga 0
    J[0 * N + 0] = 1.0 + 2 * g_dt * stiff + g_dt * exp(Y[0]);
    if (N > 1)
        J[0 * N + 1] = - g_dt * stiff;
    
    // Righe centrali
    for (int i = 1; i < N - 1; i++) {
        J[i * N + (i - 1)] = - g_dt * stiff;
        J[i * N + i]       = 1.0 + 2 * g_dt * stiff + g_dt * exp(Y[i]);
        J[i * N + (i + 1)] = - g_dt * stiff;
    }
    
    // Ultima riga
    if (N > 1) {
        J[(N - 1) * N + (N - 2)] = - g_dt * stiff;
        J[(N - 1) * N + (N - 1)] = 1.0 + 2 * g_dt * stiff + g_dt * exp(Y[N - 1]);
    }
}


// Funzione integrazionw: integra il sistema da t = 0 a T utilizzando Eulero Implicito

void integrazione(Real dt, Real T, Real y_out[N]) {
    int nSteps = static_cast<int>(T / dt); // T/dt sono real/double, static_cast tronca la parte decimale
    Real y[N];
    Real pi = 4*atan(1.0);
    // Condizioni iniziali: y_i(0) = sin(0.1*pi*(i+1))
    for (int i = 0; i < N; i++) {
        y[i] = sin(0.1 * pi * (i + 1));
    }
    
    int nit; // numero di iterazioni di Newton per ogni passo
    
    // Ciclo di integrazione nel tempo
    for (int step = 0; step < nSteps; step++) {
        // Aggiorna le variabili globali per il passo corrente:
        // g_dt (passo) e g_y_old (soluzione del passo precedente)
        g_dt = dt;
        for (int i = 0; i < N; i++) {
            g_y_old[i] = y[i];
        }
        // Risolve il sistema non lineare: Y - dt*f(Y) - y_old = 0
        newton_sist(f_res, J_res, N, y, &nit, 1e-10, 20);
    }
    
    // Copia il risultato finale in y_out
    for (int i = 0; i < N; i++) {
        y_out[i] = y[i];
    }
}
