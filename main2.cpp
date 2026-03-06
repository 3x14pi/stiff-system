#include <iostream>
#include <iomanip>
#include <cmath>
#include "funzioni2.h"
using namespace std;

int main() {
    const Real T = 100.0;  // tempo finale
    const Real dt1 = 0.1;  // passo temporale "meno fine"
    const Real dt2 = 0.05; // passo temporale "più fine"
    
    Real y_dt1[N], y_dt2[N];
    
    // Integrazione con dt1 e dt2
    integrazione(dt1, T, y_dt1);
    integrazione(dt2, T, y_dt2);
    
    // Valutazione dell'errore: norma 2 della differenza
    Real diff[N];
    for (int i = 0; i < N; i++) {
        diff[i] = y_dt1[i] - y_dt2[i];
    }
    Real errNorm2 = norm_2(diff, N);

    
    // Visualizzazione dei risultati
    cout << fixed << setprecision(10);
    cout << "Errore stimato (norma2 fra dt = " << dt1 << " e dt = " << dt2 << "): " << errNorm2 << "\n\n";
    
    cout << "Punti (i/10, y_i(100)) (utilizzando la soluzione con dt = " << dt1 << "):" << endl;
    for (int i = 0; i < N; i++) {
        cout << (i + 1) / 10.0 << " " << y_dt1[i] << endl;
    }
    cout << "Punti (i/10, y_i(100)) (utilizzando la soluzione con dt = " << dt2 << "):" << endl;
    for (int i = 0; i < N; i++) {
        cout << (i + 1) / 10.0 << " " << y_dt2[i] << endl;
    }
    
    return 0;
}
