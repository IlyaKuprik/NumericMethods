#define _CRT_SECURE_NO_DEPRECATE
#pragma once
#include <bits/stdc++.h>
using namespace std;

const long double PI = 3.1415926535;

vector<long double> get_equidistant_roots(long double a, long double b, int n){
    vector<long double> roots(n);
    for(int i = 0; i < n; i++){
        roots[i] = a + i*(b - a)/(n-1); 
    }
    return roots;
}

vector<long double> get_optimal_roots(long double a, long double b, int n){
    vector<long double> roots(n+1);
    for(int i = 0; i < n+1; i++){
        roots[i] = 0.5*((b-a)*cos((PI*(2*i + 1.0))/(2*n + 2.0)) + b + a);
    }
    sort(roots.begin(), roots.end());
    return roots;
}
