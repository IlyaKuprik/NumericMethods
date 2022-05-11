#define _CRT_SECURE_NO_DEPRECATE
#pragma once
#include "roots.cpp"
#include "../../Linalg.cpp"
#include "../../algs.cpp"
#include <bits/stdc++.h>
using namespace std;

long double lagrange_k(long double x, int k, const vector<long double> & roots){
    long double lk = 1;
    for(int i = 0; i < roots.size(); ++i){
        if(i == k) continue;
        lk *= (x - roots[i])/(roots[k] - roots[i]);
    }
    return lk;
}

long double lagrange_interpolation(long double x, const vector<long double> & roots, long double func(long double)){
    long double L = 0;
    for(int i = 0; i < roots.size(); ++i){
        L += func(roots[i])*lagrange_k(x, i, roots);
    }
    return L;
}

//разделенные разности для Ньютона
long double divided_g(int i, int k, const vector<long double> & roots){
    long double g_i_k = 1;
    for(int j = 0; j < k; ++j){
        if(i == j) continue;
        g_i_k /= (roots[i] - roots[j]);
    }
    return g_i_k;
}

long double divided_difference(int k, const vector<long double> & roots, long double func(long double)){
    vector<long double> g(k, 1);
    long double f_k = 0;
    for(int i = 0;i < k; ++i){
        g[i] = divided_g(i, k, roots);
    }
    for(int i = 0; i < k; ++i)
        f_k += g[i]*func(roots[i]);
    return f_k;
}

long double newthon_interpolation(long double x, const vector<long double> & roots, long double func(long double),
                                  vector<long double> d_v){
    long double N = func(roots[0]);
    long double w = 1;

    for(int k = 1; k < roots.size(); ++k){
        w*= (x-roots[k-1]);
        N += w*d_v[k-1];
    }
    return N;
}

//штрафное 1
vector<long double> slau_interpolation(const vector<long double> & roots, long double func(long double)){
    
    Matrix A(roots.size() + 1);
    Vector b(roots.size() + 1);
    for(int i = 0; i < roots.size()+1; ++i){
        b.vec[i] = func(roots[i]);
        for(int j = 0;j < roots.size() +1; ++j){
            A.matrix[i][j] = pow(roots[i], roots.size()-j);
        }
    }

    vector<long double> coeffs = reflection_method(A, b).vec;
    return coeffs;
}

long double slau_interpolation_at_point(long double x, const vector<long double> & roots, long double func(long double), vector<long double> coeffs){
    long double value = 0;
    for(int i = 0;i < roots.size()+1; ++i){
        value += coeffs[i]*pow(x, roots.size()-i);
    }

    return value;
}


long double interpolation_func_at_point(long double x, vector<long double> roots, long double func(long double),
                                        vector<long double> coeffs, string method = "lagrange"){
    long double func_value;
    if(method == "lagrange")
        func_value = lagrange_interpolation(x, roots, func);
    else if(method == "newton"){
        func_value = newthon_interpolation(x, roots, func, coeffs);
    }
    else if(method == "slau"){
        func_value = slau_interpolation_at_point(x, roots, func, coeffs);
    }
    return func_value;
}

pair<vector<long double>, vector<long double>> interpolation_func_on_interval(int n_dots, long double a, long double b,
                            int n_roots, long double func(long double),vector<long double> roots, string method){
   
    vector<long double> dots(n_dots);
    vector<long double> values(n_dots);

    //инициализация вектора коэфицентов
    vector<long double> coeffs(roots.size());
    if(method == "newton")
        for(int i = 0; i < roots.size()-1; ++i) coeffs[i] = divided_difference(i+2, roots, func);
    if(method == "slau")
        coeffs = slau_interpolation(roots, func);

    for(int i = 0; i < n_dots; ++i){
        long double x = a + i*((b-a)/(n_dots - 1));
        dots[i] = x;
        values[i] = interpolation_func_at_point(x, roots, func, coeffs, method);
        cout << x << " " << values[i] << '\n';
    }
    return make_pair(dots, values);
}


