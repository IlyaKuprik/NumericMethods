#define _CRT_SECURE_NO_DEPRECATE
#include <bits/stdc++.h>
#include "interpolation.cpp"
#include "spline.cpp"
using namespace std;

long double max_deviation(long double a, long double b, int n, long double func(long double),
                          bool equ_dist = true, string method = "lagrange"){
    long double ma = -1e8;
    int m = n * 100;
    vector<long double> roots;
    //инициализация корней 
    if(equ_dist)
        roots = get_equidistant_roots(a, b, n);
    else
        roots = get_optimal_roots(a, b, n);
    a = roots[0];
    b = roots[roots.size()-1];

    vector<long double> coeffs(roots.size());
    if(method == "newton")
        for(int i = 0; i < roots.size()-1; ++i) coeffs[i] = divided_difference(i+2, roots, func);
    if (method == "slau")
        coeffs = slau_interpolation(roots, func);
    // коэффиценты для сплайнов
    vector<vector<long double>> coefficents;
    if(method == "spline_10")
        coefficents = spline10_coeffs(roots, func);
    if(method == "spline_32")
        coefficents = spline32_coeffs(roots, func);
    if(method == "spline_20")
        coefficents = spline20_coeffs(roots, func);
    if(method == "spline_21")
        coefficents = spline21_coeffs(roots, func);
        
    for(int i = 0; i < m; i++){
        long double x = a + i*(b - a)/(m-1); 
        if(method == "lagrange" || method == "newton" || method == "slau")
            ma = max(ma, abs(func(x) - interpolation_func_at_point(x, roots, func, coeffs, method)));
        else
            ma = max(ma, abs(func(x) - spline_at_point(x, coefficents, roots, func, method)));
    }
    return ma;
}

void print_table(long double a, long double b, long double func(long double), string method){
    if(method == "lagrange")
        cout << "\tn\t\tm\t\tRL_n\t\tRLopt_n\n";
    if(method == "newton")
        cout << "\tn\t\tm\t\tRN_n\t\tRNopt_n\n";
    if(method == "slau")
        cout << "\tn\t\tm\t\tRSLAU_n\t\tRSLAUopt_n\n";
    if(method == "spline_10")
        cout << "\tn\t\tm\t\tRS(1,0)_n\t\tRS(1,0)opt_n\n";
    if(method == "spline_20")
        cout << "\tn\t\tm\t\tRS(2,0)_n\t\tRS(2,0)opt_n\n";
    if(method == "spline_32")
        cout << "\tn\t\tm\t\tRS(3,2)_n\t\tRS(3,2)opt_n\n";
    if(method == "spline_21")
        cout << "\tn\t\tm\t\tRS(2,1)_n\t\tRS(2,1)opt_n\n";

    for(int n_i = 3; n_i < 51; ++n_i){
        cout << setw(15) << n_i;
        cout << setw(15) << n_i*100;
        cout << setw(25) << max_deviation(a, b, n_i, func, true, method);
        cout << setw(25) << max_deviation(a, b, n_i, func, false, method) << '\n';
    }
}

long long factorial(int n){
    long long val = 1; 
    for(int i = 2; i <= n; ++i) val *= i;
    return val;
}

void methodical_err(long double x, long double a, long double b, int n_roots){
    long double M = 1e8;
    vector<long double> eq_roots = get_equidistant_roots(a, b, n_roots);
    vector<long double> op_roots = get_optimal_roots(a, b, n_roots);
    long double root_pol = 1;
    for(int i = 0;i < n_roots; ++i){
        root_pol *= abs(x - eq_roots[i]);
    }
    cout << "\tМетодическая погрешность для равноотстоящих узлов: " << M*root_pol/factorial(n_roots + 1) << endl;
    cout << "\tМетодическая погрешность для оптимальных узлов: " << M*pow(b-a, n_roots + 1)/(factorial(n_roots+1)*pow(2, 2*n_roots + 1)) << endl;
}

long double f(long double x){
    return x*x + 1 - acos(x);
}


void save_grap_data(long double a, long double b, string m_key){
    vector<int> n_roots = {3};
    for(int i = 10; i < 51; i+=10) n_roots.push_back(i);
    for(auto& n_i: n_roots){
        vector<long double> roots = get_optimal_roots(a, b, n_i);
        string path = "data/opt_" + m_key + "_" + to_string(n_i) + ".txt";
        freopen(path.c_str(), "w", stdout);
        if(m_key == "" || m_key == "lagrange" || m_key == "newton")
            interpolation_func_on_interval(1000, a, b, n_i, f, roots, m_key);
        else{
            spline_on_interval(1000, a, b, n_i, f, roots, m_key);
        }
    }
}


int main(){
    long double a = -1+1e-8, b = 1-1e-8;
    //print_table(a, b, f, "lagrange");
    // print_table(a, b, f, "newton");
    // print_table(a, b, f, "slau");
    // print_table(a, b, f, "spline_10");
    // print_table(a, b, f, "spline_21");
    //print_table(a, b, f, "spline_20");
    save_grap_data(a, b, "spline_21");
    return 0;
}