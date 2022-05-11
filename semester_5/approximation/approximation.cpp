#define _CRT_SECURE_NO_DEPRECATE
#include <bits/stdc++.h>
#include "least_squares_method.cpp"
#include "C:/vscode/5_semestr/interpolation/roots_utils.cpp"
using namespace std;

long double my_func(long double x){
    return x*x - asin(x - 0.2);
}

long double err_rate(long double value){
    return value + ((double)rand() / RAND_MAX - 0.5)/10; 
}

vector<pair<long double, long double>> print_experimental_dots(vector<long double> roots, long double func(long double)){
    vector<pair<long double, long double>> dots;
    for(int i = 0; i < roots.size(); ++i){
        for(int k = 0; k < 3; ++k){
            long double curr_value = func(roots[i]);
            dots.push_back({roots[i], err_rate(curr_value)});
            cout << roots[i] << ' ' << err_rate(curr_value) << endl;
        }
    }
    return dots;
}

void save_graph_data(long double a, long double b, string method = "normal"){
    vector<int> n_roots = {1, 2, 3, 4, 5};
    freopen("data/noise.txt", "w", stdout);
    vector<long double> roots = get_equidistant_roots(a, b, 100);
    vector<pair<long double, long double>> dots = print_experimental_dots(roots, my_func);

    for(auto& n_i: n_roots){
        string path = "data/" + method + "_" + to_string(n_i) + ".txt";
        freopen(path.c_str(), "w", stdout);
        least_squares_method_on_interval(n_i, dots, method);
    }
}


void print_deviation_table(long double a, long double b){
    vector<long double> roots = get_equidistant_roots(a, b, 100);
    vector<pair<long double, long double>> dots = print_experimental_dots(roots, my_func);
    cout << "\tn\tnormal\tortogonal\n";
    for(int deg = 1; deg < 50; deg++){
        cout << setw(10) << deg;
        cout << setw(15) << least_squares_method_on_interval(deg, dots, "normal", false);
        cout << setw(15) << least_squares_method_on_interval(deg, dots, "ortogonal", false);
        cout << '\n';
    }
}

int main(){
    long double a = -0.8 + 1e-8, b = 1.2 - 1e-8;
    print_deviation_table(a, b);;
    //save_graph_data(a, b, "normal");
    //save_graph_data(a, b, "ortogonal");
    return 0;
}