#include <iostream>
#include <cmath>

using namespace std;

double f(double x){
    return x*x - asin(x - 0.2);
}

double df(double x){
    return 2*x - 1/sqrt(1 - (x - 0.2)*(x - 0.2));
}

double ddf(double x){
    return 2 - (x - 0.2)/pow(1 - (x - 0.2)*(x - 0.2), 3/2.);
}

void root_localization(double &a, double &b, int N){
    //локализация методом последовательного перебора
    double h = (b-a)/N;
    double xk = a;
    for(int k = 0;k < N; ++k){
        double xk1 = a + k*h;
        if(f(xk)*f(xk1) < 0){
            a = xk;
            b = xk1;
            return;
        }
        xk = xk1;
    }
    root_localization(a,b,2*N);
}

double newtons_method(double eps = 1e-6){
    double a = -10, b = 10;
    root_localization(a,b,1000);
    double xk = a;
    while(true){
        double xk1;
        if(xk > a && xk < b){
            xk1 = xk;
        }else{
            xk1 = (a+b)/2;
        }
        //критерий остановки
        if(abs(xk - xk1) < eps){
            xk = xk1;
            break;
        }

        double c = f(xk1);
        if(f(a) < 0){
            if(c < 0) a = xk1; 
            else b = xk1;
        }else{
            if(c < 0) b = xk1; 
            else a = xk1;
        }
        xk = xk1;
    }
    return xk;
}

int main(){
    cout << "Решение, полученное с помощью метода Ньютона: " << newtons_method() << endl;
    return 0;
}