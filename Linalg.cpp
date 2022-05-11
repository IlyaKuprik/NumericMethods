#pragma once
#include <bits/stdc++.h>
using namespace std;

long double my_sqrt(long double u, long double eps = 1e-8){
    // Вычисление корня с помощью формулы Герона
    long double xn = 1.0;
    long double d = 1;
    while(d >= eps/3){
        long double xn1 = 0.5*(xn + u/xn);
        d = abs(xn1 - xn);
        xn = xn1;
    }
    return xn;
}

class Vector{
    public:
        vector<long double> vec;
        int size;

        Vector(int size_){
            size = size_;
            vec = vector<long double>(size);
        }

        Vector(vector<long double> arr_v){
            size = arr_v.size();
            vec = vector<long double>(size);
            for(int i  = 0; i < size; ++i){
                vec[i] = arr_v[i];
            }
        }
        
        long double get_norm2(){
            long double s = 0;
            for(int i = 0; i < size; ++i){
                s+=vec[i]*vec[i];
            }
            return my_sqrt(s);
        }

        long double get_norm1(){
            long double s = 0;
            for(int i = 0; i < size; ++i){
                s+=abs(vec[i]);
            }
            return s;
        }

        //умножение вектора на число
        const Vector operator*(const long double & rvalue) const {
            vector<long double> m0(size, 0);
            Vector p(m0);
            for(int i = 0; i < size; ++i){
                p.vec[i] = vec[i] * rvalue;
            }
            return p;
        }

        //сложение векторов
        const Vector operator+(const Vector & rvalue) const {
            assert(size == rvalue.size);
            vector<long double> m0(size, 0);
            Vector p(m0);
            for(int i = 0; i < size; ++i){
                p.vec[i] = vec[i] + rvalue.vec[i];
            }
            return p;
        }
        //разность векторов
        const Vector operator-(const Vector & rvalue) const {
            assert(size == rvalue.size);
            return *this + rvalue*(-1);
        }
        //умножениe вектор на вектор((1 x n) * (n x 1))
        const long double operator*(const Vector & r_vector) const {
            assert(size == r_vector.size);
            long double p = 0;
            for(int i = 0; i < size; ++i){
                p += vec[i]*r_vector.vec[i];
            }
            return p;
        }
        /*//умножение вектора(1 x n) на матрицу (n x n)
        const Vector operator*(const Matrix & r_matrix) const {
            assert(size == r_matrix.size);
            vector<long double> m0(size, 0);
            Vector p(m0);
            for(int j = 0; j < size; ++j){
                long double s = 0;
                for(int k = 0; k < size; ++k){
                    s += vec[k]*r_matrix.matrix[k][j];
                }
                p.vec[j] = s;
            }
            return p;
        }*/

};

class Matrix{
    public:
        vector<vector<long double>> matrix;
        int size;

        Matrix(int size_){
            size = size_;
            matrix = vector<vector<long double>>(size, vector<long double>(size));
        }

        Matrix(vector<vector<long double>> arr_m){
            size = arr_m.size();
            matrix = vector<vector<long double>>(size, vector<long double>(size));
            for(int i = 0;i < size;++i){
                for(int j = 0;j < size; ++j)
                    matrix[i][j] = arr_m[i][j];
            }
        }

        long double get_determinant(){
            if(size == 3)
                return matrix[0][0] * matrix[1][1] * matrix[2][2] + matrix[0][1] * matrix[1][2] * matrix[2][0] +
                matrix[1][0] * matrix[2][1] * matrix[0][2] - 
                matrix[2][0] * matrix[1][1] * matrix[0][2] - matrix[0][1] * matrix[1][0] * matrix[2][2] -
                matrix[0][0] * matrix[1][2] * matrix[2][1];
            return 0;
        }

        Matrix transpose(){
            vector<vector<long double>> m0(size, vector<long double>(size, 0)); // создаем вектор нулей
            Matrix t(m0);
            for(int i = 0; i < size; ++i){
                for (int j = 0; j < size; ++j){
                    t.matrix[i][j] = matrix[j][i];
                }
            }
            return t;
        }

        long double get_norm2(){
            // ||A||2
            long double su = 0;
            for(int j = 0; j < size; ++j){
                for(int i = 0; i < size; ++i){
                    su += matrix[j][i]*matrix[j][i];
                }
            }
            return my_sqrt(su);
        }

        long double get_norm1(){
            // ||A||1
            long double norm = 0;
            for(int j = 0; j < size; ++j){
                long double su = 0;
                for(int i = 0; i < size; ++i){
                    su += abs(matrix[i][j]);
                }
                norm = max(norm, su);
            }
            return norm;
        }

        long double get_norm_inf(){
            // ||A||inf
            long double norm = 0;
            for(int j = 0; j < size; ++j){
                long double su = 0;
                for(int i = 0; i < size; ++i){
                    su += abs(matrix[j][i]);
                }
                norm = max(norm, su);
            }
            return norm;
        }

        void sort(Vector &b){
            //функция, которая меняет строки СЛАУ для метода Гаусса
            for(int o = 0; o < 2; ++o){
                for(int k = 0;k < size;++k){
                    if(matrix[k][k] == 0){
                        for(int i = 0; i < size; ++i){
                            if(i == k)continue;
                            if(matrix[i][k] != 0 && (matrix[k][i] != 0 || o < 1)){
                                long double tv = b.vec[i];
                                b.vec[i] = b.vec[k];
                                b.vec[k] = tv;
                                for(int j = 0; j < size; ++j){
                                    long double t = matrix[i][j];
                                    matrix[i][j] = matrix[k][j];
                                    matrix[k][j] = t;
                                }
                            }
                        }
                    }
                }
            }
        }

        //определение суммы матриц
        const Matrix operator+(const Matrix & r_matrix) const {
            assert(size == r_matrix.size);
            vector<vector<long double>> m0(size, vector<long double>(size, 0)); 
            Matrix s(m0);
            for(int i = 0; i < size; ++i){
                for(int j = 0; j < size; ++j){
                    s.matrix[i][j] += (matrix[i][j] + r_matrix.matrix[i][j]);
                }
            }
            return s;
        }

        const Matrix operator-(const Matrix & r_matrix) const {
            assert(size == r_matrix.size);
            return *this + r_matrix*(-1);
        }
        //определение произведения матрицы
        const Matrix operator*(const Matrix & r_matrix) const {
            assert(size == r_matrix.size);
            vector<vector<long double>> m0(size, vector<long double>(size, 0));
            Matrix p(m0);
            for(int i = 0; i < size; ++i){
                for(int j = 0; j < size; ++j){
                    long double s = 0;
                    for(int k = 0; k < size; ++k){
                        s += matrix[i][k]*r_matrix.matrix[k][j];
                    }
                    p.matrix[i][j] = s;
                }
            }
            return p;
        }

        //умножение матрицы на вектор
        const Vector operator*(const Vector & r_vector) const {
            assert(size == r_vector.size);
            vector<long double> m(size,0);
            Vector p(m);
            for(int i = 0; i < size; ++i){
                for(int j = 0;j < size; ++j){
                    p.vec[i] += matrix[i][j] * r_vector.vec[j];
                }
            }
            return p;
        }
        //умножение матрицы на число
        const Matrix operator*(const long double & rvalue) const {
            vector<vector<long double>> m0(size, vector<long double>(size, 0));
            Matrix p(m0);
            for(int i = 0; i < size; ++i){
                for(int j = 0;j < size; ++j){
                    p.matrix[i][j] = matrix[i][j] * rvalue;
                }
            }
            return p;
        }
};

Vector get_basis_vector(int size, int pos){
    Vector b(size);
    b.vec[pos] = 1;
    return b;
}

ostream& operator << (ostream &out, const Matrix &matrix){
    int size = matrix.size;
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            out << matrix.matrix[i][j] << ' ';
        }
        out << '\n';
    }
    return out;
}

bool is_not_positive(Matrix & m){
    long double d1 = (m.matrix[0][0]);
    long double d2 = (m.matrix[0][0]*m.matrix[1][1] - m.matrix[0][1]*m.matrix[1][0]);
    long double d3 = (m.get_determinant());
    return !(d1 > 0 && d2 > 0 && d3 > 0);
}

ostream& operator << (ostream &out, const Vector &vector){
    int size = vector.size;
    for(int i = 0; i < size; ++i)
        out << vector.vec[i] << ' ';   
    out << '\n';
    return out;
}