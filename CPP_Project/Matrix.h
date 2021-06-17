#pragma once
#ifndef PROJECT_MATRIX_MATRIX_H
#define PROJECT_MATRIX_MATRIX_H
#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
using namespace std;
template<typename T>
class Matrix {
private:
    vector<vector<T>> matrix;
    int row;
    int col;
private:
    T get_det(vector<vector<T>> m, int n);
    void QR(Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R);
    bool check_df(const string& T_str);
public:
    //constructor
    Matrix(int row = 0, int cow = 0) {
        //initialize
        if (row > 0 && cow > 0) {
            this->row = row;
            this->col = cow;
            for (int i = 0; i < row; ++i) {
                vector<T> r;
                for (int j = 0; j < col; ++j) {
                    r.emplace_back(0);
                }
                matrix.emplace_back(r);
            }
        }
        else {
            this->row = 0;
            this->col = 0;
        }
    };
    explicit Matrix(vector<vector<T>> m) {
        if (m.size() > 0) {
            if (m[0].size() > 0) {
                this->matrix = m;
                row = m.size();
                col = m[0].size();
            }
            else {
                row = 0;
                col = 0;
            }
        }
        else {
            row = 0;
            col = 0;
        }
    }
    Matrix(const Matrix<T>& m) {
        matrix = m.matrix;
        row = m.row;
        col = m.col;
    }
    //getter & setter
    int getRow() { return row; }
    int getCol() { return col; }
    vector<vector<T>> getMatrix() { return matrix; }
    void setMatrix(vector<vector<T>> m) {
        if (m.size() > 0) {
            if (m[0].size() > 0) {
                this->matrix = m;
                row = m.size();
                col = m[0].size();
            }
            else {
                row = 0;
                col = 0;
            }
        }
        else {
            row = 0;
            col = 0;
        }
    }
    bool un_initialized();
    //Q3
    template<typename S>
    friend ostream& operator<<(ostream& os, const Matrix<S>& other);
    template<typename S>
    friend ostream& operator<<(ostream& os, const vector<S>& other);
    Matrix<T>& operator=(const Matrix<T>& next);
    Matrix<T> operator+(const Matrix<T>& next);
    template<typename S>
    friend vector<S> operator+(const vector<S>& v1, const vector<S>& v2);
    Matrix<T> operator-(const Matrix<T>& next);
    template<typename S>
    friend vector<S> operator-(vector<S>& v1, vector<S>& v2);
    Matrix<T> operator*(T val);
    template<typename S>
    friend Matrix<S> operator*(S val, const Matrix<S>& other);
    template<typename S>
    friend vector<S> operator*(vector<S>& v, S val);
    template<typename S>
    friend vector<S> operator*(S val, vector<S>& v);
    Matrix<T> operator/(T val);
    template<typename S>
    friend vector<S> operator/(vector<S>& v, S val);
    static Matrix<T> transposition(const Matrix<T>& m);
    Matrix<T> transposition();
    static Matrix<T> conjugation(const Matrix<T>& m);
    Matrix<T> conjugation();
    Matrix<T> element_wise_mul(const Matrix<T>& m);
    static Matrix<T> element_wise_mul(const Matrix<T>& m1, const Matrix<T>& m2);
    static vector<T> element_wise_mul(vector<T>& v1, vector<T>& v2);
    Matrix<T> operator*(const Matrix<T>& m);
    vector<T> operator*(vector<T>& vec);
    template<typename S>
    friend vector<S> operator*(vector<S>& vec, const Matrix<S>& other);
    template<typename S>
    friend S operator*(vector<S> v1, vector<S> v2);
    static T dot_product(Matrix<T>& m1, Matrix<T>& m2);
    T dot_product(Matrix<T>& m);
    static T dot_product(vector<T>& v1, vector<T>& v2);
    static vector<T> cross_product(vector<T>& v1, vector<T>& v2);
    //Q4
    bool check_type_axis(int type, int axis);
    T find_max();
    T find_max(int type, int axis);
    T find_min();
    T find_min(int type, int axis);
    T sum();
    T sum(int type, int axis);
    T avg();
    T avg(int type, int axis);
    //Q5
    bool check_square();
    T trace();
    T determinant();
    Matrix<T> adjugate();
    bool has_inverse();
    Matrix<T> inverse();
    T norm2();
    static T norm2(Matrix<T>& m);
    static T norm2(vector<T>& v);
    vector<T> eigenvalues();
    vector<vector<T>> eigenvectors();
    void print_eigenvectors();
    //Q6
    static Matrix<T> reshape(Matrix<T> m, int r, int c);
    Matrix<T> reshape(int r, int c);
    static Matrix<T> slicing(Matrix<T> m, int row_begin, int row_end, int col_begin, int col_end, int row_interval, int col_interval);
    static Matrix<T> row_slicing(Matrix<T> m, int row_begin, int row_end, int interval);
    static Matrix<T> col_slicing(Matrix<T> m, int col_begin, int col_end, int interval);
    Matrix<T> slicing(int row_begin, int row_end, int col_begin, int col_end, int row_interval, int col_interval);
    Matrix<T> row_slicing(int row_begin, int row_end, int interval);
    Matrix<T> col_slicing(int col_begin, int col_end, int interval);
    //Q7
    static Matrix<T> convolution(Matrix<T> m1, Matrix<T> m2);
};

//解决复数编译报错问题, 如果泛型类型是complex会自动调用该库
template<typename D>
class Matrix<complex<D>> {
private:
    int row;
    int col;
    vector<vector<complex<D>>> matrix;
private:
    complex<D> get_det(vector<vector<complex<D>>> m, int n);
public:
    //constructor
    Matrix(int row = 0, int cow = 0) {
        //initialize
        if (row > 0 && cow > 0) {
            this->row = row;
            this->col = cow;
            for (int i = 0; i < row; ++i) {
                vector<complex<D>> r;
                for (int j = 0; j < col; ++j) {
                    r.emplace_back(0);
                }
                matrix.emplace_back(r);
            }
        }
        else {
            this->row = 0;
            this->col = 0;
        }
    };
    explicit Matrix(vector<vector<complex<D>>> m) {
        if (m.size() > 0) {
            if (m[0].size() > 0) {
                this->matrix = m;
                row = m.size();
                col = m[0].size();
            }
            else {
                row = 0;
                col = 0;
            }
        }
        else {
            row = 0;
            col = 0;
        }
    }
    Matrix(const Matrix<complex<D>>& m) {
        matrix = m.matrix;
        row = m.row;
        col = m.col;
    }
    //getter & setter
    int getRow() { return row; }
    int getCol() { return col; }
    vector<vector<complex<D>>> getMatrix() { return matrix; }
    void setMatrix(vector<vector<complex<D>>> m) {
        if (m.size() > 0) {
            if (m[0].size() > 0) {
                this->matrix = m;
                row = m.size();
                col = m[0].size();
            }
            else {
                row = 0;
                col = 0;
            }
        }
        else {
            row = 0;
            col = 0;
        }
    }
    bool un_initialized();
    //Q3
    template<typename S>
    friend ostream& operator<<(ostream& os, const Matrix<complex<S>>& other);
    template<typename S>
    friend ostream& operator<<(ostream& os, const vector<complex<S>>& other);
    Matrix<complex<D>>& operator=(const Matrix<complex<D>>& next);
    Matrix<complex<D>> operator+(const Matrix<complex<D>>& next);
    template<typename S>
    friend vector<complex<S>> operator+(const vector<complex<S>>& v1, const vector<complex<S>>& v2);
    Matrix<complex<D>> operator-(const Matrix<complex<D>>& next);
    template<typename S>
    friend vector<complex<S>> operator-(vector<complex<S>>& v1, vector<complex<S>>& v2);
    Matrix<complex<D>> operator*(complex<D> val);
    template<typename S>
    friend Matrix<complex<S>> operator*(complex<S> val, const Matrix<complex<S>>& other);
    template<typename S>
    friend vector<complex<S>> operator*(vector<complex<S>>& v, complex<S> val);
    template<typename S>
    friend vector<complex<S>> operator*(complex<S> val, vector<complex<S>>& v);
    Matrix<complex<D>> operator/(complex<D> val);
    template<typename S>
    friend vector<complex<S>> operator/(vector<complex<S>>& v, complex<S> val);
    static Matrix<complex<D>> transposition(const Matrix<complex<D>>& m);
    Matrix<complex<D>> transposition();
    static Matrix<complex<D>> conjugation(const Matrix<complex<D>>& m);
    Matrix<complex<D>> conjugation();
    Matrix<complex<D>> element_wise_mul(const Matrix<complex<D>>& m);
    static Matrix<complex<D>> element_wise_mul(const Matrix<complex<D>>& m1, const Matrix<complex<D>>& m2);
    static vector<complex<D>> element_wise_mul(vector<complex<D>>& v1, vector<complex<D>>& v2);
    Matrix<complex<D>> operator*(const Matrix<complex<D>>& m);
    vector<complex<D>> operator*(vector<complex<D>>& vec);
    template<typename S>
    friend vector<complex<S>> operator*(vector<complex<S>>& vec, const Matrix<complex<S>>& other);
    template<typename S>
    friend complex<S> operator*(vector<complex<S>> v1, vector<complex<S>> v2);
    static complex<D> dot_product(Matrix<complex<D>>& m1, Matrix<complex<D>>& m2);
    complex<D> dot_product(Matrix<complex<D>>& m);
    static complex<D> dot_product(vector<complex<D>>& v1, vector<complex<D>>& v2);
    static vector<complex<D>> cross_product(vector<complex<D>>& v1, vector<complex<D>>& v2);
    //Q4
    bool check_type_axis(int type, int axis);
    complex<D> sum();
    complex<D> sum(int type, int axis);
    complex<D> avg();
    complex<D> avg(int type, int axis);
    //Q5
    bool check_square();
    complex<D> trace();
    complex<D> determinant();
    Matrix<complex<D>> adjugate();
    bool has_inverse();
    Matrix<complex<D>> inverse();
    D norm2();
    static D norm2(Matrix<complex<D>>& m);
    static D norm2(vector<complex<D>>& v);
    //Q6
    static Matrix<complex<D>> reshape(Matrix<complex<D>> m, int r, int c);
    Matrix<complex<D>> reshape(int r, int c);
    static Matrix<complex<D>> slicing(Matrix<complex<D>> m, int row_begin, int row_end, int col_begin, int col_end, int row_interval, int col_interval);
    static Matrix<complex<D>> row_slicing(Matrix<complex<D>> m, int row_begin, int row_end, int interval);
    static Matrix<complex<D>> col_slicing(Matrix<complex<D>> m, int col_begin, int col_end, int interval);
    Matrix<complex<D>> slicing(int row_begin, int row_end, int col_begin, int col_end, int row_interval, int col_interval);
    Matrix<complex<D>> row_slicing(int row_begin, int row_end, int interval);
    Matrix<complex<D>> col_slicing(int col_begin, int col_end, int interval);
    //Q7
    static Matrix<complex<D>> convolution(Matrix<complex<D>> m1, Matrix<complex<D>> m2);
};
#endif //PROJECT_MATRIX_MATRIX_H
