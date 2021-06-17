#include "Matrix.h"
//检查T是不是浮点数
template<typename T>
bool Matrix<T>::check_df(const string& T_str) {
    string s1 = typeid(double).name();
    string s2 = typeid(double_t).name();
    string s3 = typeid(float).name();
    string s4 = typeid(float_t).name();
    if ((T_str == s1) || (T_str == s2) || (T_str == s3) || (T_str == s4)) return true;
    return false;
}
//检查初始化
template<typename T>
bool Matrix<T>::un_initialized() {
    if (row == 0 || col == 0)return true;
    else return false;
}
//overload cout<<
template<typename T>
ostream& operator<<(ostream& os, const Matrix<T>& other) {
    try {
        if (other.row == 0 || other.col == 0) {
            throw 1;
        }
        else {
            for (int i = 0; i < other.row; ++i) {
                for (int j = 0; j < other.col; ++j) {
                    os << other.matrix[i][j] << " ";
                }
                os << endl;
            }
        }
    }
    catch (int i)
    {
        os << "The matrix isn't initialized." << endl;
    }
    return os;
}
template<typename T>
ostream& operator<<(ostream& os, const vector<T>& other) {
    try
    {
        if (other.size() == 0) {
            throw 0;
        }
        else {
            os << "[";
            for (int i = 0; i < other.size() - 1; ++i) {
                os << other[i] << ",";
            }
            os << other[other.size() - 1] << "]"<<endl;
        }
    }
    catch (int i)
    {
        os << "The vector isn't initialized." << endl;
    }
    return os;
}
//overload =
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& next) {
    if (this != &next) {
        this->row = next.row;
        this->col = next.col;
        this->matrix = next.matrix;
    }
    return *this;
}
//A+B
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& next) {
    try {
        if (un_initialized() || next.row == 0 || next.col == 0) {
            throw 0;
        }
        vector<vector<T>> add;
        //矩阵大小匹配
        if (this->row == next.row && this->col == next.col) {
            for (int i = 0; i < this->row; ++i) {
                vector<T> r;
                for (int j = 0; j < this->col; ++j) {
                    r.emplace_back(this->matrix[i][j] + next.matrix[i][j]);
                }
                add.emplace_back(r);
            }
            return Matrix<T>(add);
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        //矩阵大小不匹配
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
//v1+v2
template<typename T>
vector<T> operator+(const vector<T>& v1, const vector<T>& v2) {
    vector<T> v;
    try {
        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); ++i) {
                v.emplace_back(v1[i] + v2[i]);
            }
            return v;
        }
        throw 0;
    }
    catch (int i)
    {
        cout << "The sizes of vectors are mismatched." << endl;
    }
}
//A-B
template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& next) {
    try {
        if (un_initialized() || next.row == 0 || next.col == 0) {
            throw 0;
        }
        vector<vector<T>> sub;
        //矩阵大小匹配
        if (this->row == next.row && this->col == next.col) {
            for (int i = 0; i < this->row; ++i) {
                vector<T> r;
                for (int j = 0; j < this->col; ++j) {
                    r.emplace_back(this->matrix[i][j] - next.matrix[i][j]);
                }
                sub.emplace_back(r);
            }
            return Matrix<T>(sub);
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        //矩阵大小不匹配
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
//v1-v2
template<typename T>
vector<T> operator-(vector<T>& v1, vector<T>& v2) {
    try {
        if (v1.size() == 0 || v2.size() == 0) {
            throw 0;
        }
        vector<T> v;
        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); ++i) {
                v.emplace_back(v1[i] - v2[i]);
            }
            return v;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of vectors are mismatched." << endl;
    }
}
//scalar multiplication
template<typename T>
Matrix<T> Matrix<T>::operator*(T val) {
    try {
        if (un_initialized()) {
            throw 0;
        }
        vector<vector<T>> mul;
        for (int i = 0; i < this->row; ++i) {
            vector<T> r;
            for (int j = 0; j < this->col; ++j) {
                r.emplace_back(this->matrix[i][j] * val);
            }
            mul.emplace_back(r);
        }
        return Matrix<T>(mul);
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
}
template<typename T>
Matrix<T> operator*(T val, const Matrix<T>& other) {
    try {
        if (other.row == 0 || other.col == 0) {
            throw 0;
        }
        vector<vector<T>> mul;
        for (int i = 0; i < other.row; ++i) {
            vector<T> r;
            for (int j = 0; j < other.col; ++j) {
                r.emplace_back(val * other.matrix[i][j]);
            }
            mul.emplace_back(r);
        }
        return Matrix<T>(mul);
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }

}
template<typename T>
vector<T> operator*(vector<T>& v, T val) {
    try {
        if (v.size() == 0) {
            throw 0;
        }
        vector<T> mul;
        for (int i = 0; i < v.size(); ++i) {
            mul.emplace_back(v[i] * val);
        }
        return mul;
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }

}
template<typename T>
vector<T> operator*(T val, vector<T>& v) {
    try {
        if (v.size() == 0) {
            throw 0;
        }
        vector<T> mul;
        for (int i = 0; i < v.size(); ++i) {
            mul.emplace_back(v[i] * val);
        }
        return mul;
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }

}
//scalar division
template<typename T>
Matrix<T> Matrix<T>::operator/(T val) {
    try {
        if (un_initialized()) {
            throw 0;
        }
        vector<vector<T>> div;
        for (int i = 0; i < this->row; ++i) {
            vector<T> r;
            for (int j = 0; j < this->col; ++j) {
                r.emplace_back(this->matrix[i][j] / val);
            }
            div.emplace_back(r);
        }
        return Matrix<T>(div);
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }

}
template<typename T>
vector<T> operator/(vector<T>& v, T val) {
    try {
        if (v.size() == 0) {
            throw 0;
        }
        vector<T> div;
        for (int i = 0; i < v.size(); ++i) {
            div.emplace_back(v[i] / val);
        }
        return div;
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }

}
//transposition
template<typename T>
Matrix<T> Matrix<T>::transposition(const Matrix<T>& m) {
    try {
        if (m.row == 0 || m.col == 0) {
            throw 0;
        }
        vector<vector<T>> trans;
        for (int i = 0; i < m.col; ++i) {
            vector<T> trans_r;
            for (int j = 0; j < m.row; ++j) {
                trans_r.emplace_back(m.matrix[j][i]);
            }
            trans.emplace_back(trans_r);
        }
        return Matrix<T>(trans);
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
}
template<typename T>
Matrix<T> Matrix<T>::transposition() {
    try {
        if (un_initialized()) {
            throw 0;
        }
        vector<vector<T>> trans;
        for (int i = 0; i < this->col; ++i) {
            vector<T> trans_r;
            for (int j = 0; j < this->row; ++j) {
                trans_r.emplace_back(this->matrix[j][i]);
            }
            trans.emplace_back(trans_r);
        }
        return Matrix<T>(trans);
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }

};
//conjugation 共轭
template<typename T>
Matrix<T> Matrix<T>::conjugation(const Matrix<T>& m) {
    return transposition(m);
}
template<typename T>
Matrix<T> Matrix<T>::conjugation() {
    return transposition();
}
//element-wise multiplication
template<typename T>
Matrix<T> Matrix<T>::element_wise_mul(const Matrix<T>& m) {
    try {
        if (m.row == 0 || m.col == 0 || row == 0 || col == 0) {
            throw 0;
        }
        vector<vector<T>> mul;
        if (row == m.row && col == m.col) {
            for (int i = 0; i < row; ++i) {
                vector<T> r;
                for (int j = 0; j < col; ++j) {
                    r.emplace_back(matrix[i][j] * m.matrix[i][j]);
                }
                mul.emplace_back(r);
            }
            return Matrix<T>(mul);
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
template<typename T>
Matrix<T> Matrix<T>::element_wise_mul(const Matrix<T>& m1, const Matrix<T>& m2) {
    try {
        if (m1.row == 0 || m1.col == 0 || m2.row == 0 || m2.col == 0) {
            throw 0;
        }
        vector<vector<T>> mul;
        if (m1.row == m2.row && m1.col == m2.col) {
            for (int i = 0; i < m1.row; ++i) {
                vector<T> r;
                for (int j = 0; j < m1.col; ++j) {
                    r.emplace_back(m1.matrix[i][j] * m2.matrix[i][j]);
                }
                mul.emplace_back(r);
            }
            return Matrix<T>(mul);
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
template<typename T>
vector<T> Matrix<T>::element_wise_mul(vector<T>& v1, vector<T>& v2) {
    try {
        if (v1.size() == 0 || v2.size() == 0) {
            throw 0;
        }
        vector<T> mul;
        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); ++i) {
                mul.emplace_back(v1[i] * v2[i]);
            }
            return mul;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of vectors are mismatched." << endl;
    }
}
//matrix-matrix multiplication
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& m) {
    try {
        if (row == 0 || col == 0 || m.row == 0 || m.col == 0) {
            throw 0;
        }
        vector<vector<T>> mul;
        if (col == m.row) {
            //initialize
            for (int i = 0; i < row; ++i) {
                vector<T> r;
                for (int j = 0; j < m.col; ++j) {
                    r.emplace_back(0);
                }
                mul.emplace_back(r);
            }
            for (int i = 0; i < row; ++i) {
                for (int j = 0; j < m.col; ++j) {
                    for (int k = 0; k < col; ++k) {
                        mul[i][j] += matrix[i][k] * m.matrix[k][j];
                    }
                }
            }
            return Matrix<T>(mul);
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
//matrix-vector multiplication
template<typename T>
vector<T> Matrix<T>::operator*(vector<T>& vec) {
    try {
        if (row == 0 || col == 0 || vec.size() == 0) {
            throw 0;
        }
        vector<T> v;
        if (col == vec.size()) {
            //initialize
            for (int i = 0; i < row; ++i) {
                v.emplace_back(0);
            }
            for (int i = 0; i < row; ++i) {
                for (int j = 0; j < col; ++j) {
                    v[i] += matrix[i][j] * vec[j];
                }
            }
            return v;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix or the vector isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The size of matrix and the size of vector are mismatched." << endl;
    }
}
template<typename T>
vector<T> operator*(vector<T>& vec, const Matrix<T>& other) {
    try {
        if (other.row == 0 || other.col == 0 || vec.size() == 0) {
            throw 0;
        }
        vector<T> v;
        if (vec.size() == other.row) {
            //initialize
            for (int i = 0; i < other.col; ++i) {
                v.emplace_back(0);
            }
            for (int i = 0; i < other.col; ++i) {
                for (int j = 0; j < other.row; ++j) {
                    v[i] += other.matrix[j][i] * vec[j];
                }
            }
            return v;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix or the vector isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The size of matrix and the size of vector are mismatched." << endl;
    }
}
template<typename T>
T operator*(vector<T> v1, vector<T> v2) {
    try {
        if (v1.size() == 0 || v2.size() == 0) {
            throw 0;
        }
        T mul = 0;
        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); ++i) {
                mul += v1[i] * v2[i];
            }
            return mul;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of vectors are mismatched." << endl;
    }
}
//dot product 内积
template<typename T>
T Matrix<T>::dot_product(Matrix<T>& m1, Matrix<T>& m2) {
    try {
        if (m1.un_initialized() || m2.un_initialized()) {
            throw 0;
        }
        T dot_prod = 0;
        if (m1.row == m2.row && m1.col == m2.col) {
            for (int i = 0; i < m1.row; ++i) {
                for (int j = 0; j < m1.col; ++j) {
                    dot_prod += m1.matrix[i][j] * m2.matrix[i][j];
                }
            }
            return dot_prod;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
template<typename T>
T Matrix<T>::dot_product(Matrix<T>& m) {
    try {
        if (m.un_initialized() || un_initialized()) {
            throw 0;
        }
        T dot_prod = 0;
        if (row == m.row && col == m.col) {
            for (int i = 0; i < row; ++i) {
                for (int j = 0; j < col; ++j) {
                    dot_prod += matrix[i][j] * m.matrix[i][j];
                }
            }
            return dot_prod;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
template<typename T>
T Matrix<T>::dot_product(vector<T>& v1, vector<T>& v2) {
    return v1 * v2;
}
//cross product n=3
template<typename T>
vector<T> Matrix<T>::cross_product(vector<T>& v1, vector<T>& v2) {
    try {
        vector<T> prod;
        if (v1.size() == 3 && v2.size() == 3) {
            prod.emplace_back(v1[1] * v2[2] - v1[2] * v2[1]);
            prod.emplace_back(v1[2] * v2[0] - v1[0] * v2[2]);
            prod.emplace_back(v1[0] * v2[1] - v1[1] * v2[0]);
            return prod;
        }
        throw 0;
    }
    catch (int i)
    {
        cout << "The sizes of vectors are not 3 and have no cross product." << endl;
    }
}
//check type and axis
template<typename T>
bool Matrix<T>::check_type_axis(int type, int axis) {
    if (un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return false;
    }
    else {
        if (type == 0) {//row
            if (axis >= 0 && axis < row) return true;
            else {
                cout << "Invalid axis." << endl;
                return false;
            }
        }
        else if (type == 1) {//col
            if (axis >= 0 && axis < col) return true;
            else {
                cout << "Invalid axis." << endl;
                return false;
            }
        }
        else {
            cout << "Invalid type." << endl;
            return false;
        }
    }
}
//find the maximum value
template<typename T>
T Matrix<T>::find_max() {
    if (un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return 0;
    }
    T max = matrix[0][0];
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            if (max < matrix[i][j]) {
                max = matrix[i][j];
            }
        }
    }
    return max;
}
template<typename T>
T Matrix<T>::find_max(int type, int axis) {
    if (check_type_axis(type, axis)) {
        T max = 0;
        if (type == 0) {//row
            int r = axis;
            max = matrix[r][0];
            for (int i = 0; i < col; ++i) {
                if (max < matrix[r][i]) {
                    max = matrix[r][i];
                }
            }
        }
        else {//type==1, col
            int c = axis;
            max = matrix[0][c];
            for (int i = 0; i < row; ++i) {
                if (max < matrix[i][c]) {
                    max = matrix[i][c];
                }
            }
        }
        return max;
    }
    else return 0;
}
//find the minimum value
template<typename T>
T Matrix<T>::find_min() {
    if (un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return 0;
    }
    T min = matrix[0][0];
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            if (min > matrix[i][j]) {
                min = matrix[i][j];
            }
        }
    }
    return min;
}
template<typename T>
T Matrix<T>::find_min(int type, int axis) {
    if (check_type_axis(type, axis)) {
        T min = 0;
        if (type == 0) {//row
            int r = axis;
            min = matrix[r][0];
            for (int i = 0; i < col; ++i) {
                if (min > matrix[r][i]) {
                    min = matrix[r][i];
                }
            }
        }
        else {//type==1, col
            int c = axis;
            min = matrix[0][c];
            for (int i = 0; i < row; ++i) {
                if (min > matrix[i][c]) {
                    min = matrix[i][c];
                }
            }
        }
        return min;
    }
    else return 0;
}
//sum
template<typename T>
T Matrix<T>::sum() {
    T sum = 0;
    if (row > 0 && col > 0) {
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                sum += matrix[i][j];
            }
        }
    }
    else cout << "The matrix isn't initialized." << endl;
    return sum;
}
template<typename T>
T Matrix<T>::sum(int type, int axis) {
    T sum = 0;
    if (check_type_axis(type, axis)) {
        if (type == 0) {//row
            int r = axis;
            for (int i = 0; i < col; ++i) {
                sum += matrix[r][i];
            }
        }
        else {//type==1, col
            int c = axis;
            for (int i = 0; i < row; ++i) {
                sum += matrix[i][c];
            }
        }
    }
    return sum;
}
//average value
template<typename T>
T Matrix<T>::avg() {
    T s = sum();
    if (row > 0 && col > 0) {
        T ele = row * col;
        return s / ele;
    }
    else return 0;
}
template<typename T>
T Matrix<T>::avg(int type, int axis) {
    T s = sum(type, axis);
    if (row > 0 && col > 0) {
        if (type == 0) {//row
            T c = col;
            return s / c;
        }
        else {//type==1, col
            T r = row;
            return s / r;
        }
    }
    else return 0;
}

//check square
template<typename T>
bool Matrix<T>::check_square() {
    if (un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return false;
    }
    else {
        if (row == col) return true;
        else {
            cout << "The matrix isn't a square matrix." << endl;
            return false;
        }
    }
}
//trace
template<typename T>
T Matrix<T>::trace() {
    if (check_square()) {
        T tr = 0;
        for (int i = 0; i < row; ++i) {
            tr += matrix[i][i];
        }
        return tr;
    }
    else return 0;
}
//determinant
template<typename T>
T Matrix<T>::get_det(vector<vector<T>> m, int n) {
    if (n == 1) return m[0][0];
    T det = 0;
    vector<vector<T>> tmp;
    //initialize
    for (int j = 0; j < n - 1; ++j) {
        vector<T> r;
        for (int k = 0; k < n - 1; ++k) {
            r.emplace_back(0);
        }
        tmp.emplace_back(r);
    }
    int i;
    for (i = 0; i < n; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            for (int k = 0; k < n - 1; ++k) {
                if (k >= i) tmp[j][k] = m[j + 1][k + 1];
                else tmp[j][k] = m[j + 1][k];
            }
        }
        T det_tmp = get_det(tmp, n - 1);
        if (i % 2 == 0) det += m[0][i] * det_tmp;
        else det -= m[0][i] * det_tmp;
    }
    return det;
}
template<typename T>
T Matrix<T>::determinant() {
    if (check_square()) {
        T det = get_det(matrix, row);
        return det;
    }
    else return 0;
}
//adjugate matrix 伴随矩阵 A*
template<typename T>
Matrix<T> Matrix<T>::adjugate() {
    if (check_square()) {
        int n = row;
        if (n == 1) {
            vector<vector<T>> m;
            vector<T> r;
            r.emplace_back(1);
            m.emplace_back(r);
            return Matrix<T>(m);
        }
        vector<vector<T>> adj;
        vector<vector<T>> tmp;
        //initialize
        for (int i = 0; i < n - 1; ++i) {
            vector<T> r;
            for (int j = 0; j < n - 1; ++j) {
                r.emplace_back(0);
            }
            tmp.emplace_back(r);
        }
        for (int i = 0; i < n; ++i) {
            vector<T> r;
            for (int j = 0; j < n; ++j) {
                r.emplace_back(0);
            }
            adj.emplace_back(r);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n - 1; ++k) {
                    for (int l = 0; l < n - 1; ++l) {
                        int k1 = (k >= i) ? k + 1 : k;
                        int l1 = (l >= j) ? l + 1 : l;
                        tmp[k][l] = matrix[k1][l1];
                    }
                }
                adj[j][i] = get_det(tmp, n - 1);//转置
                if ((i + j) % 2 == 1) adj[j][i] = -adj[j][i];
            }
        }
        return Matrix<T>(adj);
    }
    else return Matrix<T>(0);
}
//inverse AA*=|A|I
template<typename T>
bool Matrix<T>::has_inverse() {
    if (check_square()) {
        T det = get_det(matrix, row);
        if (det == 0) {
            cout << "The determinant of the matrix is 0, and the matrix has no inverse." << endl;
            return false;
        }
        else return true;
    }
    else return false;
}
template<typename T>
Matrix<T> Matrix<T>::inverse() {
    string T_str = typeid(T).name();
    if (check_df(T_str) == false) {
        cout << "Warning: T is not a floating point number and it will cause precision loss, so it can not return correct result." << endl;
    }
    if (has_inverse()) {
        T det = get_det(matrix, row);
        Matrix<T> adj = adjugate();
        Matrix<T> inv = adj / det;
        return inv;
    }
    else return Matrix<T>(0);
}
//eigenvalue and eigenvector, QR分解
//2的范数
template<typename T>
T Matrix<T>::norm2() {
    if (un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return 0;
    }
    T norm = 0;
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            norm += matrix[i][j] * matrix[i][j];
        }
    }
    return sqrt(norm);
}
template<typename T>
T Matrix<T>::norm2(Matrix<T>& m) {
    if (m.un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return 0;
    }
    T norm = 0;
    for (int i = 0; i < m.row; ++i) {
        for (int j = 0; j < m.col; ++j) {
            norm += m.matrix[i][j] * m.matrix[i][j];
        }
    }
    return sqrt(norm);
}
template<typename T>
T Matrix<T>::norm2(vector<T>& v) {
    if (v.size() == 0) {
        cout << "The vector isn't initialized." << endl;
        return 0;
    }
    T norm = 0;
    for (int i = 0; i < v.size(); ++i) {
        norm += v[i] * v[i];
    }
    return sqrt(norm);
}
template<typename T>
void Matrix<T>::QR(Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R) {
    int n = A.row;
    vector<T> a, b;
    for (int i = 0; i < n; ++i) {
        a.emplace_back(0);
        b.emplace_back(0);
    }
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            a[i] = A.matrix[i][j];
            b[i] = a[i];//第j列
        }
        for (int i = 0; i < j; ++i) {
            R.matrix[i][j] = 0;
            for (int k = 0; k < n; ++k) {
                R.matrix[i][j] += a[k] * Q.matrix[k][i];
            }
            for (int k = 0; k < n; ++k) {
                b[k] -= R.matrix[i][j] * Q.matrix[k][i];
            }
        }
        T norm = norm2(b);
        R.matrix[j][j] = norm;
        for (int i = 0; i < n; ++i) {
            Q.matrix[i][j] = b[i] / norm;
        }
    }
}
template<typename T>
vector<T> Matrix<T>::eigenvalues() {
    string T_str = typeid(T).name();
    if (check_df(T_str) == false) {
        cout << "Warning: T is not a floating point number and it will cause precision loss, so it can not return correct result." << endl;
    }
    vector<T> evalue;
    if (check_square()) {
        if (determinant() != 0) {
            int n = row;
            const int accu = (1000 > n) ? 1000 : n;
            vector<vector<T>> tmp(n, vector<T>(n));
            vector<vector<T>> q(n, vector<T>(n));
            vector<vector<T>> r(n, vector<T>(n));
            tmp = matrix;
            Matrix<T> TMP(tmp);
            Matrix<T> Q(q);
            Matrix<T> R(r);
            for (int i = 0; i < accu; ++i) {
                QR(TMP, Q, R);
                TMP = R * Q;
            }
            for (int i = 0; i < n; ++i) {
                evalue.emplace_back(TMP.matrix[i][i]);
            }
            vector<T> evalue_copy;
            //去掉重复特征值
            sort(evalue.begin(), evalue.end());
            evalue_copy.emplace_back(evalue[0]);
            for (int i = 1; i < evalue.size(); ++i) {
                if (evalue[i] != evalue[i - 1]) {
                    evalue_copy.emplace_back(evalue[i]);
                }
            }
            return evalue_copy;
        }
        else {
            cout << "The determinant of the matrix is 0, eigenvalues can not be found by QR decomposition." << endl;
            cout << "The matrix must have eigenvalue 0." << endl;
            evalue.emplace_back(0);
            return evalue;
        }
    }
    return evalue;
}
//eigenvectors
template<typename T>
vector<vector<T>> Matrix<T>::eigenvectors() {
    vector<vector<T>> evector;
    if (check_square()) {
        vector<T> evalue = eigenvalues();
        for (int i = 0; i < evalue.size(); ++i) {
            T e = evalue[i];
            //A-eI
            int n = row;
            Matrix<T> tmp(matrix);
            int zero_flag = 0;
            for (int j = 0; j < n; ++j) {
                tmp.matrix[j][j] -= e;
                if (tmp.matrix[j][j] == 0) zero_flag = 1;
            }
            if (zero_flag == 1) {
                cout << "The elements on the main diagonal of the matrix contains 0 and can't use Gaussian method." << endl;
                return evector;
            }
            //高斯消元，将tmp化成上三角矩阵
            for (int j = 0; j < n - 1; ++j) {
                T md = tmp.matrix[j][j];
                for (int k = j; k < n; ++k) {
                    tmp.matrix[j][k] /= md;
                }
                for (int k = j + 1; k < n; ++k) {
                    T num = tmp.matrix[k][j];
                    for (int l = j; l < n; ++l) {
                        tmp.matrix[k][l] -= num * tmp.matrix[j][l];
                    }
                }
            }
            vector<T> ev;//单位向量
            //initialize 1*n
            for (int j = 0; j < n; ++j) {
                ev.emplace_back(0);
            }
            //A*ev=0 从最后一行算到第一行给ev[j]赋值
            ev[n - 1] = 1;
            T ev_sum = 1;
            for (int j = n - 2; j >= 0; --j) {
                T sum = 0;
                for (int k = j + 1; k < n; ++k) {
                    sum += tmp.matrix[j][k] * ev[k];
                }
                ev[j] = -sum / tmp.matrix[j][j];
                ev_sum += ev[j] * ev[j];
            }
            ev_sum = sqrt(ev_sum);
            for (int j = 0; j < n; ++j) {
                ev[j] /= ev_sum;
            }
            evector.emplace_back(ev);
        }
        return evector;
    }
    return evector;
}
template<typename T>
void Matrix<T>::print_eigenvectors() {
    vector<vector<T>> evector = eigenvectors();
    for (int i = 0; i < evector.size(); ++i) {
        cout << "[";
        for (int j = 0; j < evector[0].size() - 1; ++j) {
            cout << evector[i][j] << ",";
        }
        cout << evector[i][evector[0].size() - 1] << "]" << endl;
    }
}
//reshape 整型
template<typename T>
Matrix<T> Matrix<T>::reshape(Matrix<T> m, int r, int c) {
    if (m.un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return Matrix<T>(0);
    }
    vector<vector<T>> res;//r*c
    if (m.row * m.col == r * c) {
        //initialize
        for (int i = 0; i < r; ++i) {
            vector<T> ro;
            for (int j = 0; j < c; ++j) {
                ro.emplace_back(0);
            }
            res.emplace_back(ro);
        }
        for (int i = 0; i < r * c; ++i) {
            res[i % r][i / r] = m.matrix[i % m.row][i / m.row];
        }
        return Matrix<T>(res);
    }
    else {
        cout << "The sizes of row and column of new matrix are mismatched." << endl;
        return Matrix<T>(0);
    }
}
template<typename T>
Matrix<T> Matrix<T>::reshape(int r, int c) {
    return reshape(Matrix<T>(matrix), r, c);
}
//slicing 切片
template<typename T>
Matrix<T> Matrix<T>::slicing(Matrix<T> m, int row_begin, int row_end, int col_begin, int col_end, int row_interval, int col_interval) {
    if (m.un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return Matrix<T>(0);
    }
    if ((0 <= row_begin && row_begin <= row_end && row_end < m.row) && (0 <= col_begin && col_begin <= col_end && col_end < m.col)
        && (0 <= row_interval && row_interval < row_end - row_begin + 1) && (0 <= col_interval && col_interval < col_end - col_begin + 1)) {
        Matrix<T> sli_1(row_end - row_begin + 1, col_end - col_begin + 1);
        for (int i = 0; i < sli_1.row; ++i) {
            for (int j = 0; j < sli_1.col; ++j) {
                sli_1.matrix[i][j] = m.matrix[row_begin + i][col_begin + j];
            }
        }
        int row_period = row_interval + 1;
        int col_period = col_interval + 1;
        int row_time = sli_1.row / row_period;
        int col_time = sli_1.col / col_period;
        int sli_2_row = ((sli_1.row % row_period == 0)) ? row_time : (row_time + 1);
        int sli_2_col = ((sli_1.col % col_period) == 0) ? col_time : (col_time + 1);
        Matrix<T> sli_2(sli_2_row, sli_2_col);
        int i = 0, j = 0, k = 0, t = 0;
        for (i = 0, k = 0; i < sli_2_row; ++i, k = k + row_period) {
            for (j = 0, t = 0; j < sli_2_col; ++j, t = t + col_period) {
                sli_2.matrix[i][j] = sli_1.matrix[k][t];
            }
        }
        return sli_2;
    }
    else {
        cout << "The range of input is error." << endl;
        return Matrix<T>(0);
    }
}
template<typename T>
Matrix<T> Matrix<T>::row_slicing(Matrix<T> m, int row_begin, int row_end, int interval) {
    return slicing(m, row_begin, row_end, 0, m.col - 1, interval, 0);
}
template<typename T>
Matrix<T> Matrix<T>::col_slicing(Matrix<T> m, int col_begin, int col_end, int interval) {
    return slicing(m, 0, m.row - 1, col_begin, col_end, 0, interval);
}
template<typename T>
Matrix<T> Matrix<T>::slicing(int row_begin, int row_end, int col_begin, int col_end, int row_interval, int col_interval) {
    return slicing(Matrix<T>(matrix), row_begin, row_end, col_begin, col_end, row_interval, col_interval);
}
template<typename T>
Matrix<T> Matrix<T>::row_slicing(int row_begin, int row_end, int interval) {
    return row_slicing(Matrix<T>(matrix), row_begin, row_end, interval);
}
template<typename T>
Matrix<T> Matrix<T>::col_slicing(int col_begin, int col_end, int interval) {
    return col_slicing(Matrix<T>(matrix), col_begin, col_end, interval);
}
//convolution 卷积
template<typename T>
Matrix<T> Matrix<T>::convolution(Matrix<T> m1, Matrix<T> m2) {
    if (m1.un_initialized() || m2.un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return Matrix<T>(0);
    }
    Matrix<T> big;
    Matrix<T> small;
    big = (m1.row * m1.col > m2.row * m2.col) ? m1 : m2;
    small = (m1.row * m1.col <= m2.row * m2.col) ? m1 : m2;
    Matrix<T> extend(big.row + 2 * (small.row - 1), big.col + 2 * (small.col - 1));
    Matrix<T> full(big.row + small.row - 1, big.col + small.col - 1);
    Matrix<T> small_180(small.row, small.col);
    for (int i = 0; i < small.row; ++i) {
        for (int j = 0; j < small.col; ++j) {
            small_180.matrix[small.row - 1 - i][small.col - 1 - j] = small.matrix[i][j];
        }
    }
    for (int i = 0; i < big.row; ++i) {
        for (int j = 0; j < big.col; ++j) {
            extend.matrix[small.row - 1 + i][small.col - 1 + j] = big.matrix[i][j];
        }
    }
    //small_180矩阵起点在extend矩阵的位置，对应full矩阵的i,j
    for (int i = 0; i < full.row; ++i) {
        for (int j = 0; j < full.col; ++j) {
            T sum = 0;
            for (int k = 0; k < small_180.row; ++k) {
                for (int l = 0; l < small_180.col; ++l) {
                    sum += small_180.matrix[k][l] * extend.matrix[i + k][j + l];
                }
            }
            full.matrix[i][j] = sum;
        }
    }
    return full;
}

//Matrix<complex<D>>
//检查初始化
template<typename D>
bool Matrix<complex<D>>::un_initialized() {
    if (row == 0 || col == 0)return true;
    else return false;
}
//overload cout<<
template<typename D>
ostream& operator<<(ostream& os, const Matrix<complex<D>>& other) {
    if (other.row == 0 || other.col == 0) {
        os << "The matrix isn't initialized." << endl;
    }
    else {
        for (int i = 0; i < other.row; ++i) {
            for (int j = 0; j < other.col; ++j) {
                os << other.matrix[i][j] << " ";
            }
            os << endl;
        }
    }
    return os;
}
template<typename D>
ostream& operator<<(ostream& os, const vector<complex<D>>& other) {
    if (other.size() == 0) {
        os << "The vector isn't initialized." << endl;
    }
    else {
        os << "[";
        for (int i = 0; i < other.size() - 1; ++i) {
            os << other[i] << ",";
        }
        os << other[other.size() - 1] << "]";
    }
    return os;
}
//overload =
template<typename D>
Matrix<complex<D>>& Matrix<complex<D>>::operator=(const Matrix<complex<D>>& next) {
    if (this != &next) {
        this->row = next.row;
        this->col = next.col;
        this->matrix = next.matrix;
    }
    return *this;
}
//A+B
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::operator+(const Matrix<complex<D>>& next) {
    try {
        if (un_initialized() || next.row == 0 || next.col == 0) {
            throw 0;
        }
        vector<vector<complex<D>>> add;
        //矩阵大小匹配
        if (this->row == next.row && this->col == next.col) {
            for (int i = 0; i < this->row; ++i) {
                vector<complex<D>> r;
                for (int j = 0; j < this->col; ++j) {
                    r.emplace_back(this->matrix[i][j] + next.matrix[i][j]);
                }
                add.emplace_back(r);
            }
            return Matrix<complex<D>>(add);
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        //矩阵大小不匹配
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
//v1+v2
template<typename D>
vector<complex<D>> operator+(const vector<complex<D>>& v1, const vector<complex<D>>& v2) {
    try {
        vector<complex<D>> v;
        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); ++i) {
                v.emplace_back(v1[i] + v2[i]);
            }
            return v;
        }
        throw 0;
    }
    catch (int i)
    {
        cout << "The sizes of vectors are mismatched." << endl;
    }
}
//A-B
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::operator-(const Matrix<complex<D>>& next) {
    try {
        if (un_initialized() || next.row == 0 || next.col == 0) {
            throw 0;
        }
        vector<vector<complex<D>>> sub;
        //矩阵大小匹配
        if (this->row == next.row && this->col == next.col) {
            for (int i = 0; i < this->row; ++i) {
                vector<complex<D>> r;
                for (int j = 0; j < this->col; ++j) {
                    r.emplace_back(this->matrix[i][j] - next.matrix[i][j]);
                }
                sub.emplace_back(r);
            }
            return Matrix<complex<D>>(sub);
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        //矩阵大小不匹配
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
//v1-v2
template<typename D>
vector<complex<D>> operator-(vector<complex<D>>& v1, vector<complex<D>>& v2) {
    try {
        if (v1.size() == 0 || v2.size() == 0) {
            throw 0;
        }
        vector<complex<D>> v;
        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); ++i) {
                v.emplace_back(v1[i] - v2[i]);
            }
            return v;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of vectors are mismatched." << endl;
    }
}
//scalar multiplication
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::operator*(complex<D> val) {
    try {
        if (un_initialized()) {
            throw 0;
        }
        vector<vector<complex<D>>> mul;
        for (int i = 0; i < this->row; ++i) {
            vector<complex<D>> r;
            for (int j = 0; j < this->col; ++j) {
                this->matrix[i][j] *= val;
                r.emplace_back(this->matrix[i][j] * val);
            }
            mul.emplace_back(r);
        }
        return Matrix<complex<D>>(mul);
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
}
template<typename D>
Matrix<complex<D>> operator*(complex<D> val, const Matrix<complex<D>>& other) {
    try {
        if (other.row == 0 || other.col == 0) {
            throw 0;
        }
        vector<vector<complex<D>>> mul;
        for (int i = 0; i < other.row; ++i) {
            vector<complex<D>> r;
            for (int j = 0; j < other.col; ++j) {
                r.emplace_back(val * other.matrix[i][j]);
            }
            mul.emplace_back(r);
        }
        return Matrix<complex<D>>(mul);
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }

}
template<typename D>
vector<complex<D>> operator*(vector<complex<D>>& v, complex<D> val) {
    try {
        if (v.size() == 0) {
            throw 0;
        }
        vector<complex<D>> mul;
        for (int i = 0; i < v.size(); ++i) {
            mul.emplace_back(v[i] * val);
        }
        return mul;
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }

}
template<typename D>
vector<complex<D>> operator*(complex<D> val, vector<complex<D>>& v) {
    try {
        if (v.size() == 0) {
            throw 0;
        }
        vector<complex<D>> mul;
        for (int i = 0; i < v.size(); ++i) {
            mul.emplace_back(v[i] * val);
        }
        return mul;
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }

}
//scalar division
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::operator/(complex<D> val) {
    try {
        if (un_initialized()) {
            throw 0;
        }
        vector<vector<complex<D>>> div;
        for (int i = 0; i < this->row; ++i) {
            vector<complex<D>> r;
            for (int j = 0; j < this->col; ++j) {
                r.emplace_back(this->matrix[i][j] / val);
            }
            div.emplace_back(r);
        }
        return Matrix<complex<D>>(div);
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }

}
template<typename D>
vector<complex<D>> operator/(vector<complex<D>>& v, complex<D> val) {
    try {
        if (v.size() == 0) {
            throw 0;
        }
        vector<complex<D>> div;
        for (int i = 0; i < v.size(); ++i) {
            div.emplace_back(v[i] / val);
        }
        return div;
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }

}
//transposition
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::transposition(const Matrix<complex<D>>& m) {
    try {
        if (m.row == 0 || m.col == 0) {
            throw 0;
        }
        vector<vector<complex<D>>> trans;
        for (int i = 0; i < m.col; ++i) {
            vector<complex<D>> trans_r;
            for (int j = 0; j < m.row; ++j) {
                trans_r.emplace_back(m.matrix[j][i]);
            }
            trans.emplace_back(trans_r);
        }
        return Matrix<complex<D>>(trans);
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
}
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::transposition() {
    try {
        if (un_initialized()) {
            throw 0;
        }
        vector<vector<complex<D>>> trans;
        for (int i = 0; i < this->col; ++i) {
            vector<complex<D>> trans_r;
            for (int j = 0; j < this->row; ++j) {
                trans_r.emplace_back(this->matrix[j][i]);
            }
            trans.emplace_back(trans_r);
        }
        return Matrix<complex<D>>(trans);
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }

};
//conjugation 共轭
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::conjugation(const Matrix<complex<D>>& m) {
    Matrix<complex<D>> coj = transposition(m);
    for (int i = 0; i < coj.row; ++i) {
        for (int j = 0; j < coj.col; ++j) {
            coj.matrix[i][j] = conj(coj.matrix[i][j]);
        }
    }
    return coj;
}
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::conjugation() {
    Matrix<complex<D>> coj = transposition();
    for (int i = 0; i < coj.row; ++i) {
        for (int j = 0; j < coj.col; ++j) {
            coj.matrix[i][j] = conj(coj.matrix[i][j]);
        }
    }
    return coj;
}
//element-wise multiplication
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::element_wise_mul(const Matrix<complex<D>>& m) {
    try {
        if (m.row == 0 || m.col == 0 || row == 0 || col == 0) {
            throw 0;
        }
        vector<vector<complex<D>>> mul;
        if (row == m.row && col == m.col) {
            for (int i = 0; i < row; ++i) {
                vector<complex<D>> r;
                for (int j = 0; j < col; ++j) {
                    r.emplace_back(matrix[i][j] * m.matrix[i][j]);
                }
                mul.emplace_back(r);
            }
            return Matrix<complex<D>>(mul);
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of matrices are mismatched." << endl;
    }

}
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::element_wise_mul(const Matrix<complex<D>>& m1, const Matrix<complex<D>>& m2) {
    try {
        if (m1.row == 0 || m1.col == 0 || m2.row == 0 || m2.col == 0) {
            throw 0;
        }
        vector<vector<complex<D>>> mul;
        if (m1.row == m2.row && m1.col == m2.col) {
            for (int i = 0; i < m1.row; ++i) {
                vector<complex<D>> r;
                for (int j = 0; j < m1.col; ++j) {
                    r.emplace_back(m1.matrix[i][j] * m2.matrix[i][j]);
                }
                mul.emplace_back(r);
            }
            return Matrix<complex<D>>(mul);
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
template<typename D>
vector<complex<D>> Matrix<complex<D>>::element_wise_mul(vector<complex<D>>& v1, vector<complex<D>>& v2) {
    try {
        if (v1.size() == 0 || v2.size() == 0) {
            throw 0;
        }
        vector<complex<D>> mul;
        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); ++i) {
                mul.emplace_back(v1[i] * v2[i]);
            }
            return mul;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of vectors are mismatched." << endl;
    }
}
//matrix-matrix multiplication
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::operator*(const Matrix<complex<D>>& m) {
    try {
        if (row == 0 || col == 0 || m.row == 0 || m.col == 0) {
            throw 0;
        }
        vector<vector<complex<D>>> mul;
        if (col == m.row) {
            //initialize
            for (int i = 0; i < row; ++i) {
                vector<complex<D>> r;
                for (int j = 0; j < m.col; ++j) {
                    r.emplace_back(0);
                }
                mul.emplace_back(r);
            }
            for (int i = 0; i < row; ++i) {
                for (int j = 0; j < m.col; ++j) {
                    for (int k = 0; k < col; ++k) {
                        mul[i][j] += matrix[i][k] * m.matrix[k][j];
                    }
                }
            }
            return Matrix<complex<D>>(mul);
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
//matrix-vector multiplication
template<typename D>
vector<complex<D>> Matrix<complex<D>>::operator*(vector<complex<D>>& vec) {
    try {
        if (row == 0 || col == 0 || vec.size() == 0) {
            throw 0;
        }
        vector<complex<D>> v;
        if (col == vec.size()) {
            //initialize
            for (int i = 0; i < row; ++i) {
                v.emplace_back(0);
            }
            for (int i = 0; i < row; ++i) {
                for (int j = 0; j < col; ++j) {
                    v[i] += matrix[i][j] * vec[j];
                }
            }
            return v;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix or the vector isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The size of matrix and the size of vector are mismatched." << endl;
    }
}
template<typename D>
vector<complex<D>> operator*(vector<complex<D>>& vec, const Matrix<complex<D>>& other) {
    try {
        if (other.row == 0 || other.col == 0 || vec.size() == 0) {
            throw 0;
        }
        vector<complex<D>> v;
        if (vec.size() == other.row) {
            //initialize
            for (int i = 0; i < other.col; ++i) {
                v.emplace_back(0);
            }
            for (int i = 0; i < other.col; ++i) {
                for (int j = 0; j < other.row; ++j) {
                    v[i] += other.matrix[j][i] * vec[j];
                }
            }
            return v;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix or the vector isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The size of matrix and the size of vector are mismatched." << endl;
    }
}
template<typename D>
complex<D> operator*(vector<complex<D>> v1, vector<complex<D>> v2) {
    try {
        if (v1.size() == 0 || v2.size() == 0) {
            throw 0;
        }
        complex<D> mul = 0;
        if (v1.size() == v2.size()) {
            for (int i = 0; i < v1.size(); ++i) {
                mul += v1[i] * v2[i];
            }
            return mul;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The vector isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of vectors are mismatched." << endl;
    }
}
//dot product 内积
template<typename D>
complex<D> Matrix<complex<D>>::dot_product(Matrix<complex<D>>& m1, Matrix<complex<D>>& m2) {
    try {
        if (m1.un_initialized() || m2.un_initialized()) {
            throw 0;
        }
        complex<D> dot_prod = 0;
        if (m1.row == m2.row && m1.col == m2.col) {
            for (int i = 0; i < m1.row; ++i) {
                for (int j = 0; j < m1.col; ++j) {
                    dot_prod += m1.matrix[i][j] * m2.matrix[i][j];
                }
            }
            return dot_prod;
        }
        throw 'a';
    }
    catch (int i) {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a) {
        cout << "The sizes of matrices are mismatched." << endl;
    }

}
template<typename D>
complex<D> Matrix<complex<D>>::dot_product(Matrix<complex<D>>& m) {
    try {
        if (m.un_initialized() || un_initialized()) {
            throw 0;
        }
        complex<D> dot_prod = 0;
        if (row == m.row && col == m.col) {
            for (int i = 0; i < row; ++i) {
                for (int j = 0; j < col; ++j) {
                    dot_prod += matrix[i][j] * m.matrix[i][j];
                }
            }
            return dot_prod;
        }
        throw 'a';
    }
    catch (int i)
    {
        cout << "The matrix isn't initialized." << endl;
    }
    catch (char a)
    {
        cout << "The sizes of matrices are mismatched." << endl;
    }
}
template<typename D>
complex<D> Matrix<complex<D>>::dot_product(vector<complex<D>>& v1, vector<complex<D>>& v2) {
    return v1 * v2;
}
//cross product n=3
template<typename D>
vector<complex<D>> Matrix<complex<D>>::cross_product(vector<complex<D>>& v1, vector<complex<D>>& v2) {
    try {
        vector<complex<D>> prod;
        if (v1.size() == 3 && v2.size() == 3) {
            prod.emplace_back(v1[1] * v2[2] - v1[2] * v2[1]);
            prod.emplace_back(v1[2] * v2[0] - v1[0] * v2[2]);
            prod.emplace_back(v1[0] * v2[1] - v1[1] * v2[0]);
            return prod;
        }
        throw 0;
    }
    catch (int i)
    {
        cout << "The sizes of vectors are not 3 and have no cross product." << endl;
    }
}
//check type and axis
template<typename D>
bool Matrix<complex<D>>::check_type_axis(int type, int axis) {
    if (un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return false;
    }
    else {
        if (type == 0) {//row
            if (axis >= 0 && axis < row) return true;
            else {
                cout << "Invalid axis." << endl;
                return false;
            }
        }
        else if (type == 1) {//col
            if (axis >= 0 && axis < col) return true;
            else {
                cout << "Invalid axis." << endl;
                return false;
            }
        }
        else {
            cout << "Invalid type." << endl;
            return false;
        }
    }
}
//sum
template<typename D>
complex<D> Matrix<complex<D>>::sum() {
    complex<D> sum = 0;
    if (row > 0 && col > 0) {
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                sum += matrix[i][j];
            }
        }
    }
    else cout << "The matrix isn't initialized." << endl;
    return sum;
}
template<typename D>
complex<D> Matrix<complex<D>>::sum(int type, int axis) {
    complex<D> sum = 0;
    if (check_type_axis(type, axis)) {
        if (type == 0) {//row
            int r = axis;
            for (int i = 0; i < col; ++i) {
                sum += matrix[r][i];
            }
        }
        else {//type==1, col
            int c = axis;
            for (int i = 0; i < row; ++i) {
                sum += matrix[i][c];
            }
        }
    }
    return sum;
}
//average value
template<typename D>
complex<D> Matrix<complex<D>>::avg() {
    complex<D> s = sum();
    if (row > 0 && col > 0) {
        complex<D> ele = row * col;
        return s / ele;
    }
    else return 0;
}
template<typename D>
complex<D> Matrix<complex<D>>::avg(int type, int axis) {
    complex<D> s = sum(type, axis);
    if (row > 0 && col > 0) {
        if (type == 0) {//row
            complex<D> c = col;
            return s / c;
        }
        else {//type==1, col
            complex<D> r = row;
            return s / r;
        }
    }
    else return 0;
}

//check square
template<typename D>
bool Matrix<complex<D>>::check_square() {
    if (un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return false;
    }
    else {
        if (row == col) return true;
        else {
            cout << "The matrix isn't a square matrix." << endl;
            return false;
        }
    }
}
//trace
template<typename D>
complex<D> Matrix<complex<D>>::trace() {
    if (check_square()) {
        complex<D> tr = 0;
        for (int i = 0; i < row; ++i) {
            tr += matrix[i][i];
        }
        return tr;
    }
    else return 0;
}
//determinant
template<typename D>
complex<D> Matrix<complex<D>>::get_det(vector<vector<complex<D>>> m, int n) {
    if (n == 1) return m[0][0];
    complex<D> det = 0;
    vector<vector<complex<D>>> tmp;
    //initialize
    for (int j = 0; j < n - 1; ++j) {
        vector<complex<D>> r;
        for (int k = 0; k < n - 1; ++k) {
            r.emplace_back(0);
        }
        tmp.emplace_back(r);
    }
    int i;
    for (i = 0; i < n; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            for (int k = 0; k < n - 1; ++k) {
                if (k >= i) tmp[j][k] = m[j + 1][k + 1];
                else tmp[j][k] = m[j + 1][k];
            }
        }
        complex<D> det_tmp = get_det(tmp, n - 1);
        if (i % 2 == 0) det += m[0][i] * det_tmp;
        else det -= m[0][i] * det_tmp;
    }
    return det;
}
template<typename D>
complex<D> Matrix<complex<D>>::determinant() {
    if (check_square()) {
        complex<D> det = get_det(matrix, row);
        return det;
    }
    else return 0;
}
//adjugate matrix 伴随矩阵 A*
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::adjugate() {
    if (check_square()) {
        int n = row;
        if (n == 1) {
            vector<vector<complex<D>>> m;
            vector<complex<D>> r;
            r.emplace_back(1);
            m.emplace_back(r);
            return Matrix<complex<D>>(m);
        }
        vector<vector<complex<D>>> adj;
        vector<vector<complex<D>>> tmp;
        //initialize
        for (int i = 0; i < n - 1; ++i) {
            vector<complex<D>> r;
            for (int j = 0; j < n - 1; ++j) {
                r.emplace_back(0);
            }
            tmp.emplace_back(r);
        }
        for (int i = 0; i < n; ++i) {
            vector<complex<D>> r;
            for (int j = 0; j < n; ++j) {
                r.emplace_back(0);
            }
            adj.emplace_back(r);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n - 1; ++k) {
                    for (int l = 0; l < n - 1; ++l) {
                        int k1 = (k >= i) ? k + 1 : k;
                        int l1 = (l >= j) ? l + 1 : l;
                        tmp[k][l] = matrix[k1][l1];
                    }
                }
                adj[j][i] = get_det(tmp, n - 1);//转置
                if ((i + j) % 2 == 1) adj[j][i] = -adj[j][i];
            }
        }
        return Matrix<complex<D>>(adj);
    }
    else return Matrix<complex<D>>(0);
}
//inverse AA*=|A|I
template<typename D>
bool Matrix<complex<D>>::has_inverse() {
    if (check_square()) {
        complex<D> det = get_det(matrix, row);
        if (det == 0) {
            cout << "The determinant of the matrix is 0, and the matrix has no inverse." << endl;
            return false;
        }
        else return true;
    }
    else return false;
}
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::inverse() {
    string T_str = typeid(complex<D>).name();
    if (has_inverse()) {
        complex<D> det = get_det(matrix, row);
        Matrix<complex<D>> adj = adjugate();
        Matrix<complex<D>> inv = adj / det;
        return inv;
    }
    else return Matrix<complex<D>>(0);
}
//2的范数
template<typename D>
D Matrix<complex<D>>::norm2() {
    if (un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return 0;
    }
    D norm = 0;
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            norm += matrix[i][j].imag() * matrix[i][j].imag() + matrix[i][j].real() * matrix[i][j].real();//|x|^2
        }
    }
    return sqrt(norm);
}
template<typename D>
D Matrix<complex<D>>::norm2(Matrix<complex<D>>& m) {
    if (m.un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return 0;
    }
    D norm = 0;
    for (int i = 0; i < m.row; ++i) {
        for (int j = 0; j < m.col; ++j) {
            norm += m.matrix[i][j].imag() * m.matrix[i][j].imag() + m.matrix[i][j].real() * m.matrix[i][j].real();
        }
    }
    return sqrt(norm);
}
template<typename D>
D Matrix<complex<D>>::norm2(vector<complex<D>>& v) {
    if (v.size() == 0) {
        cout << "The vector isn't initialized." << endl;
        return 0;
    }
    D norm = 0;
    for (int i = 0; i < v.size(); ++i) {
        norm += v[i].imag() * v[i].imag() + v[i].real() * v[i].real();
    }
    return sqrt(norm);
}
//reshape 整型
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::reshape(Matrix<complex<D>> m, int r, int c) {
    if (m.un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return Matrix<complex<D>>(0);
    }
    vector<vector<complex<D>>> res;//r*c
    if (m.row * m.col == r * c) {
        //initialize
        for (int i = 0; i < r; ++i) {
            vector<complex<D>> ro;
            for (int j = 0; j < c; ++j) {
                ro.emplace_back(0);
            }
            res.emplace_back(ro);
        }
        for (int i = 0; i < r * c; ++i) {
            res[i % r][i / r] = m.matrix[i % m.row][i / m.row];
        }
        return Matrix<complex<D>>(res);
    }
    else {
        cout << "The sizes of row and column of new matrix are mismatched." << endl;
        return Matrix<complex<D>>(0);
    }
}
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::reshape(int r, int c) {
    return reshape(Matrix<complex<D>>(matrix), r, c);
}
//slicing 切片
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::slicing(Matrix<complex<D>> m, int row_begin, int row_end, int col_begin, int col_end, int row_interval, int col_interval) {
    if (m.un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return Matrix<complex<D>>(0);
    }
    if ((0 <= row_begin && row_begin <= row_end && row_end < m.row) && (0 <= col_begin && col_begin <= col_end && col_end < m.col)
        && (0 <= row_interval && row_interval < row_end - row_begin + 1) && (0 <= col_interval && col_interval < col_end - col_begin + 1)) {
        Matrix<complex<D>> sli_1(row_end - row_begin + 1, col_end - col_begin + 1);
        for (int i = 0; i < sli_1.row; ++i) {
            for (int j = 0; j < sli_1.col; ++j) {
                sli_1.matrix[i][j] = m.matrix[row_begin + i][col_begin + j];
            }
        }
        int row_period = row_interval + 1;
        int col_period = col_interval + 1;
        int row_time = sli_1.row / row_period;
        int col_time = sli_1.col / col_period;
        int sli_2_row = ((sli_1.row % row_period == 0)) ? row_time : (row_time + 1);
        int sli_2_col = ((sli_1.col % col_period) == 0) ? col_time : (col_time + 1);
        Matrix<complex<D>> sli_2(sli_2_row, sli_2_col);
        int i = 0, j = 0, k = 0, t = 0;
        for (i = 0, k = 0; i < sli_2_row; ++i, k = k + row_period) {
            for (j = 0, t = 0; j < sli_2_col; ++j, t = t + col_period) {
                sli_2.matrix[i][j] = sli_1.matrix[k][t];
            }
        }
        return sli_2;
    }
    else {
        cout << "The range of input is error." << endl;
        return Matrix<complex<D>>(0);
    }
}
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::row_slicing(Matrix<complex<D>> m, int row_begin, int row_end, int interval) {
    return slicing(m, row_begin, row_end, 0, m.col - 1, interval, 0);
}
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::col_slicing(Matrix<complex<D>> m, int col_begin, int col_end, int interval) {
    return slicing(m, 0, m.row - 1, col_begin, col_end, 0, interval);
}
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::slicing(int row_begin, int row_end, int col_begin, int col_end, int row_interval, int col_interval) {
    return slicing(Matrix<complex<D>>(matrix), row_begin, row_end, col_begin, col_end, row_interval, col_interval);
}
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::row_slicing(int row_begin, int row_end, int interval) {
    return row_slicing(Matrix<complex<D>>(matrix), row_begin, row_end, interval);
}
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::col_slicing(int col_begin, int col_end, int interval) {
    return col_slicing(Matrix<complex<D>>(matrix), col_begin, col_end, interval);
}
//convolution 卷积
template<typename D>
Matrix<complex<D>> Matrix<complex<D>>::convolution(Matrix<complex<D>> m1, Matrix<complex<D>> m2) {
    if (m1.un_initialized() || m2.un_initialized()) {
        cout << "The matrix isn't initialized." << endl;
        return Matrix<complex<D>>(0);
    }
    Matrix<complex<D>> big;
    Matrix<complex<D>> small;
    big = (m1.row * m1.col > m2.row * m2.col) ? m1 : m2;
    small = (m1.row * m1.col <= m2.row * m2.col) ? m1 : m2;
    Matrix<complex<D>> extend(big.row + 2 * (small.row - 1), big.col + 2 * (small.col - 1));
    Matrix<complex<D>> full(big.row + small.row - 1, big.col + small.col - 1);
    Matrix<complex<D>> small_180(small.row, small.col);
    for (int i = 0; i < small.row; ++i) {
        for (int j = 0; j < small.col; ++j) {
            small_180.matrix[small.row - 1 - i][small.col - 1 - j] = small.matrix[i][j];
        }
    }
    for (int i = 0; i < big.row; ++i) {
        for (int j = 0; j < big.col; ++j) {
            extend.matrix[small.row - 1 + i][small.col - 1 + j] = big.matrix[i][j];
        }
    }
    //small_180矩阵起点在extend矩阵的位置，对应full矩阵的i,j
    for (int i = 0; i < full.row; ++i) {
        for (int j = 0; j < full.col; ++j) {
            complex<D> sum = 0;
            for (int k = 0; k < small_180.row; ++k) {
                for (int l = 0; l < small_180.col; ++l) {
                    sum += small_180.matrix[k][l] * extend.matrix[i + k][j + l];
                }
            }
            full.matrix[i][j] = sum;
        }
    }
    return full;
}