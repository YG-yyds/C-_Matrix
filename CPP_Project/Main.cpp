#include <iostream>
#include <vector>
#include "Matrix.cpp"
using namespace std;

int main() {
	int n, m;
	n = 3;
	m = 2;
	vector<vector<double>> vv1(n, vector<double>(n));
	vector<vector<double>> vv2(n, vector<double>(n));
	vector<vector<double>> vv3(m, vector<double>(m));
	vector<vector<double>> vv4(n, vector<double>(m));
	vector<double> v1(n);
	vector<double> v2(n);
	vector<double> v3(m);
	vector<double> v_null;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			vv1[i][j] = (double)i * n + j;
		}
		v1[i] = (double)i;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			vv2[i][j] = (double)i * n + j + 1;
		}
		v2[i] = (double)i + 1;
	}

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			vv3[i][j] = (double)i * m + j;
		}
		v3[i] = (double)i + 1;
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			vv4[i][j] = (double)i * n + j;
		}
	}

	Matrix<double> m1(vv1);
	Matrix<double> m2(vv2);
	Matrix<double> m3(vv3);
	Matrix<double> m4(vv4);
	Matrix<double> m_null;

	// <<
	cout << "The empty matrix m_null:" << endl;
	cout << m_null;
	cout << "The first matrix m1:" << endl;
	cout << m1;
	cout << "The second matrix m2:" << endl;
	cout << m2;
	cout << "The third matrix m3:" << endl;
	cout << m3;
	cout << "The third matrix m4:" << endl;
	cout << m4;
	cout << endl;
	cout << "The empty vector v_null:" << endl;
	cout << v_null;
	cout << "The first vector v1:" << endl;
	cout << v1;
	cout << "The second vector v2:" << endl;
	cout << v2;
	cout << "The third vector v3:" << endl;
	cout << v3;

	cout << endl;

	// +
	cout << "m1 + m_null:" << endl;
	Matrix<double> m_plus1 = m1 + m_null;
	cout << "m1 + m3:" << endl;
	Matrix<double> m_plus2 = m1 + m3;
	cout << "m1 + m2:" << endl;
	Matrix<double> m_plus3 = m1 + m2;
	cout << m_plus3;
	cout << endl;
	cout << "v1 + v_null:" << endl;
	vector<double> v_plus1 = v1 + v_null;
	cout << "v1 + v3:" << endl;
	vector<double> v_plus2 = v1 + v3;
	cout << "v1 + v2:" << endl;
	vector<double> v_plus3 = v1 + v2;
	cout << v_plus3;

	cout << endl;

	// -
	cout << "m1 - m_null:" << endl;
	Matrix<double> m_minus1 = m1 - m_null;
	cout << "m1 - m3:" << endl;
	Matrix<double> m_minus2 = m1 - m3;
	cout << "m1 - m2:" << endl;
	Matrix<double> m_minus3 = m1 - m2;
	cout << m_minus3;
	cout << endl;
	cout << "v1 - v_null:" << endl;
	vector<double> v_minus1 = v1 - v_null;
	cout << "v1 - v3:" << endl;
	vector<double> v_minus2 = v1 - v3;
	cout << "v1 - v2:" << endl;
	vector<double> v_minus3 = v1 - v2;
	cout << v_minus3;

	cout << endl;

	// *
	cout << "m_null * 2:" << endl;
	Matrix<double> m_mul1 = m_null * 2.0;
	cout << "m1 * 2:" << endl;
	Matrix<double> m_mul2 = m1 * 2.0;
	cout << m_mul2;
	cout << "2 * m1:" << endl;
	Matrix<double> m_mul3 = 2.0 * m1;
	cout << m_mul3;
	cout << endl;
	cout << "v_null * 2:" << endl;
	vector<double> v_mul1 = v_null * 2.0;
	cout << "v1 * 2:" << endl;
	vector<double> v_mul2 = v1 * 2.0;
	cout << v_mul2;
	cout << "2 * v1:" << endl;
	vector<double> v_mul3 = 2.0 * v1;
	cout << v_mul3;

	cout << endl;

	// /
	cout << "m_null / 2:" << endl;
	Matrix<double> m_div1 = m_null / 2.0;
	cout << "m1 / 2" << endl;
	Matrix<double> m_div2 = m1 / 2.0;
	cout << m_div2;
	cout << endl;
	cout << "v_null / 2:" << endl;
	vector<double> v_div1 = v_null / 2.0;
	cout << "v1 / 2" << endl;
	vector<double> v_div2 = v1 / 2.0;
	cout << v_div2;

	cout << endl;

	// transposition
	cout << "transposition of m_null:" << endl;
	Matrix<double> m_trans1 = m_null.transposition();
	cout << "transposition of m1:" << endl;
	Matrix<double> m_trans2 = m1.transposition();
	cout << m_trans2;
	cout << "transposition of m3:" << endl;
	Matrix<double> m_trans3 = m1.transposition(m3);
	cout << m_trans3;

	cout << endl;

	// conjugation
	cout << "conjugation of m_null:" << endl;
	Matrix<double> m_conj1 = m_null.conjugation();
	cout << "conjugation of m1:" << endl;
	Matrix<double> m_conj2 = m1.conjugation();
	cout << m_conj2;
	cout << "conjugation of m3:" << endl;
	Matrix<double> m_conj3 = m1.conjugation(m3);
	cout << m_conj3;

	cout << endl;

	// element_wise_mul
	cout << "element_wise_mul of m1 and m_mul:" << endl;
	Matrix<double> m_element_wise_mul1 = m1.element_wise_mul(m_null);
	cout << "element_wise_mul of m1 and m3:" << endl;
	Matrix<double> m_element_wise_mul2 = m1.element_wise_mul(m3);
	cout << "element_wise_mul of m1 and m2:" << endl;
	Matrix<double> m_element_wise_mul3 = m1.element_wise_mul(m1, m2);
	cout << m_element_wise_mul3;
	cout << endl;
	cout << "element_wise_mul of v1 and v_null:" << endl;
	vector<double> v_element_wise_mul1 = m1.element_wise_mul(v1, v_null);
	cout << "element_wise_mul of v1 and v3:" << endl;
	vector<double> v_element_wise_mul2 = m1.element_wise_mul(v1, v3);
	cout << "element_wise_mul of v1 and v2:" << endl;
	vector<double> v_element_wise_mul3 = m1.element_wise_mul(v1, v2);
	cout << v_element_wise_mul3;

	cout << endl;

	// *
	cout << "m1 * m_null:" << endl;
	Matrix<double> m_mulm1 = m1 * m_null;
	cout << "m1 * m3:" << endl;
	Matrix<double> m_mulm2 = m1 * m3;
	cout << "m1 * m2:" << endl;
	Matrix<double> m_mulm3 = m1 * m2;
	cout << m_mulm3;
	cout << endl;
	cout << "v1 * v_null:" << endl;
	vector<double> m_mulv1 = m1 * v_null;
	cout << "v1 * v3:" << endl;
	vector<double> m_mulv2 = m1 * v3;
	cout << "v1 * v2:" << endl;
	vector<double> m_mulv3 = m1 * v2;
	cout << m_mulv3;

	cout << endl;

	// dot_product
	cout << "m1 dot_product m_null:" << endl;
	double m_dot_product1 = m1.dot_product(m_null);
	cout << "m1 dot_product m3:" << endl;
	double m_dot_product2 = m1.dot_product(m3);
	cout << "m1 dot_product m2:" << endl;
	double m_dot_product3 = m1.dot_product(m1, m2);
	cout << m_dot_product3 << endl;
	cout << endl;
	cout << "v1 dot_product v_null:" << endl;
	double v_dot_product1 = m1.dot_product(v1, v_null);
	cout << "v1 dot_product v3:" << endl;
	double v_dot_product2 = m1.dot_product(v1, v3);
	cout << "v1 dot_product v2:" << endl;
	double v_dot_product3 = m1.dot_product(v1, v2);
	cout << v_dot_product3 << endl;


	cout << endl;

	// cross_product
	cout << "v1 cross_product v3:" << endl;
	vector<double> v_cross_product1 = m1.cross_product(v1, v3);
	cout << "v1 cross_product v2:" << endl;
	vector<double> v_cross_product2 = m1.cross_product(v1, v2);
	cout << v_cross_product2;

	cout << endl;

	// check_type_axis
	cout << "check_type_axis of m_null:" << endl;
	bool check_type_axis1 = m_null.check_type_axis(0, 0);
	cout << "check_type_axis of m1:" << endl;
	bool check_type_axis2 = m1.check_type_axis(0, n + 1);
	bool check_type_axis3 = m1.check_type_axis(0, n - 1);
	cout << "valid axis of m1:" << n - 1 << endl;

	cout << endl;

	// find_max
	cout << "find_max of m_null:" << endl;
	double m_max1 = m_null.find_max();
	cout << "find_max of m1:" << endl;
	double m_max2 = m1.find_max();
	cout << m_max2 << endl;
	cout << "find_max of m1 row1:" << endl;
	double m_max3 = m1.find_max(0, 0);
	cout << m_max3 << endl;
	cout << "find_max of m1 col1:" << endl;
	double m_max4 = m1.find_max(1, 0);
	cout << m_max4 << endl;

	cout << endl;

	// find_min
	cout << "find_min of m_null:" << endl;
	double m_min1 = m_null.find_min();
	cout << "find_min of m1:" << endl;
	double m_min2 = m1.find_min();
	cout << m_min2 << endl;
	cout << "find_min of m1 row2:" << endl;
	double m_min3 = m1.find_min(0, 1);
	cout << m_min3 << endl;
	cout << "find_min of m1 col2:" << endl;
	double m_min4 = m1.find_min(1, 1);
	cout << m_min4 << endl;

	cout << endl;

	// sum
	cout << "sum of m_null:" << endl;
	double m_sum1 = m_null.sum();
	cout << "sum of m1:" << endl;
	double m_sum2 = m1.sum();
	cout << m_sum2 << endl;
	cout << "sum of m1 row1:" << endl;
	double m_sum3 = m1.sum(0, 0);
	cout << m_sum3 << endl;
	cout << "sum of m1 col1:" << endl;
	double m_sum4 = m1.sum(1, 0);
	cout << m_sum4 << endl;

	cout << endl;

	// avg
	cout << "avg of m_null:" << endl;
	double m_avg1 = m_null.avg();
	cout << "avg of m1:" << endl;
	double m_avg2 = m1.avg();
	cout << m_avg2 << endl;
	cout << "avg of m1 row1:" << endl;
	double m_avg3 = m1.avg(0, 0);
	cout << m_avg3 << endl;
	cout << "avg of m1 col1:" << endl;
	double m_avg4 = m1.avg(1, 0);
	cout << m_avg4 << endl;

	cout << endl;

	// check_square
	cout << "check_square of m_null:" << endl;
	bool check_square1 = m_null.check_square();
	cout << "check_square of m4:" << endl;
	bool check_square2 = m4.check_square();
	cout << "check_square of m1:" << endl;
	bool check_square3 = m1.check_square();
	cout << check_square3 << endl;

	cout << endl;

	// trace
	cout << "trace of m_null:" << endl;
	double trace1 = m_null.trace();
	cout << "trace of m4:" << endl;
	double trace2 = m4.trace();
	cout << "trace of m1:" << endl;
	double trace3 = m1.trace();
	cout << trace3 << endl;

	cout << endl;

	//determinant
	cout << "determinant of m1:" << endl;
	cout << m1.determinant() << endl;

	cout << endl;

	// inverse
	cout << "inverse of m1:" << endl;
	Matrix<double> inverse1 = m1.inverse();

	// norm2
	cout << "norm2 of m_null:" << endl;
	double norm21 = m_null.norm2();
	cout << "norm2 of m1:" << endl;
	double norm22 = m1.norm2();
	cout << norm22 << endl;

	cout << endl;

	// eigenvalues
	cout << "eigenvalues of m3:" << endl;
	vector<double> eigenvalues1 = m3.eigenvalues();
	cout << eigenvalues1 << endl;

	cout << endl;

	// print_eigenvectors
	cout << "print_eigenvectors of m3:" << endl;
	m3.print_eigenvectors();

	cout << endl;

	// reshape
	cout << "reshape of m1:" << endl;
	Matrix<double> reshape1 = m1.reshape(1, m1.getRow() * m1.getCol());
	cout << reshape1;

	cout << endl;

	// slicing
	cout << "testcase of slicing:" << endl;
	vector<vector<double>> s(4, vector<double>(4));
	s = { {1,2,3,4},{5,6,7,8},{1,3,5,7},{2,4,6,8} };
	Matrix<double> S(s);
	cout << S;
	cout << Matrix<double>::slicing(S, 0, 3, 1, 2, 1, 0);
	cout << S.slicing(2, 3, 0, 3, 0, 0);

	cout << endl;

	// convolution
	cout << "testcase of convolution:" << endl;
	vector<vector<double>> a(5, vector<double>(5));
	a = { {17,24,1,8,15},{23,5,7,14,16},{4,6,13,20,22},{10,12,19,21,3},{11,18,25,2,9} };
	Matrix<double> A(a);
	vector<vector<double>> b(3, vector<double>(3));
	b = { {1,3,1},{0,5,0},{2,1,2} };
	Matrix<double> B(b);
	cout << Matrix<double>::convolution(A, B);
	cout << endl;

	//complex number
	vector<vector<complex<double>>> p(2, vector<complex<double>>(2));
	complex<double> c(1, 1);
	p = { {1,c},{c,2} };
	Matrix<complex<double>> P(p);

	//Matrix<complex<T>>
	complex<double> c1(1, 1);
	complex<double> c2(0, 1);
	vector<vector<complex<double>>> cm;
	vector<complex<double>> r1, r2;
	r1.emplace_back(c1);
	r1.emplace_back(c2);
	r2.emplace_back(c1);
	r2.emplace_back(c1);
	// <<
	cout << "The first complex c1:" << endl;
	cout << c1 << endl;
	cout << "The second complex c2:" << endl;
	cout << c2 << endl;
	cout << endl;
	cout << "The first complex vector r1:" << endl;
	cout << r1 << endl;
	cout << "The second complex vector r2:" << endl;
	cout << r2 << endl;
	cout << endl;

	//vector +
	cout << "r1 + r2:" << endl;
	cout << r1 + r2 << endl;
	//vector -
	cout << "r1 - r2:" << endl;
	cout << r1 - r2 << endl;
	//vector *
	cout << "r1 * r2:" << endl;
	cout << r1 * r2 << endl;
	//vector / val
	cout << "r1 / c2:" << endl;
	cout << r1 / c2 << endl;
	cout << endl;
	cm.emplace_back(r1);
	cm.emplace_back(r2);
	Matrix<complex<double>> C(cm);
	//print C
	cout << "The complex matrix C:" << endl;
	cout << C << endl;
	//matrix calculation
	cout << "The calculation C + 2 * c1 * C:" << endl;
	cout << C + (complex<double>)2 * c1 * C << endl;
	//conjugation
	cout << "The conjugation of complex matrix C:" << endl;
	cout << C.conjugation();
	cout << endl;
	//norm2
	cout << "The norm of complex matrix C:" << endl;	
	cout << C.norm2() << endl;
	cout << "The norm of complex vector r1:" << endl;
	cout << Matrix<complex<double>>::norm2(r1) << endl;
	cout << endl;

	return 0;
}