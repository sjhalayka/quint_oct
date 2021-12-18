#include <cmath>
#include <vector>
#include <iostream>	
using namespace std;


template<class T, size_t N>
class vertex
{
public:
	vertex(void)
	{

	}

	bool operator==(const vertex& rhs)
	{
		bool all_equal = true;

		for (size_t i = 0; i < N; i++)
		{
			T f = fabs(vertex_data[i] - rhs.vertex_data[i]);

			if (f > 0.0001)
			{
				all_equal = false;
				break;
			}
		}

		return all_equal;
	}

	bool operator!=(const vertex& rhs)
	{
		return !(*this == rhs);
	}

	T magnitude(void)
	{
		T all_self_dot = 0;

		for (size_t i = 0; i < N; i++)
			all_self_dot += (vertex_data[i] * vertex_data[i]);

		return sqrt(all_self_dot);
	}

	vertex operator+(const vertex& right) const
	{
		vertex out;

		for (size_t i = 0; i < N; i++)
			out.vertex_data[i] = vertex_data[i] + right.vertex_data[i];

		return out;
	}

	void normalize(void)
	{
		T mag = magnitude();

		if (mag != 0)
		{
			for (size_t i = 0; i < N; i++)
				vertex_data[i] = vertex_data[i] / mag;
		}
	}

	float vertex_data[N];
};


template<class T, size_t N>
vertex<T, N> pow(const vertex<T, N>& in, T beta)
{
	T all_self_dot = 0;
	T imag_self_dot = 0;
	vertex<T, N> out;

	for (size_t i = 0; i < N; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < N; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < N; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const T all_len = sqrt(all_self_dot);
	const T imag_len = sqrt(imag_self_dot);
	const T self_dot_beta = pow(all_self_dot, beta / 2.0f);

	out.vertex_data[0] = self_dot_beta * cos(beta * acos(in.vertex_data[0] / all_len));

	if (imag_len != 0)
	{
		for (size_t i = 1; i < N; i++)
			out.vertex_data[i] = in.vertex_data[i] * self_dot_beta * sin(beta * acos(in.vertex_data[0] / all_len)) / imag_len;
	}

	return out;
}


template<class T, size_t N = 8>
vertex<T, 8> traditional_mul(const vertex<T, 8>& in_a, const vertex<T, 8>& in_b)
{
	vertex<T, 8> out;

	out.vertex_data[0] = in_a.vertex_data[0] * in_b.vertex_data[0] - in_a.vertex_data[1] * in_b.vertex_data[1] - in_a.vertex_data[2] * in_b.vertex_data[2] - in_a.vertex_data[3] * in_b.vertex_data[3] - in_a.vertex_data[4] * in_b.vertex_data[4] - in_a.vertex_data[5] * in_b.vertex_data[5] - in_a.vertex_data[6] * in_b.vertex_data[6] - in_a.vertex_data[7] * in_b.vertex_data[7];
	out.vertex_data[1] = in_a.vertex_data[0] * in_b.vertex_data[1] + in_a.vertex_data[1] * in_b.vertex_data[0] + in_a.vertex_data[2] * in_b.vertex_data[3] - in_a.vertex_data[3] * in_b.vertex_data[2] + in_a.vertex_data[4] * in_b.vertex_data[5] - in_a.vertex_data[5] * in_b.vertex_data[4] - in_a.vertex_data[6] * in_b.vertex_data[7] + in_a.vertex_data[7] * in_b.vertex_data[6];
	out.vertex_data[2] = in_a.vertex_data[0] * in_b.vertex_data[2] - in_a.vertex_data[1] * in_b.vertex_data[3] + in_a.vertex_data[2] * in_b.vertex_data[0] + in_a.vertex_data[3] * in_b.vertex_data[1] + in_a.vertex_data[4] * in_b.vertex_data[6] + in_a.vertex_data[5] * in_b.vertex_data[7] - in_a.vertex_data[6] * in_b.vertex_data[4] - in_a.vertex_data[7] * in_b.vertex_data[5];
	out.vertex_data[3] = in_a.vertex_data[0] * in_b.vertex_data[3] + in_a.vertex_data[1] * in_b.vertex_data[2] - in_a.vertex_data[2] * in_b.vertex_data[1] + in_a.vertex_data[3] * in_b.vertex_data[0] + in_a.vertex_data[4] * in_b.vertex_data[7] - in_a.vertex_data[5] * in_b.vertex_data[6] + in_a.vertex_data[6] * in_b.vertex_data[5] - in_a.vertex_data[7] * in_b.vertex_data[4];
	out.vertex_data[4] = in_a.vertex_data[0] * in_b.vertex_data[4] - in_a.vertex_data[1] * in_b.vertex_data[5] - in_a.vertex_data[2] * in_b.vertex_data[6] - in_a.vertex_data[3] * in_b.vertex_data[7] + in_a.vertex_data[4] * in_b.vertex_data[0] + in_a.vertex_data[5] * in_b.vertex_data[1] + in_a.vertex_data[6] * in_b.vertex_data[2] + in_a.vertex_data[7] * in_b.vertex_data[3];
	out.vertex_data[5] = in_a.vertex_data[0] * in_b.vertex_data[5] + in_a.vertex_data[1] * in_b.vertex_data[4] - in_a.vertex_data[2] * in_b.vertex_data[7] + in_a.vertex_data[3] * in_b.vertex_data[6] - in_a.vertex_data[4] * in_b.vertex_data[1] + in_a.vertex_data[5] * in_b.vertex_data[0] - in_a.vertex_data[6] * in_b.vertex_data[3] + in_a.vertex_data[7] * in_b.vertex_data[2];
	out.vertex_data[6] = in_a.vertex_data[0] * in_b.vertex_data[6] + in_a.vertex_data[1] * in_b.vertex_data[7] + in_a.vertex_data[2] * in_b.vertex_data[4] - in_a.vertex_data[3] * in_b.vertex_data[5] - in_a.vertex_data[4] * in_b.vertex_data[2] + in_a.vertex_data[5] * in_b.vertex_data[3] + in_a.vertex_data[6] * in_b.vertex_data[0] - in_a.vertex_data[7] * in_b.vertex_data[1];
	out.vertex_data[7] = in_a.vertex_data[0] * in_b.vertex_data[7] - in_a.vertex_data[1] * in_b.vertex_data[6] + in_a.vertex_data[2] * in_b.vertex_data[5] + in_a.vertex_data[3] * in_b.vertex_data[4] - in_a.vertex_data[4] * in_b.vertex_data[3] - in_a.vertex_data[5] * in_b.vertex_data[2] + in_a.vertex_data[6] * in_b.vertex_data[1] + in_a.vertex_data[7] * in_b.vertex_data[0];

	return out;
}


template<class T, size_t N>
vertex<T, N> exp(const vertex<T, N>& in)
{
	T all_self_dot = 0;
	T imag_self_dot = 0;
	vertex<T, N> out;

	for (size_t i = 0; i < N; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < N; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < N; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const T l_d = sqrt(all_self_dot);
	const T l_e = sqrt(imag_self_dot);

	out.vertex_data[0] = exp(in.vertex_data[0]) * cos(l_e);

	if (l_e != 0)
	{
		for (size_t i = 1; i < N; i++)
			out.vertex_data[i] = in.vertex_data[i] / l_e * exp(in.vertex_data[0]) * sin(l_e);
	}

	return out;
}

template<class T, size_t N>
vertex<T, N> ln(const vertex<T, N>& in)
{
	T all_self_dot = 0;
	T imag_self_dot = 0;
	vertex<T, N> out;

	for (size_t i = 0; i < N; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < N; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < N; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const T l_d = sqrt(all_self_dot);
	const T l_e = sqrt(imag_self_dot);

	if (in.vertex_data[0] != 0)
	{
		out.vertex_data[0] = log(l_d);
	}

	if (l_e != 0)
	{
		for (size_t i = 1; i < N; i++)
			out.vertex_data[i] = in.vertex_data[i] / l_e * acos(in.vertex_data[0] / l_d);
	}

	return out;
}

template<class T, size_t N>
vertex<T, N> mul(const vertex<T, N>& in_a, const vertex<T, N>& in_b)
{
	// A*B == exp(ln(A) + ln(B))
	return exp(ln(in_a) + ln(in_b));
}

int main(void)
{
	// Compare quintonion pow to mul

	//vertex<float, 5> a;

	//a.vertex_data[0] = 0.1f;
	//a.vertex_data[1] = 0.2f;
	//a.vertex_data[2] = 0.3f;
	//a.vertex_data[3] = 0.4f;
	//a.vertex_data[4] = 0.5f;

	//vertex<float, 5> x = pow(a, 2.0f);

	//vertex<float, 5> y = mul(a, a);

	//cout << x.vertex_data[0] << " " << x.vertex_data[1] << " " << x.vertex_data[2] << " " << x.vertex_data[3] << " " << x.vertex_data[4] << endl;
	//cout << y.vertex_data[0] << " " << y.vertex_data[1] << " " << y.vertex_data[2] << " " << y.vertex_data[3] << " " << y.vertex_data[4] << endl;

	//return 0;



	// Test quintonions for various attributes

	//vertex<float, 5> a;
	//a.vertex_data[0] = 0.1f;
	//a.vertex_data[1] = 0.2f;
	//a.vertex_data[2] = 0.3f;
	//a.vertex_data[3] = 0.4f;
	//a.vertex_data[4] = 0.5f;

	//vertex<float, 5> b;
	//b.vertex_data[0] = 1.0f;
	//b.vertex_data[1] = 0.9f;
	//b.vertex_data[2] = 0.8f;
	//b.vertex_data[3] = 0.7f;
	//b.vertex_data[4] = 0.6f;

	//vertex<float, 5> c;
	//c.vertex_data[0] = 10.0f;
	//c.vertex_data[1] = 9.0f;
	//c.vertex_data[2] = 8.0f;
	//c.vertex_data[3] = 7.0f;
	//c.vertex_data[4] = 6.0f;

	//vertex<float, 5> x = mul(a, b);
	//vertex<float, 5> y = mul(b, a);

	//if (x != y)
	//{		
	//	cout << "commutativity failure" << endl;

	//	cout << x.vertex_data[0] << " " << x.vertex_data[1] << " " << x.vertex_data[2] << " " << x.vertex_data[3] << " " << x.vertex_data[4] << endl;
	//	cout << y.vertex_data[0] << " " << y.vertex_data[1] << " " << y.vertex_data[2] << " " << y.vertex_data[3] << " " << y.vertex_data[4] << endl;
	//}

	//x = mul(mul(a, b), c);
	//y = mul(a, mul(b, c));

	//if (x != y)
	//{
	//	cout << "associativity failure" << endl;

	//	cout << x.vertex_data[0] << " " << x.vertex_data[1] << " " << x.vertex_data[2] << " " << x.vertex_data[3] << " " << x.vertex_data[4] << endl;
	//	cout << y.vertex_data[0] << " " << y.vertex_data[1] << " " << y.vertex_data[2] << " " << y.vertex_data[3] << " " << y.vertex_data[4] << endl;
	//}

	//x = mul(a, b + c);
	//y = mul(a, b) + mul(a, c);

	//if (x != y)
	//{
	//	cout << "distributive failure" << endl;

	//	cout << x.vertex_data[0] << " " << x.vertex_data[1] << " " << x.vertex_data[2] << " " << x.vertex_data[3] << " " << x.vertex_data[4] << endl;
	//	cout << y.vertex_data[0] << " " << y.vertex_data[1] << " " << y.vertex_data[2] << " " << y.vertex_data[3] << " " << y.vertex_data[4] << endl;
	//}

	//return 0;



	// Test octonion for various attributes

	//vertex<float, 8> a;
	//a.vertex_data[0] = 0.1f;
	//a.vertex_data[1] = 0.2f;
	//a.vertex_data[2] = 0.3f;
	//a.vertex_data[3] = 0.4f;
	//a.vertex_data[4] = 0.5f;
	//a.vertex_data[5] = 0.6f;
	//a.vertex_data[6] = 0.7f;
	//a.vertex_data[7] = 0.8f;

	//vertex<float, 8> b;
	//b.vertex_data[0] = 1.0f;
	//b.vertex_data[1] = 0.9f;
	//b.vertex_data[2] = 0.8f;
	//b.vertex_data[3] = 0.7f;
	//b.vertex_data[4] = 0.6f;
	//b.vertex_data[5] = 0.5f;
	//b.vertex_data[6] = 0.4f;
	//b.vertex_data[7] = 0.3f;

	//vertex<float, 8> c;
	//c.vertex_data[0] = 10.0f;
	//c.vertex_data[1] = 9.0f;
	//c.vertex_data[2] = 8.0f;
	//c.vertex_data[3] = 7.0f;
	//c.vertex_data[4] = 6.0f;
	//c.vertex_data[5] = 5.0f;
	//c.vertex_data[6] = 4.0f;
	//c.vertex_data[7] = 3.0f;

	//vertex<float, 8> x = mul(a, b);
	//vertex<float, 8> y = mul(b, a);

	//if (x != y)
	//	cout << "commutativity failure" << endl;

	//x = mul(mul(a, b), c);
	//y = mul(a, mul(b, c));

	//if (x != y)
	//	cout << "associativity failure" << endl;

	//x = mul(a, b + c);
	//y = mul(a, b) + mul(a, c);

	//if (x != y)
	//	cout << "distributive failure" << endl;

	//return 0;



	// Test octonion multiplication where A != B

	vertex<float, 8> a;
	a.vertex_data[0] = 0.1f;
	a.vertex_data[1] = 0.2f;
	a.vertex_data[2] = 0.3f;
	a.vertex_data[3] = 0.4f;
	a.vertex_data[4] = 0.5f;
	a.vertex_data[5] = 0.6f;
	a.vertex_data[6] = 0.7f;
	a.vertex_data[7] = 0.8f;

	vertex<float, 8> b;
	b.vertex_data[0] = 10.0f;
	b.vertex_data[1] = 9.0f;
	b.vertex_data[2] = 8.0f;
	b.vertex_data[3] = 7.0f;
	b.vertex_data[4] = 6.0f;
	b.vertex_data[5] = 5.0f;
	b.vertex_data[6] = 4.0f;
	b.vertex_data[7] = 3.0f;

	vertex<float, 8> P = traditional_mul(a, b);
	vertex<float, 8> P2 = mul(a, b);

	for (size_t i = 0; i < 8; i++)
		cout << P.vertex_data[i] << " ";

	cout << endl;

	for (size_t i = 0; i < 8; i++)
		cout << P2.vertex_data[i] << " ";

	cout << endl;

	cout << P.magnitude() << " " << P2.magnitude() << endl;

	return 0;



	// Test for subalgebra	

	//srand(time(0));

	//for (size_t num_tries = 0; num_tries < 10000; num_tries++)
	//{
	//	float a0 = (rand() % RAND_MAX) / static_cast<float>(RAND_MAX - 1);

	//	if (rand() % 2)
	//		a0 = -a0;

	//	float a1 = (rand() % RAND_MAX) / static_cast<float>(RAND_MAX - 1);

	//	if (rand() % 2)
	//		a1 = -a1;

	//	float a2 = (rand() % RAND_MAX) / static_cast<float>(RAND_MAX - 1);

	//	if (rand() % 2)
	//		a2 = -a2;

	//	float a3 = (rand() % RAND_MAX) / static_cast<float>(RAND_MAX - 1);

	//	if (rand() % 2)
	//		a3 = -a3;

	//	float a4 = (rand() % RAND_MAX) / static_cast<float>(RAND_MAX - 1);

	//	if (rand() % 2)
	//		a4 = -a4;

	//	vertex<float, 8> a;
	//	a.vertex_data[0] = a0;
	//	a.vertex_data[1] = a1;
	//	a.vertex_data[2] = a2;
	//	a.vertex_data[3] = a3;
	//	a.vertex_data[4] = a4;
	//	a.vertex_data[5] = 0;
	//	a.vertex_data[6] = 0;
	//	a.vertex_data[7] = 0;

	//	vertex<float, 8> b = a;

	//	vertex<float, 8> P = traditional_mul(a, b);

	//	if (P.vertex_data[5] || P.vertex_data[6] || P.vertex_data[7])
	//	{
	//		cout << "Error: non-zero components!" << endl;
	//	}
	//}

	return 0;
}