#include <cmath>
#include <vector>
#include <iostream>
#include <complex>
using namespace std;



template<class T, size_t N>
class vertex
{
public:
	vector<T> vd; // vertex data

	vertex(void)
	{
		vd.resize(N, 0);
	}

	bool operator==(const vertex& rhs) const
	{
		bool all_equal = true;

		for (size_t i = 0; i < N; i++)
		{
			T f = fabs(vd[i] - rhs.vd[i]);

			if (f > 0.0001)
			{
				all_equal = false;
				break;
			}
		}

		return all_equal;
	}

	bool operator!=(const vertex& rhs) const
	{
		return !(*this == rhs);
	}

	vertex operator+(const vertex& rhs) const
	{
		vertex out;

		for (size_t i = 0; i < N; i++)
			out.vd[i] = vd[i] + rhs.vd[i];

		return out;
	}

	T magnitude(void) const
	{
		return sqrt(all_dot(*this));
	}

	T all_dot(const vertex& rhs) const
	{
		T all_self_dot = 0;

		for (size_t i = 0; i < N; i++)
			all_self_dot += (vd[i] * rhs.vd[i]);

		return all_self_dot;
	}

	T imag_dot(const vertex& rhs)
	{
		T imag_self_dot = 0;

		for (size_t i = 1; i < N; i++)
			imag_self_dot += (vd[i] * rhs.vd[i]);

		return imag_self_dot;
	}
};

template<class T, size_t N>
vertex<T, N> pow(const vertex<T, N>& in, T beta)
{
	T all_self_dot = 0;
	T imag_self_dot = 0;
	vertex<T, N> out;

	for (size_t i = 1; i < N; i++)
		imag_self_dot += (in.vd[i] * in.vd[i]);

	all_self_dot = imag_self_dot + (in.vd[0] * in.vd[0]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < N; i++)
			out.vd[i] = 0;

		return out;
	}

	const T l_d = sqrt(all_self_dot);
	const T l_e = sqrt(imag_self_dot);
	const T self_dot_beta = pow(all_self_dot, beta / 2.0f);

	out.vd[0] = self_dot_beta * cos(beta * acos(in.vd[0] / l_d));

	if (l_e != 0)
	{
		const T x = 1 / l_e;
		const T y = self_dot_beta * sin(beta * acos(in.vd[0] / l_d));
		const T z = x * y;

		for (size_t i = 1; i < N; i++)
			out.vd[i] = in.vd[i] * z;
	}

	return out;
}

template<class T, size_t N = 1>
vertex<T, 1> traditional_mul(const vertex<T, 1>& in_a, const vertex<T, 1>& in_b)
{
	vertex<T, 1> out;

	out.vd[0] = in_a.vd[0] * in_b.vd[0];

	return out;
}

template<class T, size_t N = 2>
vertex<T, 2> traditional_mul(const vertex<T, 2>& in_a, const vertex<T, 2>& in_b)
{
	vertex<T, 2> out;

	out.vd[0] = in_a.vd[0] * in_b.vd[0] - in_a.vd[1] * in_b.vd[1];
	out.vd[1] = in_a.vd[0] * in_b.vd[1] + in_a.vd[1] * in_b.vd[0];

	return out;
}

template<class T, size_t N = 4>
vertex<T, 4> traditional_mul(const vertex<T, 4>& in_a, const vertex<T, 4>& in_b)
{
	vertex<T, 4> out;

	out.vd[0] = in_a.vd[0] * in_b.vd[0] - in_a.vd[1] * in_b.vd[1] - in_a.vd[2] * in_b.vd[2] - in_a.vd[3] * in_b.vd[3];
	out.vd[1] = in_a.vd[0] * in_b.vd[1] + in_a.vd[1] * in_b.vd[0] + in_a.vd[2] * in_b.vd[3] - in_a.vd[3] * in_b.vd[2];
	out.vd[2] = in_a.vd[0] * in_b.vd[2] - in_a.vd[1] * in_b.vd[3] + in_a.vd[2] * in_b.vd[0] + in_a.vd[3] * in_b.vd[1];
	out.vd[3] = in_a.vd[0] * in_b.vd[3] + in_a.vd[1] * in_b.vd[2] - in_a.vd[2] * in_b.vd[1] + in_a.vd[3] * in_b.vd[0];

	return out;
}

template<class T, size_t N = 8>
vertex<T, 8> traditional_mul(const vertex<T, 8>& in_a, const vertex<T, 8>& in_b)
{
	vertex<T, 8> out;

	out.vd[0] = in_a.vd[0] * in_b.vd[0] - in_a.vd[1] * in_b.vd[1] - in_a.vd[2] * in_b.vd[2] - in_a.vd[3] * in_b.vd[3] - in_a.vd[4] * in_b.vd[4] - in_a.vd[5] * in_b.vd[5] - in_a.vd[6] * in_b.vd[6] - in_a.vd[7] * in_b.vd[7];
	out.vd[1] = in_a.vd[0] * in_b.vd[1] + in_a.vd[1] * in_b.vd[0] + in_a.vd[2] * in_b.vd[3] - in_a.vd[3] * in_b.vd[2] + in_a.vd[4] * in_b.vd[5] - in_a.vd[5] * in_b.vd[4] - in_a.vd[6] * in_b.vd[7] + in_a.vd[7] * in_b.vd[6];
	out.vd[2] = in_a.vd[0] * in_b.vd[2] - in_a.vd[1] * in_b.vd[3] + in_a.vd[2] * in_b.vd[0] + in_a.vd[3] * in_b.vd[1] + in_a.vd[4] * in_b.vd[6] + in_a.vd[5] * in_b.vd[7] - in_a.vd[6] * in_b.vd[4] - in_a.vd[7] * in_b.vd[5];
	out.vd[3] = in_a.vd[0] * in_b.vd[3] + in_a.vd[1] * in_b.vd[2] - in_a.vd[2] * in_b.vd[1] + in_a.vd[3] * in_b.vd[0] + in_a.vd[4] * in_b.vd[7] - in_a.vd[5] * in_b.vd[6] + in_a.vd[6] * in_b.vd[5] - in_a.vd[7] * in_b.vd[4];
	out.vd[4] = in_a.vd[0] * in_b.vd[4] - in_a.vd[1] * in_b.vd[5] - in_a.vd[2] * in_b.vd[6] - in_a.vd[3] * in_b.vd[7] + in_a.vd[4] * in_b.vd[0] + in_a.vd[5] * in_b.vd[1] + in_a.vd[6] * in_b.vd[2] + in_a.vd[7] * in_b.vd[3];
	out.vd[5] = in_a.vd[0] * in_b.vd[5] + in_a.vd[1] * in_b.vd[4] - in_a.vd[2] * in_b.vd[7] + in_a.vd[3] * in_b.vd[6] - in_a.vd[4] * in_b.vd[1] + in_a.vd[5] * in_b.vd[0] - in_a.vd[6] * in_b.vd[3] + in_a.vd[7] * in_b.vd[2];
	out.vd[6] = in_a.vd[0] * in_b.vd[6] + in_a.vd[1] * in_b.vd[7] + in_a.vd[2] * in_b.vd[4] - in_a.vd[3] * in_b.vd[5] - in_a.vd[4] * in_b.vd[2] + in_a.vd[5] * in_b.vd[3] + in_a.vd[6] * in_b.vd[0] - in_a.vd[7] * in_b.vd[1];
	out.vd[7] = in_a.vd[0] * in_b.vd[7] - in_a.vd[1] * in_b.vd[6] + in_a.vd[2] * in_b.vd[5] + in_a.vd[3] * in_b.vd[4] - in_a.vd[4] * in_b.vd[3] - in_a.vd[5] * in_b.vd[2] + in_a.vd[6] * in_b.vd[1] + in_a.vd[7] * in_b.vd[0];

	return out;
}

template<class T, size_t N = 16>
vertex<T, 16> traditional_mul(const vertex<T, 16>& in_a, const vertex<T, 16>& in_b)
{
	vertex<T, 16> out;

	out.vd[0] = in_a.vd[0] * in_b.vd[0] - in_a.vd[1] * in_b.vd[1] - in_a.vd[2] * in_b.vd[2] - in_a.vd[3] * in_b.vd[3] - in_a.vd[4] * in_b.vd[4] - in_a.vd[5] * in_b.vd[5] - in_a.vd[6] * in_b.vd[6] - in_a.vd[7] * in_b.vd[7] - in_a.vd[8] * in_b.vd[8] - in_a.vd[9] * in_b.vd[9] - in_a.vd[10] * in_b.vd[10] - in_a.vd[11] * in_b.vd[11] - in_a.vd[12] * in_b.vd[12] - in_a.vd[13] * in_b.vd[13] - in_a.vd[14] * in_b.vd[14] - in_a.vd[15] * in_b.vd[15];
	out.vd[1] = in_a.vd[0] * in_b.vd[1] + in_a.vd[1] * in_b.vd[0] + in_a.vd[2] * in_b.vd[3] - in_a.vd[3] * in_b.vd[2] + in_a.vd[4] * in_b.vd[5] - in_a.vd[5] * in_b.vd[4] - in_a.vd[6] * in_b.vd[7] + in_a.vd[7] * in_b.vd[6] + in_a.vd[8] * in_b.vd[9] - in_a.vd[9] * in_b.vd[8] - in_a.vd[10] * in_b.vd[11] + in_a.vd[11] * in_b.vd[10] - in_a.vd[12] * in_b.vd[13] + in_a.vd[13] * in_b.vd[12] + in_a.vd[14] * in_b.vd[15] - in_a.vd[15] * in_b.vd[14];
	out.vd[2] = in_a.vd[0] * in_b.vd[2] - in_a.vd[1] * in_b.vd[3] + in_a.vd[2] * in_b.vd[0] + in_a.vd[3] * in_b.vd[1] + in_a.vd[4] * in_b.vd[6] + in_a.vd[5] * in_b.vd[7] - in_a.vd[6] * in_b.vd[4] - in_a.vd[7] * in_b.vd[5] + in_a.vd[8] * in_b.vd[10] + in_a.vd[9] * in_b.vd[11] - in_a.vd[10] * in_b.vd[8] - in_a.vd[11] * in_b.vd[9] - in_a.vd[12] * in_b.vd[14] - in_a.vd[13] * in_b.vd[15] + in_a.vd[14] * in_b.vd[12] + in_a.vd[15] * in_b.vd[13];
	out.vd[3] = in_a.vd[0] * in_b.vd[3] + in_a.vd[1] * in_b.vd[2] - in_a.vd[2] * in_b.vd[1] + in_a.vd[3] * in_b.vd[0] + in_a.vd[4] * in_b.vd[7] - in_a.vd[5] * in_b.vd[6] + in_a.vd[6] * in_b.vd[5] - in_a.vd[7] * in_b.vd[4] + in_a.vd[8] * in_b.vd[11] - in_a.vd[9] * in_b.vd[10] + in_a.vd[10] * in_b.vd[9] - in_a.vd[11] * in_b.vd[8] - in_a.vd[12] * in_b.vd[15] + in_a.vd[13] * in_b.vd[14] - in_a.vd[14] * in_b.vd[13] + in_a.vd[15] * in_b.vd[12];
	out.vd[4] = in_a.vd[0] * in_b.vd[4] - in_a.vd[1] * in_b.vd[5] - in_a.vd[2] * in_b.vd[6] - in_a.vd[3] * in_b.vd[7] + in_a.vd[4] * in_b.vd[0] + in_a.vd[5] * in_b.vd[1] + in_a.vd[6] * in_b.vd[2] + in_a.vd[7] * in_b.vd[3] + in_a.vd[8] * in_b.vd[12] + in_a.vd[9] * in_b.vd[13] + in_a.vd[10] * in_b.vd[14] + in_a.vd[11] * in_b.vd[15] - in_a.vd[12] * in_b.vd[8] - in_a.vd[13] * in_b.vd[9] - in_a.vd[14] * in_b.vd[10] - in_a.vd[15] * in_b.vd[11];
	out.vd[5] = in_a.vd[0] * in_b.vd[5] + in_a.vd[1] * in_b.vd[4] - in_a.vd[2] * in_b.vd[7] + in_a.vd[3] * in_b.vd[6] - in_a.vd[4] * in_b.vd[1] + in_a.vd[5] * in_b.vd[0] - in_a.vd[6] * in_b.vd[3] + in_a.vd[7] * in_b.vd[2] + in_a.vd[8] * in_b.vd[13] - in_a.vd[9] * in_b.vd[12] + in_a.vd[10] * in_b.vd[15] - in_a.vd[11] * in_b.vd[14] + in_a.vd[12] * in_b.vd[9] - in_a.vd[13] * in_b.vd[8] + in_a.vd[14] * in_b.vd[11] - in_a.vd[15] * in_b.vd[10];
	out.vd[6] = in_a.vd[0] * in_b.vd[6] + in_a.vd[1] * in_b.vd[7] + in_a.vd[2] * in_b.vd[4] - in_a.vd[3] * in_b.vd[5] - in_a.vd[4] * in_b.vd[2] + in_a.vd[5] * in_b.vd[3] + in_a.vd[6] * in_b.vd[0] - in_a.vd[7] * in_b.vd[1] + in_a.vd[8] * in_b.vd[14] - in_a.vd[9] * in_b.vd[15] - in_a.vd[10] * in_b.vd[12] + in_a.vd[11] * in_b.vd[13] + in_a.vd[12] * in_b.vd[10] - in_a.vd[13] * in_b.vd[11] - in_a.vd[14] * in_b.vd[8] + in_a.vd[15] * in_b.vd[9];
	out.vd[7] = in_a.vd[0] * in_b.vd[7] - in_a.vd[1] * in_b.vd[6] + in_a.vd[2] * in_b.vd[5] + in_a.vd[3] * in_b.vd[4] - in_a.vd[4] * in_b.vd[3] - in_a.vd[5] * in_b.vd[2] + in_a.vd[6] * in_b.vd[1] + in_a.vd[7] * in_b.vd[0] + in_a.vd[8] * in_b.vd[15] + in_a.vd[9] * in_b.vd[14] - in_a.vd[10] * in_b.vd[13] - in_a.vd[11] * in_b.vd[12] + in_a.vd[12] * in_b.vd[11] + in_a.vd[13] * in_b.vd[10] - in_a.vd[14] * in_b.vd[9] - in_a.vd[15] * in_b.vd[8];
	out.vd[8] = in_a.vd[0] * in_b.vd[8] - in_a.vd[1] * in_b.vd[9] - in_a.vd[2] * in_b.vd[10] - in_a.vd[3] * in_b.vd[11] - in_a.vd[4] * in_b.vd[12] - in_a.vd[5] * in_b.vd[13] - in_a.vd[6] * in_b.vd[14] - in_a.vd[7] * in_b.vd[15] + in_a.vd[8] * in_b.vd[0] + in_a.vd[9] * in_b.vd[1] + in_a.vd[10] * in_b.vd[2] + in_a.vd[11] * in_b.vd[3] + in_a.vd[12] * in_b.vd[4] + in_a.vd[13] * in_b.vd[5] + in_a.vd[14] * in_b.vd[6] + in_a.vd[15] * in_b.vd[7];
	out.vd[9] = in_a.vd[0] * in_b.vd[9] + in_a.vd[1] * in_b.vd[8] - in_a.vd[2] * in_b.vd[11] + in_a.vd[3] * in_b.vd[10] - in_a.vd[4] * in_b.vd[13] + in_a.vd[5] * in_b.vd[12] + in_a.vd[6] * in_b.vd[15] - in_a.vd[7] * in_b.vd[14] - in_a.vd[8] * in_b.vd[1] + in_a.vd[9] * in_b.vd[0] - in_a.vd[10] * in_b.vd[3] + in_a.vd[11] * in_b.vd[2] - in_a.vd[12] * in_b.vd[5] + in_a.vd[13] * in_b.vd[4] + in_a.vd[14] * in_b.vd[7] - in_a.vd[15] * in_b.vd[6];
	out.vd[10] = in_a.vd[0] * in_b.vd[10] + in_a.vd[1] * in_b.vd[11] + in_a.vd[2] * in_b.vd[8] - in_a.vd[3] * in_b.vd[9] - in_a.vd[4] * in_b.vd[14] - in_a.vd[5] * in_b.vd[15] + in_a.vd[6] * in_b.vd[12] + in_a.vd[7] * in_b.vd[13] - in_a.vd[8] * in_b.vd[2] + in_a.vd[9] * in_b.vd[3] + in_a.vd[10] * in_b.vd[0] - in_a.vd[11] * in_b.vd[1] - in_a.vd[12] * in_b.vd[6] - in_a.vd[13] * in_b.vd[7] + in_a.vd[14] * in_b.vd[4] + in_a.vd[15] * in_b.vd[5];
	out.vd[11] = in_a.vd[0] * in_b.vd[11] - in_a.vd[1] * in_b.vd[10] + in_a.vd[2] * in_b.vd[9] + in_a.vd[3] * in_b.vd[8] - in_a.vd[4] * in_b.vd[15] + in_a.vd[5] * in_b.vd[14] - in_a.vd[6] * in_b.vd[13] + in_a.vd[7] * in_b.vd[12] - in_a.vd[8] * in_b.vd[3] - in_a.vd[9] * in_b.vd[2] + in_a.vd[10] * in_b.vd[1] + in_a.vd[11] * in_b.vd[0] - in_a.vd[12] * in_b.vd[7] + in_a.vd[13] * in_b.vd[6] - in_a.vd[14] * in_b.vd[5] + in_a.vd[15] * in_b.vd[4];
	out.vd[12] = in_a.vd[0] * in_b.vd[12] + in_a.vd[1] * in_b.vd[13] + in_a.vd[2] * in_b.vd[14] + in_a.vd[3] * in_b.vd[15] + in_a.vd[4] * in_b.vd[8] - in_a.vd[5] * in_b.vd[9] - in_a.vd[6] * in_b.vd[10] - in_a.vd[7] * in_b.vd[11] - in_a.vd[8] * in_b.vd[4] + in_a.vd[9] * in_b.vd[5] + in_a.vd[10] * in_b.vd[6] + in_a.vd[11] * in_b.vd[7] + in_a.vd[12] * in_b.vd[0] - in_a.vd[13] * in_b.vd[1] - in_a.vd[14] * in_b.vd[2] - in_a.vd[15] * in_b.vd[3];
	out.vd[13] = in_a.vd[0] * in_b.vd[13] - in_a.vd[1] * in_b.vd[12] + in_a.vd[2] * in_b.vd[15] - in_a.vd[3] * in_b.vd[14] + in_a.vd[4] * in_b.vd[9] + in_a.vd[5] * in_b.vd[8] + in_a.vd[6] * in_b.vd[11] - in_a.vd[7] * in_b.vd[10] - in_a.vd[8] * in_b.vd[5] - in_a.vd[9] * in_b.vd[4] + in_a.vd[10] * in_b.vd[7] - in_a.vd[11] * in_b.vd[6] + in_a.vd[12] * in_b.vd[1] + in_a.vd[13] * in_b.vd[0] + in_a.vd[14] * in_b.vd[3] - in_a.vd[15] * in_b.vd[2];
	out.vd[14] = in_a.vd[0] * in_b.vd[14] - in_a.vd[1] * in_b.vd[15] - in_a.vd[2] * in_b.vd[12] + in_a.vd[3] * in_b.vd[13] + in_a.vd[4] * in_b.vd[10] - in_a.vd[5] * in_b.vd[11] + in_a.vd[6] * in_b.vd[8] + in_a.vd[7] * in_b.vd[9] - in_a.vd[8] * in_b.vd[6] - in_a.vd[9] * in_b.vd[7] - in_a.vd[10] * in_b.vd[4] + in_a.vd[11] * in_b.vd[5] + in_a.vd[12] * in_b.vd[2] - in_a.vd[13] * in_b.vd[3] + in_a.vd[14] * in_b.vd[0] + in_a.vd[15] * in_b.vd[1];
	out.vd[15] = in_a.vd[0] * in_b.vd[15] + in_a.vd[1] * in_b.vd[14] - in_a.vd[2] * in_b.vd[13] - in_a.vd[3] * in_b.vd[12] + in_a.vd[4] * in_b.vd[11] + in_a.vd[5] * in_b.vd[10] - in_a.vd[6] * in_b.vd[9] + in_a.vd[7] * in_b.vd[8] - in_a.vd[8] * in_b.vd[7] + in_a.vd[9] * in_b.vd[6] - in_a.vd[10] * in_b.vd[5] - in_a.vd[11] * in_b.vd[4] + in_a.vd[12] * in_b.vd[3] + in_a.vd[13] * in_b.vd[2] - in_a.vd[14] * in_b.vd[1] + in_a.vd[15] * in_b.vd[0];

	return out;
}



template<class T, size_t N>
vertex<T, N> exp(const vertex<T, N>& in)
{
	T all_self_dot = 0;
	T imag_self_dot = 0;
	vertex<T, N> out;

	for (size_t i = 1; i < N; i++)
		imag_self_dot += (in.vd[i] * in.vd[i]);

	all_self_dot = imag_self_dot + (in.vd[0] * in.vd[0]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < N; i++)
			out.vd[i] = 0;

		return out;
	}

//	const T l_d = sqrt(all_self_dot); // not needed
	const T l_e = sqrt(imag_self_dot);

	out.vd[0] = exp(in.vd[0]) * cos(l_e);

	if (l_e != 0)
	{
		const T x = 1 / l_e;
		const T y = exp(in.vd[0]) * sin(l_e);
		const T z = x * y;

		for (size_t i = 1; i < N; i++)
			out.vd[i] = in.vd[i] * z;
	}

	return out;
}

template<class T, size_t N>
vertex<T, N> log(const vertex<T, N>& in)
{
	T all_self_dot = 0;
	T imag_self_dot = 0;
	vertex<T, N> out;

	for (size_t i = 1; i < N; i++)
		imag_self_dot += (in.vd[i] * in.vd[i]);

	all_self_dot = imag_self_dot + (in.vd[0] * in.vd[0]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < N; i++)
			out.vd[i] = 0;

		return out;
	}

	const T l_d = sqrt(all_self_dot);
	const T l_e = sqrt(imag_self_dot);

	if (l_d != 0)
	{
		out.vd[0] = log(l_d);
	}

	if (l_e != 0)
	{
		const T x = 1 / l_e;
		const T y = acos(in.vd[0] / l_d);
		const T z = x * y;

		for (size_t i = 1; i < N; i++)
			out.vd[i] = in.vd[i] * z;
	}

	return out;
}

template<class T, size_t N>
vertex<T, N> mul(const vertex<T, N>& in_a, const vertex<T, N>& in_b)
{
	return exp(log(in_a) + log(in_b));
}



int main(void)
{
	// Compare real numbers 

	//vertex<float, 1> a;
	//a.vd[0] = 0.1234f;

	//vertex<float, 1> b;
	//b.vd[0] = 4.567f;

	//vertex<float, 1> x = mul(a, b);
	//vertex<float, 1> y = traditional_mul(a, b);

	//cout << x.vd[0] << endl;
	//cout << y.vd[0] << endl;

	//return 0;



	// Compare complex numbers 

	//vertex<float, 2> a;
	//a.vd[0] = 0.1f;
	//a.vd[1] = 0.2f;

	//vertex<float, 2> b;
	//b.vd[0] = 1.0f;
	//b.vd[1] = 0.9f;

	//vertex<float, 2> x = mul(a, b);
	//vertex<float, 2> y = traditional_mul(a, b);

	//cout << x.vd[0] << " " << x.vd[1] << endl;
	//cout << y.vd[0] << " " << y.vd[1] << endl;

	//complex<float> cf_a(a.vd[0], a.vd[1]);
	//complex<float> cf_b(b.vd[0], b.vd[1]);

	//complex<float> cf_x = cf_a * cf_b;

	//cout << cf_x.real() << " " << cf_x.imag() << endl;

	//return 0;



	// Compare quaternion numbers 

	//vertex<float, 4> a;
	//a.vd[0] = 0.1f;
	//a.vd[1] = 0.2f;
	//a.vd[2] = 0.3f;
	//a.vd[3] = 0.4f;

	//vertex<float, 4> b;
	//b.vd[0] = 1.0f;
	//b.vd[1] = 0.9f;
	//b.vd[2] = 0.8f;
	//b.vd[3] = 0.7f;

	//vertex<float, 4> x = mul(a, b);
	//vertex<float, 4> y = traditional_mul(a, b);

	//cout << x.vd[0] << " " << x.vd[1] << " " << x.vd[2] << " " << x.vd[3] << endl;
	//cout << y.vd[0] << " " << y.vd[1] << " " << y.vd[2] << " " << y.vd[3] << endl;
	//
	//cout << x.magnitude() << " " << y.magnitude() << endl;

	//return 0;



	// Compare quintonion pow to mul

	//vertex<float, 5> a;

	//a.vd[0] = 0.1f;
	//a.vd[1] = 0.2f;
	//a.vd[2] = 0.3f;
	//a.vd[3] = 0.4f;
	//a.vd[4] = 0.5f;

	//vertex<float, 5> x = pow(a, 2.0f);

	//vertex<float, 5> y = mul(a, a);

	//cout << x.vd[0] << " " << x.vd[1] << " " << x.vd[2] << " " << x.vd[3] << " " << x.vd[4] << endl;
	//cout << y.vd[0] << " " << y.vd[1] << " " << y.vd[2] << " " << y.vd[3] << " " << y.vd[4] << endl;

	//return 0;



	// Test quintonions for various attributes

	//vertex<float, 5> a;
	//a.vd[0] = 0.1f;
	//a.vd[1] = 0.2f;
	//a.vd[2] = 0.3f;
	//a.vd[3] = 0.4f;
	//a.vd[4] = 0.5f;

	//vertex<float, 5> b;
	//b.vd[0] = 1.0f;
	//b.vd[1] = 0.9f;
	//b.vd[2] = 0.8f;
	//b.vd[3] = 0.7f;
	//b.vd[4] = 0.6f;

	//vertex<float, 5> c;
	//c.vd[0] = 10.0f;
	//c.vd[1] = 9.0f;
	//c.vd[2] = 8.0f;
	//c.vd[3] = 7.0f;
	//c.vd[4] = 6.0f;

	//vertex<float, 5> x = mul(a, b);
	//vertex<float, 5> y = mul(b, a);

	//if (x != y)
	//{		
	//	cout << "commutativity failure" << endl;

	//	cout << x.vd[0] << " " << x.vd[1] << " " << x.vd[2] << " " << x.vd[3] << " " << x.vd[4] << endl;
	//	cout << y.vd[0] << " " << y.vd[1] << " " << y.vd[2] << " " << y.vd[3] << " " << y.vd[4] << endl;
	//}

	//x = mul(mul(a, b), c);
	//y = mul(a, mul(b, c));

	//if (x != y)
	//{
	//	cout << "associativity failure" << endl;

	//	cout << x.vd[0] << " " << x.vd[1] << " " << x.vd[2] << " " << x.vd[3] << " " << x.vd[4] << endl;
	//	cout << y.vd[0] << " " << y.vd[1] << " " << y.vd[2] << " " << y.vd[3] << " " << y.vd[4] << endl;
	//}

	//x = mul(a, b + c);
	//y = mul(a, b) + mul(a, c);

	//if (x != y)
	//{
	//	cout << "distributive failure" << endl;

	//	cout << x.vd[0] << " " << x.vd[1] << " " << x.vd[2] << " " << x.vd[3] << " " << x.vd[4] << endl;
	//	cout << y.vd[0] << " " << y.vd[1] << " " << y.vd[2] << " " << y.vd[3] << " " << y.vd[4] << endl;
	//}

	//return 0;



	// Test octonion new multiplication for various attributes

	//vertex<float, 8> a;
	//a.vd[0] = 0.1f;
	//a.vd[1] = 0.2f;
	//a.vd[2] = 0.3f;
	//a.vd[3] = 0.4f;
	//a.vd[4] = 0.5f;
	//a.vd[5] = 0.6f;
	//a.vd[6] = 0.7f;
	//a.vd[7] = 0.8f;

	//vertex<float, 8> b;
	//b.vd[0] = 1.0f;
	//b.vd[1] = 0.9f;
	//b.vd[2] = 0.8f;
	//b.vd[3] = 0.7f;
	//b.vd[4] = 0.6f;
	//b.vd[5] = 0.5f;
	//b.vd[6] = 0.4f;
	//b.vd[7] = 0.3f;

	//vertex<float, 8> c;
	//c.vd[0] = 10.0f;
	//c.vd[1] = 9.0f;
	//c.vd[2] = 8.0f;
	//c.vd[3] = 7.0f;
	//c.vd[4] = 6.0f;
	//c.vd[5] = 5.0f;
	//c.vd[6] = 4.0f;
	//c.vd[7] = 3.0f;

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



	// Test octonion traditional multiplication for various attributes

	//vertex<float, 8> a;
	//a.vd[0] = 0.1f;
	//a.vd[1] = 0.2f;
	//a.vd[2] = 0.3f;
	//a.vd[3] = 0.4f;
	//a.vd[4] = 0.5f;
	//a.vd[5] = 0.6f;
	//a.vd[6] = 0.7f;
	//a.vd[7] = 0.8f;

	//vertex<float, 8> b;
	//b.vd[0] = 1.0f;
	//b.vd[1] = 0.9f;
	//b.vd[2] = 0.8f;
	//b.vd[3] = 0.7f;
	//b.vd[4] = 0.6f;
	//b.vd[5] = 0.5f;
	//b.vd[6] = 0.4f;
	//b.vd[7] = 0.3f;

	//vertex<float, 8> c;
	//c.vd[0] = 10.0f;
	//c.vd[1] = 9.0f;
	//c.vd[2] = 8.0f;
	//c.vd[3] = 7.0f;
	//c.vd[4] = 6.0f;
	//c.vd[5] = 5.0f;
	//c.vd[6] = 4.0f;
	//c.vd[7] = 3.0f;

	//vertex<float, 8> x = traditional_mul(a, b);
	//vertex<float, 8> y = traditional_mul(b, a);

	//if (x != y)
	//	cout << "commutativity failure" << endl;

	//x = traditional_mul(traditional_mul(a, b), c);
	//y = traditional_mul(a, traditional_mul(b, c));

	//if (x != y)
	//	cout << "associativity failure" << endl;

	//x = traditional_mul(a, b + c);
	//y = traditional_mul(a, b) + traditional_mul(a, c);

	//if (x != y)
	//	cout << "distributive failure" << endl;

	//return 0;



	// Test octonion multiplication where A != B

	//vertex<float, 8> a;
	//a.vd[0] = 0.1f;
	//a.vd[1] = 0.2f;
	//a.vd[2] = 0.3f;
	//a.vd[3] = 0.4f;
	//a.vd[4] = 0.5f;
	//a.vd[5] = 0.6f;
	//a.vd[6] = 0.7f;
	//a.vd[7] = 0.8f;

	//vertex<float, 8> b;
	//b.vd[0] = 10.0f;
	//b.vd[1] = 9.0f;
	//b.vd[2] = 8.0f;
	//b.vd[3] = 7.0f;
	//b.vd[4] = 6.0f;
	//b.vd[5] = 5.0f;
	//b.vd[6] = 4.0f;
	//b.vd[7] = 3.0f;

	//vertex<float, 8> P = traditional_mul(a, b);
	//vertex<float, 8> P2 = mul(a, b);

	//for (size_t i = 0; i < 8; i++)
	//	cout << P.vd[i] << " ";

	//cout << endl;

	//for (size_t i = 0; i < 8; i++)
	//	cout << P2.vd[i] << " ";

	//cout << endl;

	//cout << P.magnitude() << " " << P2.magnitude() << endl;

	//return 0;



	// Test for 5D subalgebra	

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
	//	a.vd[0] = a0;
	//	a.vd[1] = a1;
	//	a.vd[2] = a2;
	//	a.vd[3] = a3;
	//	a.vd[4] = a4;
	//	a.vd[5] = 0;
	//	a.vd[6] = 0;
	//	a.vd[7] = 0;

	//	vertex<float, 8> b = a;

	//	vertex<float, 8> P = traditional_mul(a, b);

	//	if (P.vd[5] || P.vd[6] || P.vd[7])
	//	{
	//		cout << "Error: non-zero components!" << endl;
	//	}
	//}



	// Test sedonion numbers

	vertex<float, 16> a;
	a.vd[0] = 0.1f;
	a.vd[1] = 0.2f;
	a.vd[2] = 0.3f;
	a.vd[3] = 0.4f;
	a.vd[4] = 0.5f;
	a.vd[5] = 0.6f;
	a.vd[6] = 0.7f;
	a.vd[7] = 0.8f;
	a.vd[8] = 0.9f;
	a.vd[9] = 1.0f;
	a.vd[10] = 1.1f;
	a.vd[11] = 1.2f;
	a.vd[12] = 1.3f;
	a.vd[13] = 1.4f;
	a.vd[14] = 1.5f;
	a.vd[15] = 1.6f;

	vertex<float, 16> b;
	b.vd[0] = 10.0f;
	b.vd[1] = 9.0f;
	b.vd[2] = 8.0f;
	b.vd[3] = 7.0f;
	b.vd[4] = 6.0f;
	b.vd[5] = 5.0f;
	b.vd[6] = 4.0f;
	b.vd[7] = 3.0f;
	b.vd[8] = 2.0f;
	b.vd[9] = 1.0f;
	b.vd[10] = 0.0f;
	b.vd[11] = -1.0f;
	b.vd[12] = -2.0f;
	b.vd[13] = -3.0f;
	b.vd[14] = -4.0f;
	b.vd[15] = -5.0f;
	
	vertex<float, 16> P = traditional_mul(a, b);
	vertex<float, 16> P2 = mul(a, b);

	for (size_t i = 0; i < 16; i++)
		cout << P.vd[i] << " ";

	cout << endl;

	for (size_t i = 0; i < 16; i++)
		cout << P2.vd[i] << " ";

	cout << endl;

	cout << P.magnitude() << " " << P2.magnitude() << endl;

	return 0;

}