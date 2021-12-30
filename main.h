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

