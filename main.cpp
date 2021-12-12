#include <cmath>
#include <vector>
#include <iostream>

using namespace std;

class quintonion
{
public:

	quintonion(void)
	{
		vertex_data.resize(vertex_length, 0);
	}

	bool operator==(const quintonion& rhs)
	{
		bool all_equal = true;

		for (size_t i = 0; i < vertex_length; i++)
		{
			float f = fabsf(vertex_data[i] - rhs.vertex_data[i]);

			if (f > 0.0001)
			{
				all_equal = false;
				break;
			}
		}

		return all_equal;
	}

	bool operator!=(const quintonion& rhs)
	{
		return !(*this == rhs);
	}

	float magnitude(void)
	{
		float all_self_dot = 0;

		for (size_t i = 0; i < vertex_length; i++)
			all_self_dot += (vertex_data[i] * vertex_data[i]);

		return sqrtf(all_self_dot);
	}

	quintonion operator+(const quintonion& right) const
	{
		quintonion out;

		for (size_t i = 0; i < right.vertex_length; i++)
			out.vertex_data[i] = vertex_data[i] + right.vertex_data[i];

		return out;
	}

	quintonion operator/(const float& right) const
	{
		quintonion out;

		for (size_t i = 0; i < vertex_length; i++)
			out.vertex_data[i] = vertex_data[i] / right;

		return out;
	}

	size_t vertex_length = 5;
	vector<float> vertex_data;
};

quintonion conj_number_type(quintonion& in)
{
	quintonion out;

	out.vertex_data[0] = in.vertex_data[0];

	for (size_t i = 1; i < in.vertex_length; i++)
		out.vertex_data[i] = -in.vertex_data[i];

	return out;
}

quintonion pow_number_type(quintonion& in, float exponent)
{
	const float beta = exponent;

	float all_self_dot = 0;
	float imag_self_dot = 0;
	quintonion out;

	for (size_t i = 0; i < in.vertex_length; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < in.vertex_length; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < out.vertex_length; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float all_len = sqrtf(all_self_dot);
	const float imag_len = sqrtf(imag_self_dot);
	const float self_dot_beta = powf(all_self_dot, beta / 2.0f);

	out.vertex_data[0] = self_dot_beta * std::cos(beta * std::acos(in.vertex_data[0] / all_len));

	if (imag_len != 0)
	{
		for (size_t i = 1; i < out.vertex_length; i++)
			out.vertex_data[i] = in.vertex_data[i] * self_dot_beta * sin(beta * acos(in.vertex_data[0] / all_len)) / imag_len;
	}

	return out;
}


quintonion exp(const quintonion& in)
{
	float all_self_dot = 0;
	float imag_self_dot = 0;
	quintonion out;

	for (size_t i = 0; i < in.vertex_length; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < in.vertex_length; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < out.vertex_length; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float l_d = sqrtf(all_self_dot);
	const float l_e = sqrtf(imag_self_dot);

	out.vertex_data[0] = std::exp(in.vertex_data[0]) * cos(l_e);

	if (l_e != 0)
	{
		for (size_t i = 1; i < out.vertex_length; i++)
			out.vertex_data[i] = in.vertex_data[i] / l_e * std::exp(in.vertex_data[0]) * std::sin(l_e);
	}

	return out;
}

quintonion ln(const quintonion& in)
{
	float all_self_dot = 0;
	float imag_self_dot = 0;
	quintonion out;

	for (size_t i = 0; i < in.vertex_length; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < in.vertex_length; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < out.vertex_length; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float l_d = sqrtf(all_self_dot);
	const float l_e = sqrtf(imag_self_dot);

	if (in.vertex_data[0] != 0)
	{
		out.vertex_data[0] = log(l_d);
	}

	if (l_e != 0)
	{
		for (size_t i = 1; i < out.vertex_length; i++)
			out.vertex_data[i] = in.vertex_data[i] / l_e * acos(in.vertex_data[0] / l_d);
	}

	return out;
}

quintonion mul(const quintonion& in_a, const quintonion& in_b)
{
	// A*B == exp(ln(A) + ln(B))
	return exp(ln(in_a) + ln(in_b));
}

class octonion
{
public:
	inline octonion(void)
	{
		vertex_data.resize(vertex_length, 0);
	}

	bool operator==(const octonion& rhs)
	{
		bool all_equal = true;

		for (size_t i = 0; i < vertex_length; i++)
		{
			float f = fabsf(vertex_data[i] - rhs.vertex_data[i]);

			if (f > 0.0001)
			{
				all_equal = false;
				break;
			}
		}

		return all_equal;
	}

	bool operator!=(const octonion& rhs)
	{
		return !(*this == rhs);
	}

	float magnitude(void)
	{
		float all_self_dot = 0;

		for (size_t i = 0; i < vertex_length; i++)
			all_self_dot += (vertex_data[i] * vertex_data[i]);

		return sqrtf(all_self_dot);
	}

	octonion operator+(const octonion& rhs)
	{
		octonion result;

		for (size_t i = 0; i < vertex_length; i++)
			result.vertex_data[i] = vertex_data[i] + rhs.vertex_data[i];

		return result;
	}

	inline octonion(
		float src_r,
		float src_i,
		float src_j,
		float src_k,
		float src_u1,
		float src_i1,
		float src_j1,
		float src_k1)
	{
		vertex_data.resize(vertex_length, 0);

		vertex_data[0] = src_r;
		vertex_data[1] = src_i;
		vertex_data[2] = src_j;
		vertex_data[3] = src_k;
		vertex_data[4] = src_u1;
		vertex_data[5] = src_i1;
		vertex_data[6] = src_j1;
		vertex_data[7] = src_k1;
	}

	size_t vertex_length = 8;
	vector<float> vertex_data;
};

octonion traditional_mul(const octonion& qA, const octonion& qB)
{
	octonion out;

	out.vertex_data[0] = qA.vertex_data[0] * qB.vertex_data[0] - qA.vertex_data[1] * qB.vertex_data[1] - qA.vertex_data[2] * qB.vertex_data[2] - qA.vertex_data[3] * qB.vertex_data[3] - qA.vertex_data[4] * qB.vertex_data[4] - qA.vertex_data[5] * qB.vertex_data[5] - qA.vertex_data[6] * qB.vertex_data[6] - qA.vertex_data[7] * qB.vertex_data[7];
	out.vertex_data[1] = qA.vertex_data[0] * qB.vertex_data[1] + qA.vertex_data[1] * qB.vertex_data[0] + qA.vertex_data[2] * qB.vertex_data[3] - qA.vertex_data[3] * qB.vertex_data[2] + qA.vertex_data[4] * qB.vertex_data[5] - qA.vertex_data[5] * qB.vertex_data[4] - qA.vertex_data[6] * qB.vertex_data[7] + qA.vertex_data[7] * qB.vertex_data[6];
	out.vertex_data[2] = qA.vertex_data[0] * qB.vertex_data[2] - qA.vertex_data[1] * qB.vertex_data[3] + qA.vertex_data[2] * qB.vertex_data[0] + qA.vertex_data[3] * qB.vertex_data[1] + qA.vertex_data[4] * qB.vertex_data[6] + qA.vertex_data[5] * qB.vertex_data[7] - qA.vertex_data[6] * qB.vertex_data[4] - qA.vertex_data[7] * qB.vertex_data[5];
	out.vertex_data[3] = qA.vertex_data[0] * qB.vertex_data[3] + qA.vertex_data[1] * qB.vertex_data[2] - qA.vertex_data[2] * qB.vertex_data[1] + qA.vertex_data[3] * qB.vertex_data[0] + qA.vertex_data[4] * qB.vertex_data[7] - qA.vertex_data[5] * qB.vertex_data[6] + qA.vertex_data[6] * qB.vertex_data[5] - qA.vertex_data[7] * qB.vertex_data[4];
	out.vertex_data[4] = qA.vertex_data[0] * qB.vertex_data[4] - qA.vertex_data[1] * qB.vertex_data[5] - qA.vertex_data[2] * qB.vertex_data[6] - qA.vertex_data[3] * qB.vertex_data[7] + qA.vertex_data[4] * qB.vertex_data[0] + qA.vertex_data[5] * qB.vertex_data[1] + qA.vertex_data[6] * qB.vertex_data[2] + qA.vertex_data[7] * qB.vertex_data[3];
	out.vertex_data[5] = qA.vertex_data[0] * qB.vertex_data[5] + qA.vertex_data[1] * qB.vertex_data[4] - qA.vertex_data[2] * qB.vertex_data[7] + qA.vertex_data[3] * qB.vertex_data[6] - qA.vertex_data[4] * qB.vertex_data[1] + qA.vertex_data[5] * qB.vertex_data[0] - qA.vertex_data[6] * qB.vertex_data[3] + qA.vertex_data[7] * qB.vertex_data[2];
	out.vertex_data[6] = qA.vertex_data[0] * qB.vertex_data[6] + qA.vertex_data[1] * qB.vertex_data[7] + qA.vertex_data[2] * qB.vertex_data[4] - qA.vertex_data[3] * qB.vertex_data[5] - qA.vertex_data[4] * qB.vertex_data[2] + qA.vertex_data[5] * qB.vertex_data[3] + qA.vertex_data[6] * qB.vertex_data[0] - qA.vertex_data[7] * qB.vertex_data[1];
	out.vertex_data[7] = qA.vertex_data[0] * qB.vertex_data[7] - qA.vertex_data[1] * qB.vertex_data[6] + qA.vertex_data[2] * qB.vertex_data[5] + qA.vertex_data[3] * qB.vertex_data[4] - qA.vertex_data[4] * qB.vertex_data[3] - qA.vertex_data[5] * qB.vertex_data[2] + qA.vertex_data[6] * qB.vertex_data[1] + qA.vertex_data[7] * qB.vertex_data[0];

	return out;
}


octonion exp(const octonion& in)
{
	float all_self_dot = 0;
	float imag_self_dot = 0;
	octonion out;

	for (size_t i = 0; i < in.vertex_length; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < in.vertex_length; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < out.vertex_length; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float l_d = sqrtf(all_self_dot);
	const float l_e = sqrtf(imag_self_dot);

	out.vertex_data[0] = std::exp(in.vertex_data[0]) * cos(l_e);

	if (l_e != 0)
	{
		for (size_t i = 1; i < out.vertex_length; i++)
			out.vertex_data[i] = in.vertex_data[i] / l_e * std::exp(in.vertex_data[0]) * std::sin(l_e);
	}

	return out;
}

octonion ln(const octonion& in)
{
	float all_self_dot = 0;
	float imag_self_dot = 0;
	octonion out;

	for (size_t i = 0; i < in.vertex_length; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < in.vertex_length; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < out.vertex_length; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float l_d = sqrtf(all_self_dot);
	const float l_e = sqrtf(imag_self_dot);

	if (in.vertex_data[0] != 0)
	{
		out.vertex_data[0] = log(l_d);
	}

	if (l_e != 0)
	{
		for (size_t i = 1; i < out.vertex_length; i++)
			out.vertex_data[i] = in.vertex_data[i] / l_e * acos(in.vertex_data[0] / l_d);
	}

	return out;
}

octonion mul(const octonion& in_a, const octonion& in_b)
{
	// A*B == exp(ln(A) + ln(B))
	return exp(ln(in_a) + ln(in_b));
}

int main(void)
{
	// Compare pow to mul

	//quintonion a;

	//a.vertex_data[0] = 0.1f;
	//a.vertex_data[1] = 0.2f;
	//a.vertex_data[2] = 0.3f;
	//a.vertex_data[3] = 0.4f;
	//a.vertex_data[4] = 0.5f;

	//quintonion x = pow_number_type(a, 2.0f);

	//quintonion y = mul(a, a);

	//cout << x.vertex_data[0] << " " << x.vertex_data[1] << " " << x.vertex_data[2] << " " << x.vertex_data[3] << " " << x.vertex_data[4] << endl;
	//cout << y.vertex_data[0] << " " << y.vertex_data[1] << " " << y.vertex_data[2] << " " << y.vertex_data[3] << " " << y.vertex_data[4] << endl;

	//return 0;







	// Test for various attributes

	//quintonion a;
	//a.vertex_data[0] = 0.1f;
	//a.vertex_data[1] = 0.2f;
	//a.vertex_data[2] = 0.3f;
	//a.vertex_data[3] = 0.4f;
	//a.vertex_data[4] = 0.5f;

	//quintonion b;
	//b.vertex_data[0] = 1.0f;
	//b.vertex_data[1] = 0.9f;
	//b.vertex_data[2] = 0.8f;
	//b.vertex_data[3] = 0.7f;
	//b.vertex_data[4] = 0.6f;

	//quintonion c;
	//b.vertex_data[0] = 10.0f;
	//b.vertex_data[1] = 9.0f;
	//b.vertex_data[2] = 8.0f;
	//b.vertex_data[3] = 7.0f;
	//b.vertex_data[4] = 6.0f;

	//quintonion x = mul(a, b);
	//quintonion y = mul(b, a);

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






	// Test octonion multiplication where A != B
	octonion A(0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f);
	octonion B(10.0f, 9.0f, 8.0f, 7.0f, 6.0f, 0.5f, 0.4f, 0.3f);

	octonion P = traditional_mul(A, B);
	octonion P2 = mul(A, B);

	for (size_t i = 0; i < P.vertex_length; i++)
		cout << P.vertex_data[i] << " ";

	cout << endl;

	for (size_t i = 0; i < P2.vertex_length; i++)
		cout << P2.vertex_data[i] << " ";

	cout << endl;

	cout << P.magnitude() << " " << P2.magnitude() << endl;

	return 0;
	





	// Test 3) Test for subalgebra	
	//
	//srand(time(0));

	//for (size_t num_tries = 0; num_tries < 10000; num_tries++)
	//{
	//	float a0 = (rand() % RAND_MAX);

	//	if (rand() % 2)
	//		a0 = -a0;

	//	float a1 = (rand() % RAND_MAX);

	//	if (rand() % 2)
	//		a1 = -a1;

	//	float a2 = (rand() % RAND_MAX);

	//	if (rand() % 2)
	//		a2 = -a2;

	//	float a3 = (rand() % RAND_MAX);

	//	if (rand() % 2)
	//		a3 = -a3;

	//	float a4 = (rand() % RAND_MAX);

	//	if (rand() % 2)
	//		a4 = -a4;

	//	octonion A(a0, a1, a2, a3, a4, 0, 0, 0);
	//	octonion B = A;

	//	octonion P = traditional_mul(A, B);

	//	if (P.vertex_data[5] || P.vertex_data[6] || P.vertex_data[7])
	//	{
	//		cout << "Error: non-zero components!" << endl;
	//	}
	//}


	return 0;
}