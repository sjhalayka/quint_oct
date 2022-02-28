#include "main.h"

float rand_num(void)
{
	return static_cast<float>(rand()) / static_cast<float>(RAND_MAX); // 1;
}

long double rand_num_ld(void)
{
	return static_cast<long double>(rand()) / static_cast<long double>(RAND_MAX); // 1;
}

// Function declarations (see definitions below main())
void compare_square_function(void);
void compare_real_numbers(void);
void compare_complex_numbers(void);
void compare_quaternion_numbers(void);
void compare_quintonion_pow_to_mul(void);
void test_quintonions(void);
void test_octonion_new_multiplication(void);
void test_octonion_traditional_multiplication(void);
void test_octonion_multiplication(void);
void test_for_5D_subalgebra(void);
void test_octonion_multiplication_speed(void);
void test_sedonion_multiplication(void);
void test_sedonion_traditional_multiplication(void);
void test_sedenion_multiplication_speed(void);
void test_pathion_multiplication_speed(void);

void test_octonion_pow_speed(void);
void test_sedenion_pow_speed(void);
void test_pathion_pow_speed(void);

void test_pathion_multiplication(void);

void test_power_associativity(void);

// See:
// https://theworld.com/~sweetser/quaternions/intro/tools/tools.html
// https://ece.uwaterloo.ca/~dwharder/C++/CQOST/src/
// https://math.stackexchange.com/a/2554/745219
// http://paulbourke.net/fractals/quatjulia/
// https://math.stackexchange.com/a/1047246/745219
// https://www.mathworks.com/help/nav/ref/quaternion.log.html
// https://www.mathworks.com/help/nav/ref/quaternion.exp.html
// https://www.mathworks.com/help/robotics/ref/quaternion.power.html



int main(void)
{
	srand(static_cast<unsigned>(time(0)));

	cout << fixed << endl;
	cout << setprecision(10) << endl;
	 
//	test_power_associativity();
	test_octonion_traditional_multiplication();

	return 0;



	//compare_real_numbers();

	//compare_complex_numbers();

	//compare_quaternion_numbers();

	//test_octonion_multiplication();

	//test_sedonion_multiplication();

	// test_pathion_multiplication();

//	return 0;



	// compare_square_function();
	//compare_real_numbers();
	//compare_complex_numbers();
	
	//compare_quaternion_numbers();
	//compare_quintonion_pow_to_mul();

	//test_octonion_multiplication();

	//compare_quintonion_pow_to_mul();
	
	//test_quintonions();
	
	//test_octonion_new_multiplication();
	
	//test_octonion_traditional_multiplication();
	
	//test_octonion_multiplication();

	//test_octonion_multiplication_speed();
	
	//test_sedonion_multiplication();
	
	//test_sedonion_traditional_multiplication();
	
	//test_sedenion_multiplication_speed();
	//test_pathion_multiplication_speed();

//	return 0;
}

void compare_square_function(void)
{
	cout << "Comparing square function" << endl;

	vertex<float, 5> x;
	x.vd[0] = 1.0f;
	x.vd[1] = 2.0f;
	x.vd[2] = 3.0f;
	x.vd[3] = 4.0f;
	x.vd[4] = 5.0f;

	vertex<float, 5> y = mul(x, x); // y = pow(x, 2.0f);
	vertex<float, 5> z = square(x);

	for (size_t i = 0; i < 5; i++)
		cout << y.vd[i] << " ";

	cout << endl;

	for (size_t i = 0; i < 5; i++)
		cout << z.vd[i] << " ";

	cout << endl;
}

void compare_real_numbers(void)
{
	cout << "Comparing real numbers" << endl;

	vertex<float, 1> a;
	a.vd[0] = -rand_num() * 0.1234f;

	vertex<float, 1> b;// = a;
	b.vd[0] = rand_num() * 0.1234f;

	vertex<float, 1> x = mul(a, b);
	vertex<float, 1> y = traditional_mul(a, b);

	cout << x.vd[0] << endl;
	cout << y.vd[0] << endl;

	cout << x.magnitude() << " " << y.magnitude() << endl;


	cout << x.magnitude() / y.magnitude() << endl;

	cout << endl;
}

void compare_complex_numbers(void)
{
	cout << "Comparing complex numbers" << endl;

	vertex<float, 2> a;
	a.vd[0] = rand_num() * 0.1f;
	a.vd[1] = -rand_num() * 0.2f;

	vertex<float, 2> b = a;
	b.vd[0] = -rand_num() * 1.0f;
	b.vd[1] = rand_num() * 0.9f;

	vertex<float, 2> x = mul(a, b);
	vertex<float, 2> y = traditional_mul(a, b);

	cout << x.vd[0] << " " << x.vd[1] << endl;
	cout << y.vd[0] << " " << y.vd[1] << endl;

	cout << x.magnitude() << " " << y.magnitude() << endl;

	cout << x.magnitude() / y.magnitude() << endl;

	cout << endl;



	//complex<float> cf_a(a.vd[0], a.vd[1]);
	//complex<float> cf_b(b.vd[0], b.vd[1]);

	//complex<float> cf_x = cf_a * cf_b;

	//cout << cf_x.real() << " " << cf_x.imag() << endl;

	//cout << "Comparing complex numbers exp()" << endl;

	//a = exp(a);
	//cf_a = exp(cf_a);

	//cout << a.vd[0] << " " << a.vd[1] << endl;
	//cout << cf_a.real() << " " << cf_a.imag() << endl;

	//cout << "Comparing complex numbers log()" << endl;

	//complex<float> cf_b(b.vd[0], b.vd[1]);

	//b = log(b);
	//cf_b = log(cf_b);

	//cout << b.vd[0] << " " << b.vd[1] << endl;
	//cout << cf_b.real() << " " << cf_b.imag() << endl;

	//cout << endl;
}

void compare_quaternion_numbers(void)
{
	cout << "Comparing quaternion numbers" << endl;

	vertex<float, 4> a;
	a.vd[0] = -rand_num() * 0.1f;
	a.vd[1] = rand_num() * 0.2f;
	a.vd[2] = rand_num() * 0.3f;
	a.vd[3] = -rand_num() * 0.4f;

	vertex<float, 4> b;// = a;
	b.vd[0] = -rand_num() * 1.0f;
	b.vd[1] = rand_num() * 0.9f;
	b.vd[2] = -rand_num() * 0.8f;
	b.vd[3] = rand_num() * 0.7f;

	vertex<float, 4> x = mul(a, b);
	vertex<float, 4> y = traditional_mul(a, b);

	cout << x.vd[0] << " " << x.vd[1] << " " << x.vd[2] << " " << x.vd[3] << endl;
	cout << y.vd[0] << " " << y.vd[1] << " " << y.vd[2] << " " << y.vd[3] << endl;

	cout << "Magnitudes:" << endl;
	cout << x.magnitude() << " " << y.magnitude() << endl;


	cout << x.magnitude() / y.magnitude() << endl;

	cout << endl;
}

void compare_quintonion_pow_to_mul(void)
{
	cout << "Comparing quintonion pow to new multiplication" << endl;

	vertex<float, 5> a;

	a.vd[0] = 0.1f;
	a.vd[1] = 0.2f;
	a.vd[2] = 0.3f;
	a.vd[3] = 0.4f;
	a.vd[4] = 0.5f;

	vertex<float, 5> x = pow(a, 2.0f);
	vertex<float, 5> y = mul(a, a);

	cout << x.vd[0] << " " << x.vd[1] << " " << x.vd[2] << " " << x.vd[3] << " " << x.vd[4] << endl;
	cout << y.vd[0] << " " << y.vd[1] << " " << y.vd[2] << " " << y.vd[3] << " " << y.vd[4] << endl;

	cout << endl;
}

void test_quintonions(void)
{
	cout << "Test quintonion attributes:" << endl;

	vertex<float, 5> a;
	a.vd[0] = 0.1f;
	a.vd[1] = 0.2f;
	a.vd[2] = 0.3f;
	a.vd[3] = 0.4f;
	a.vd[4] = 0.5f;

	vertex<float, 5> b;
	b.vd[0] = 1.0f;
	b.vd[1] = 0.9f;
	b.vd[2] = 0.8f;
	b.vd[3] = 0.7f;
	b.vd[4] = 0.6f;

	vertex<float, 5> c;
	c.vd[0] = 10.0f;
	c.vd[1] = 9.0f;
	c.vd[2] = 8.0f;
	c.vd[3] = 7.0f;
	c.vd[4] = 6.0f;

	vertex<float, 5> x = mul(a, b);
	vertex<float, 5> y = mul(b, a);

	if (x != y)
		cout << "commutativity failure" << endl;
	else
		cout << "commutativity OK" << endl;

	x = mul(mul(a, b), c);
	y = mul(a, mul(b, c));

	if (x != y)
		cout << "associativity failure" << endl;
	else
		cout << "associativity OK" << endl;

	x = mul(a, b + c);
	y = mul(a, b) + mul(a, c);

	if (x != y)
		cout << "distributativity failure" << endl;
	else
		cout << "distributativity OK" << endl;

	cout << endl;
}

void test_octonion_new_multiplication(void)
{
	cout << "Test octonion new multiplication attributes:" << endl;

	vertex<float, 8> a;
	a.vd[0] = 0.1f;
	a.vd[1] = 0.2f;
	a.vd[2] = 0.3f;
	a.vd[3] = 0.4f;
	a.vd[4] = 0.5f;
	a.vd[5] = 0.6f;
	a.vd[6] = 0.7f;
	a.vd[7] = 0.8f;

	vertex<float, 8> b;
	b.vd[0] = 1.0f;
	b.vd[1] = 0.9f;
	b.vd[2] = 0.8f;
	b.vd[3] = 0.7f;
	b.vd[4] = 0.6f;
	b.vd[5] = 0.5f;
	b.vd[6] = 0.4f;
	b.vd[7] = 0.3f;

	vertex<float, 8> c;
	c.vd[0] = 10.0f;
	c.vd[1] = 9.0f;
	c.vd[2] = 8.0f;
	c.vd[3] = 7.0f;
	c.vd[4] = 6.0f;
	c.vd[5] = 5.0f;
	c.vd[6] = 4.0f;
	c.vd[7] = 3.0f;

	vertex<float, 8> x = mul(a, b);
	vertex<float, 8> y = mul(b, a);

	if (x != y)
		cout << "commutativity failure" << endl;
	else
		cout << "commutativity OK" << endl;

	x = mul(mul(a, b), c);
	y = mul(a, mul(b, c));

	if (x != y)
		cout << "associativity failure" << endl;
	else
		cout << "associativity OK" << endl;

	x = mul(a, b + c);
	y = mul(a, b) + mul(a, c);

	if (x != y)
		cout << "distributativity failure" << endl;
	else
		cout << "distributativity OK" << endl;

	cout << endl;
}

void test_octonion_traditional_multiplication(void)
{
	cout << "Test octonion traditional multiplication attributes:" << endl;

	vertex<float, 8> a;
	a.vd[0] = 0.1f;
	a.vd[1] = 0.2f;
	a.vd[2] = 0.3f;
	a.vd[3] = 0.4f;
	a.vd[4] = 0.5f;
	a.vd[5] = 0.6f;
	a.vd[6] = 0.7f;
	a.vd[7] = 0.8f;

	vertex<float, 8> b;
	b.vd[0] = 1.0f;
	b.vd[1] = 0.9f;
	b.vd[2] = 0.8f;
	b.vd[3] = 0.7f;	
	b.vd[4] = 0.6f;
	b.vd[5] = 0.5f;
	b.vd[6] = 0.4f;
	b.vd[7] = 0.3f;

	vertex<float, 8> c;
	c.vd[0] = 10.0f;
	c.vd[1] = 9.0f;
	c.vd[2] = 8.0f;
	c.vd[3] = 7.0f;
	c.vd[4] = 6.0f;
	c.vd[5] = 5.0f;
	c.vd[6] = 4.0f;
	c.vd[7] = 3.0f;

	vertex<float, 8> x = traditional_mul(a, b);
	vertex<float, 8> y = traditional_mul(b, a);

	if (x != y)
		cout << "commutativity failure" << endl;
	else
		cout << "commutativity OK" << endl;

	x = traditional_mul(traditional_mul(a, b), c);
	y = traditional_mul(a, traditional_mul(b, c));

	if (x != y)
		cout << "associativity failure" << endl;
	else
		cout << "associativity OK" << endl;

	x = traditional_mul(a, b + c);
	y = traditional_mul(a, b) + traditional_mul(a, c);

	if (x != y)
		cout << "distributativity failure" << endl;
	else
		cout << "distributativity OK" << endl;

	cout << endl;
}

void test_octonion_multiplication(void)
{
	cout << "Test octonion multiplications:" << endl;

	vertex<float, 8> a;
	a.vd[0] = rand_num() * 0.1f;
	a.vd[1] = rand_num() * 0.2f;
	a.vd[2] = rand_num() * 0.3f;
	a.vd[3] = rand_num() * 0.4f;
	a.vd[4] = rand_num() * 0.5f;
	a.vd[5] = rand_num() * 0.6f;
	a.vd[6] = rand_num() * 0.7f;
	a.vd[7] = rand_num() * 0.8f;

	vertex<float, 8> b;// = a;
	b.vd[0] = rand_num() * 10.0f;
	b.vd[1] = rand_num() * 9.0f;
	b.vd[2] = rand_num() * 8.0f;
	b.vd[3] = rand_num() * 7.0f;
	b.vd[4] = rand_num() * 6.0f;
	b.vd[5] = rand_num() * 5.0f;
	b.vd[6] = rand_num() * 4.0f;
	b.vd[7] = rand_num() * 3.0f;

	vertex<float, 8> P = traditional_mul(a, b);
	vertex<float, 8> P2 = mul(a, b);

	for (size_t i = 0; i < 8; i++)
		cout << P.vd[i] << " ";

	cout << endl;

	for (size_t i = 0; i < 8; i++)
		cout << P2.vd[i] << " ";

	cout << endl;

	cout << "Magnitudes:" << endl;

	cout << P.magnitude() << " " << P2.magnitude() << endl;


	cout << P2.magnitude() / P.magnitude() << endl;

	cout << endl;
}



void test_octonion_multiplication_speed(void)
{
	std::chrono::high_resolution_clock::time_point start_time, end_time;
	std::chrono::duration<float, std::milli> elapsed;

	start_time = std::chrono::high_resolution_clock::now();

	const size_t num_iterations = 10000000;

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 8> a;

		for (size_t i = 0; i < 8; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 8> b;

		for (size_t i = 0; i < 8; i++)
			b.vd[i] = -1.0f * (i + 1);

		vertex<float, 8> x = traditional_mul(a, b);
	}


	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Octonion O(n^2) duration: " << elapsed.count() / 1000.0f << " seconds" << endl;


	start_time = std::chrono::high_resolution_clock::now();

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 8> a;

		for (size_t i = 0; i < 8; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 8> b;

		for (size_t i = 0; i < 8; i++)
			b.vd[i] = -1.0f * (i + 1);

		vertex<float, 8> x = mul(a, b);
	}

	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Octonion O(n) duration: " << elapsed.count() / 1000.0f << " seconds" << endl;
}

void test_sedonion_multiplication(void)
{
	cout << "Test sedonion multiplications:" << endl;

	vertex<float, 16> a;
	a.vd[0] = rand_num() * 0.1f;
	a.vd[1] = rand_num() * 0.2f;
	a.vd[2] = rand_num() * 0.3f;
	a.vd[3] = rand_num() * 0.4f;
	a.vd[4] = rand_num() * 0.5f;
	a.vd[5] = rand_num() * 0.6f;
	a.vd[6] = rand_num() * 0.7f;
	a.vd[7] = rand_num() * 0.8f;
	a.vd[8] = rand_num() * 0.9f;
	a.vd[9] = rand_num() * 1.0f;
	a.vd[10] = rand_num() * 1.1f;
	a.vd[11] = rand_num() * 1.2f;
	a.vd[12] = rand_num() * 1.3f;
	a.vd[13] = rand_num() * 1.4f;
	a.vd[14] = rand_num() * 1.5f;
	a.vd[15] = rand_num() * 1.6f;

	vertex<float, 16> b;// = a;
	b.vd[0] = rand_num() * 10.0f;
	b.vd[1] = rand_num() * 9.0f;
	b.vd[2] = rand_num() * 8.0f;
	b.vd[3] = rand_num() * 7.0f;
	b.vd[4] = rand_num() * 6.0f;
	b.vd[5] = rand_num() * 5.0f;
	b.vd[6] = rand_num() * 4.0f;
	b.vd[7] = rand_num() * 3.0f;
	b.vd[8] = rand_num() * 2.0f;
	b.vd[9] = rand_num() * 1.0f;
	b.vd[10] = rand_num() * 0.0f;
	b.vd[11] = rand_num() * -1.0f;
	b.vd[12] = rand_num() * -2.0f;
	b.vd[13] = rand_num() * -3.0f;
	b.vd[14] = rand_num() * -4.0f;
	b.vd[15] = rand_num() * -5.0f;

	vertex<float, 16> P = traditional_mul(a, b);
	vertex<float, 16> P2 = mul(a, b);

	for (size_t i = 0; i < 16; i++)
		cout << P.vd[i] << " ";

	cout << endl;

	for (size_t i = 0; i < 16; i++)
		cout << P2.vd[i] << " ";

	cout << endl;

	cout << "Magnitudes:" << endl;

	cout << P.magnitude() << " " << P2.magnitude() << endl;


	cout << P2.magnitude() / P.magnitude() << endl;

	cout << endl;
}

void test_sedonion_traditional_multiplication(void)
{
	cout << "Test sedenion traditional multiplication attributes:" << endl;

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

	vertex<float, 16> c;
	c.vd[0] = 1.0f;
	c.vd[1] = 2.0f;
	c.vd[2] = 3.0f;
	c.vd[3] = 4.0f;
	c.vd[4] = 5.0f;
	c.vd[5] = 6.0f;
	c.vd[6] = 7.0f;
	c.vd[7] = 8.0f;
	c.vd[8] = 9.0f;
	c.vd[9] = 10.0f;
	c.vd[10] = 11.0f;
	c.vd[11] = 12.0f;
	c.vd[12] = 13.0f;
	c.vd[13] = 14.0f;
	c.vd[14] = 15.0f;
	c.vd[15] = 16.0f;

	vertex<float, 16> x = traditional_mul(a, b);
	vertex<float, 16> y = traditional_mul(b, a);

	if (x != y)
		cout << "commutativity failure" << endl;
	else
		cout << "commutativity OK" << endl;

	x = traditional_mul(traditional_mul(a, b), c);
	y = traditional_mul(a, traditional_mul(b, c));

	if (x != y)
		cout << "associativity failure" << endl;
	else
		cout << "associativity OK" << endl;

	x = traditional_mul(a, b + c);
	y = traditional_mul(a, b) + traditional_mul(a, c);

	if (x != y)
		cout << "distributativity failure" << endl;
	else
		cout << "distributativity OK" << endl;

	cout << endl;
}

void test_sedenion_multiplication_speed(void)
{
	std::chrono::high_resolution_clock::time_point start_time, end_time;
	std::chrono::duration<float, std::milli> elapsed;

	start_time = std::chrono::high_resolution_clock::now();

	const size_t num_iterations = 10000000;

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 16> a;

		for (size_t i = 0; i < 16; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 16> b;

		for (size_t i = 0; i < 16; i++)
			b.vd[i] = -1.0f * (i + 1);

		vertex<float, 16> x = traditional_mul(a, b);
	}


	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Sedenion O(n^2) duration: " << elapsed.count() / 1000.0f << " seconds" << endl;


	start_time = std::chrono::high_resolution_clock::now();

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 16> a;

		for (size_t i = 0; i < 16; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 16> b;

		for (size_t i = 0; i < 16; i++)
			b.vd[i] = -1.0f * (i + 1);

		vertex<float, 16> x = mul(a, b);
	}

	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Sedenion O(n) duration: " << elapsed.count() / 1000.0f << " seconds" << endl;
}

void test_pathion_multiplication_speed(void)
{
	std::chrono::high_resolution_clock::time_point start_time, end_time;
	std::chrono::duration<float, std::milli> elapsed;

	start_time = std::chrono::high_resolution_clock::now();

	const size_t num_iterations = 10000000;

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 32> a;

		for (size_t i = 0; i < 32; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 32> b;

		for (size_t i = 0; i < 32; i++)
			b.vd[i] = -1.0f * (i + 1);

		vertex<float, 32> x = traditional_mul(a, b);
	}

	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Pathion O(n^2) duration: " << elapsed.count() / 1000.0f << " seconds" << endl;


	start_time = std::chrono::high_resolution_clock::now();

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 32> a;

		for (size_t i = 0; i < 32; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 32> b;

		for (size_t i = 0; i < 32; i++)
			b.vd[i] = -1.0f * (i + 1);

		vertex<float, 32> x = mul(a, b);
	}

	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Pathion O(n) duration: " << elapsed.count() / 1000.0f << " seconds" << endl;
}


void test_octonion_pow_speed(void)
{
	std::chrono::high_resolution_clock::time_point start_time, end_time;
	std::chrono::duration<float, std::milli> elapsed;

	start_time = std::chrono::high_resolution_clock::now();

	const size_t num_iterations = 10000000;

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 8> a;

		for (size_t i = 0; i < 8; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 8> b = a;

		vertex<float, 8> x = mul(a, b);
	}

	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Octonion multiplication duration: " << elapsed.count() / 1000.0f << " seconds" << endl;


	start_time = std::chrono::high_resolution_clock::now();

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 8> a;

		for (size_t i = 0; i < 8; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 8> b = a;

		vertex<float, 8> x = pow(a, 2.0f);
	}

	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Octonion pow duration: " << elapsed.count() / 1000.0f << " seconds" << endl;
}

void test_sedenion_pow_speed(void)
{
	std::chrono::high_resolution_clock::time_point start_time, end_time;
	std::chrono::duration<float, std::milli> elapsed;

	start_time = std::chrono::high_resolution_clock::now();

	const size_t num_iterations = 10000000;

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 16> a;

		for (size_t i = 0; i < 16; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 16> b = a;

		vertex<float, 16> x = mul(a, b);
	}

	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Sedenion multiplication duration: " << elapsed.count() / 1000.0f << " seconds" << endl;


	start_time = std::chrono::high_resolution_clock::now();

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 16> a;

		for (size_t i = 0; i < 16; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 16> b = a;

		vertex<float, 16> x = pow(a, 2.0f);
	}

	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Sedenion pow duration: " << elapsed.count() / 1000.0f << " seconds" << endl;
}

void test_pathion_pow_speed(void)
{
	std::chrono::high_resolution_clock::time_point start_time, end_time;
	std::chrono::duration<float, std::milli> elapsed;

	start_time = std::chrono::high_resolution_clock::now();

	const size_t num_iterations = 10000000;

	size_t exponent = 25;

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 32> a;

		for (size_t i = 0; i < 32; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 32> b = a;

		for(size_t i = 0; i < (exponent - 1); i++)
			a = mul(a, b);
	}



	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Pathion multiplication duration: " << elapsed.count() / 1000.0f << " seconds" << endl;


	start_time = std::chrono::high_resolution_clock::now();

	for (size_t i = 0; i < num_iterations; i++)
	{
		vertex<float, 32> a;

		for (size_t i = 0; i < 32; i++)
			a.vd[i] = 0.1f * (i + 1);

		vertex<float, 32> b = a;

		vertex<float, 32> x = pow(a, static_cast<float>( exponent));
	}

	end_time = std::chrono::high_resolution_clock::now();

	elapsed = end_time - start_time;

	cout << "Pathion pow duration: " << elapsed.count() / 1000.0f << " seconds" << endl;
}


void test_pathion_multiplication(void)
{
	cout << "Test pathion multiplications:" << endl;

	vertex<float, 32> a;

	for (size_t i = 0; i < 32; i++)
		a.vd[i] = rand_num()/* / static_cast<float>(RAND_MAX)*/ * 0.1f * (i + 1);

	vertex<float, 32> b;// = a;

	for (size_t i = 0; i < 32; i++)
		b.vd[i] = rand_num()/* / static_cast<float>(RAND_MAX)*/ * 0.1f * (i + 1);

	vertex<float, 32> P = traditional_mul(a, b);
	vertex<float, 32> P2 = mul(a, b);

	for (size_t i = 0; i < 32; i++)
		cout << P.vd[i] << " ";

	cout << endl;

	for (size_t i = 0; i < 32; i++)
		cout << P2.vd[i] << " ";

	cout << endl;

	cout << "Magnitudes:" << endl;

	cout << P.magnitude() << " " << P2.magnitude() << endl;

	cout << P2.magnitude() / P.magnitude() << endl;

	cout << endl;
}

void test_power_associativity(void)
{
	// Test power associativity for traditional multiplication
	const size_t n = 32; // choose any n from 2, 4, 8, 16, 32

	vertex<long double, n> a_base;

	for (size_t i = 0; i < n; i++)
		a_base.vd[i] = rand_num_ld() * 10 * static_cast<long double>(i + 1);

	vertex<long double, n> a = traditional_mul(a_base, a_base);

	cout << get_traditional_commutator(a, a_base).magnitude() << endl;
}