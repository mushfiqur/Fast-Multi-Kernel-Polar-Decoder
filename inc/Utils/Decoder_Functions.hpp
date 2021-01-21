#ifndef DECODER_FUNCTIONS_HPP
#define DECODER_FUNCTIONS_HPP

#include <algorithm>

template <typename B>
constexpr B bit_init()
{
	return (B)(((B)1) << (sizeof(B) * 8 - 1));
}

template <typename R>
constexpr R init_LLR()
{
	return (R)0;
}


template <typename B, typename R>
B sgn(R val)
{
	return (B)((R(0) < val) - (val < R(0)));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <typename R>
inline R F(const R& lambda_a, const R& lambda_b)
{
	auto sign_lambda_a_b = sgn<int, R>(lambda_a * lambda_b);
	auto abs_lambda_a = (lambda_a >= 0) ? lambda_a : -lambda_a;
	auto abs_lambda_b = (lambda_b >= 0) ? lambda_b : -lambda_b;

	return (R)sign_lambda_a_b * std::min(abs_lambda_a, abs_lambda_b);
}

template <typename R>
inline R F_Ternary(const R& lambda_a, const R& lambda_b, const R& lambda_c)
{
	auto sign_lambda_a_b_c = sgn<int, R>(lambda_a * lambda_b * lambda_c);
	auto abs_lambda_a = (lambda_a >= 0) ? lambda_a : -lambda_a;
	auto abs_lambda_b = (lambda_b >= 0) ? lambda_b : -lambda_b;
	auto abs_lambda_c = (lambda_c >= 0) ? lambda_c : -lambda_c;

	if ((abs_lambda_a <= abs_lambda_b) && (abs_lambda_a <= abs_lambda_c))
		return (R)sign_lambda_a_b_c * abs_lambda_a;

	else if((abs_lambda_b <= abs_lambda_a) && (abs_lambda_b <= abs_lambda_c))
		return (R)sign_lambda_a_b_c * abs_lambda_b;

	else if ((abs_lambda_c <= abs_lambda_a) && (abs_lambda_c <= abs_lambda_b))
		return (R)sign_lambda_a_b_c * abs_lambda_c;
}

template <typename B, typename R>
inline R G(const R& lambda_a, const R& lambda_b, const B& u)
{
	return ((u == 0) ? lambda_a : -lambda_a) + lambda_b;
}

template <typename B, typename R>
inline R G_1(const R& lambda_a, const R& lambda_b, const R& lambda_c, const B& u)
{
	return ((u == 0) ? lambda_a : -lambda_a) + F(lambda_b, lambda_c);
}

template <typename B, typename R>
inline R G_2(const R& lambda_a, const R& lambda_b, const R& lambda_c, const B& u0, const B& u1)
{
	return ((u0 == 0) ? (lambda_b) : (-lambda_b)) + (((u0 ^ u1) == 0) ? (lambda_c) : (-lambda_c));
}
template <typename B, typename R>
inline B H(const R& lambda_a)
{
	return (lambda_a > (R)0) ? (B)0 : (B)1;
	// return ((lambda_a < init_LLR<R>()) * bit_init<B>());
}


#endif