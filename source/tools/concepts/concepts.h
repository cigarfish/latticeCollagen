/*
 * concepts.h
 *
 * Andreas Buttenschoen
 * 2014
 *
 */

#ifndef CONCEPTS_H
#define CONCEPTS_H

#define CONCEPTS_H_VERSION 0.1

#include <type_traits>

template<bool B, typename T = void>
using Enable_if = typename std::enable_if<B, T>::type;

constexpr bool All() { return true; }

template<typename... Args>
constexpr bool All(bool b, Args... args)
{
	return b && All(args...);
}

constexpr bool Some() { return false; }

template<typename... Args>
constexpr bool Some(bool b, Args... args)
{
	return b || Some(args...);
}

template<typename T, typename U>
constexpr bool Convertible()
{
	return std::is_convertible<T, U>::value;
}

template<typename T, typename U>
constexpr bool Same()
{
	return std::is_same<T, U>::value;
}

template<typename T, typename U>
constexpr bool Different()
{
	return !Same<T, U>();
}

template<typename T>
constexpr bool Scalar()
{
	return std::is_scalar<T>::value;
}

template<typename T>
constexpr bool Pointer()
{
	return std::is_pointer<T>::value;
}

template<typename T>
constexpr bool Floating_Point()
{
	return std::is_floating_point<T>::value;
}

template<class... T>
using Common_type = typename std::common_type<T...>::type;

// p 800 28.4.4
struct substitution_failure {}; // represent a failure to declare sth

template<typename T>
struct substitution_succeeded : std::true_type
{};

template<>
struct substitution_succeeded<substitution_failure> : std::false_type
{};

template<typename T>
struct get_value_type {
private:
	template<typename X>
		static typename X::value_type check(const X&); 	 // has value_type
	static substitution_failure check(...);			 // doesn't have value_type
public:
	using type = decltype(check(std::declval<T>()));
};

template<typename T>
struct has_value_type : substitution_succeeded<typename get_value_type<T>::type>
{};

template<typename T>
constexpr bool Has_value_type()
{
	return has_value_type<T>::value;
}

// TODO ATTENTION
// can we make a template for the function?
//template<typename T>
//struct get_extent_result {
//private:
//	template<typename X>
//		static auto check(X const& x) -> decltype(extent(x));
//	static substitution_failure check(...);
//public:
//	using type = decltype(check(std::declval<T>()));
//};
//
//template<typename T>
//struct has_extent : substitution_succeeded<typename get_extent_result<T>::type>
//{};
//
//template<typename T>
//constexpr bool Has_extent()
//{
//	return has_extent<T>::value;
//}

template<typename E>
constexpr auto to_integral(E e) -> typename std::underlying_type<E>::type
{
	return static_cast<typename std::underlying_type<E>::type>(e);
}

#endif
