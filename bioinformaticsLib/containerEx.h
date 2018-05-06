// ***************************************************************************
// ContainerEx.h (c) 2018 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 13 January 2018
// ---------------------------------------------------------------------------
// Various additions to the standard template library containers.
//
//  #defines that add names to the member variables of the classes used for iteration, 
//    to get rid of the .first and.second
//
// Additions to <map> and <vector> classes to add extra useful functions
// ***************************************************************************

#ifndef CONTAINEREX_HEADER
#define CONTAINEREX_HEADER

#include <string>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <numeric>

/*********************************************
	MAKE_NAMED_PAIR: An extension for pair where a class _name is created where instead of using .first and .second to access the two parts
	.<_first> and .<second> are used

	Need to override constructors so that the reference variables are initialised in each case.  Includes a && move constructor
*/


#define MAKE_NAMED_PAIR(_name,type_f,_first,type_s,_second) \
struct _name  : std::pair<type_f,type_s> { \
	type_f & _first; type_s & _second ; \
	_name(const type_f & _a,const type_s & _b) :  std::pair<type_f,type_s>(_a,_b), _first(first) ,_second(second) {}; \
	_name(const _name & nit) :  std::pair<type_f,type_s>(nit), _first(first) ,_second(second) {}; \
	_name(_name && nit) :  std::pair<type_f,type_s>(std::move(nit)), _first(first) ,_second(second) {}; \
	_name() :  _first(first) ,_second(second) {}; \
		};

/*********************************************
	Some macros for use when using C++ 2011 range based iteration of maps.  The default is the value returned is a pair, 
	  which is then accessed using .first and .second
	ADD_MAPPAIR creates a 'Pair' class such that the key is returned using _first() and a reference to the value with _second();
	ADD_MAPPAIR2 creates a 'Pair' class such that the key is returned using _first and a reference to the value with _second;
		
	ADD_CONSTMAPPAIR is for use when the map is const = unmodifiable.
*/


#define ADD_MAPPAIR(_first,_second) \
struct Pair  { \
	std::pair<const key_type,mapped_type> & _p; \
	Pair(std::pair<const key_type,mapped_type> & _a) : _p(_a){}; \
	key_type _first () {return _p.first;}; \
	mapped_type & _second ()  {return _p.second;};  \
}; 

#define ADD_MAPPAIR_VAR(_first,_second) \
struct Pair { \
	const key_type & _first; mapped_type & _second ; \
	Pair(std::pair<const key_type,mapped_type> & _a) :  _first(_a.first) ,_second(_a.second) {}; \
};

#define ADD_CONSTMAPPAIR_VAR(_first,_second) \
struct Pair { \
	const key_type & _first; const mapped_type & _second ; \
	Pair(const std::pair<const key_type,mapped_type> & _a) :  _first(_a.first) ,_second(_a.second) {}; \
};


#define ADD_CONSTMAPPAIR(_first,_second) \
struct Pair  { \
	const std::pair<const key_type,mapped_type> & _p; \
	Pair(const std::pair<const key_type,mapped_type> & _a) : _p(_a){}; \
	key_type _first () const {return _p.first;}; \
	const mapped_type &_second () const {return _p.second;};  \
}; 

/*********************************************
	For use when iterating maps.   Defines a new type 'Iterator' and 'ReverseIterator' to be used 
	in place of 'iterator' and reverse_iterator, so that the
	key and the value can be found using <_first>() and <_second>();

*/


#define ADD_ITER(_first,_second) \
struct Iterator : public iterator { \
	Iterator() {}; \
	Iterator(const iterator & _a) : iterator(_a){}; \
	key_type _first() { return (*this)->first;}; \
	mapped_type & _second() { return (*this)->second;}; \
	Iterator operator ++ (int) { Iterator x = *this; iterator::operator ++(1); return x;};\
	Iterator & operator ++ () { iterator::operator ++(); return *this;};\
}; \
Iterator Begin(){ return begin();}; \
Iterator End() { return end(); }; \
Iterator Find(const key_type& key) { return find(key); };

#define ADD_CONST_ITER(_first,_second) \
struct Iterator : public const_iterator { \
	Iterator() {}; \
	Iterator(const const_iterator & _a) : const_iterator(_a){}; \
	key_type _first() { return (*this)->first;}; \
	const mapped_type & _second() { return (*this)->second;}; \
	Iterator operator ++ (int) { Iterator x = *this; const_iterator::operator ++(1); return x;};\
	Iterator & operator ++ () { const_iterator::operator ++(); return *this;};\
};  \
Iterator Begin() { return begin(); }; \
Iterator End() { return end(); };

#define ADD_REVITER(_first,_second) \
struct ReverseIterator : public reverse_iterator { \
	ReverseIterator() {}; \
	ReverseIterator(const reverse_iterator & _a) : reverse_iterator(_a){}; \
	key_type _first() { return (*this)->first;}; \
	mapped_type & _second() { return (*this)->second;}; \
}; 

/*********************************************
	An extension to the map class so that the [] operator initialises the value to zero if the 
	entry did not previously exist

*/


template<class _Kty,typename _Ty> 
class  mapInitDef : public std::map<_Kty,_Ty>
{
	_Ty initValue;
public:
	mapInitDef(_Ty initValue = 0): initValue(initValue){};
	_Ty & operator [] (_Kty n)
	{
		//	Need explicit this -> when inheriting from a template class because otherwise the compiler
		//	may not be able to unambiguously determine the function to use.
		//	MSVCC does not need the this ->
		auto i = this -> find(n);
		if (i != this -> end())
			return i->second;
		return (this -> emplace(n,initValue).first->second);
	}
};

template<class _Kty,typename _Ty> 
class mapZeroDef : public mapInitDef<_Kty,_Ty>
{
public:
	mapZeroDef() : mapInitDef<_Kty,_Ty>(0){};
};

template<typename _Ty> class  vectorZeroDef : public std::vector<_Ty>
{
public:
	vectorZeroDef(size_t s) : std::vector<_Ty>(s)
	{
		for (auto & i: (*this))
			i = 0;
	}
};

//	An extension of the set class that adds a constructor and an add function that allows multiple entries to be added at a time
//	In the case of sets of strings, strings, character arrays and const char * can all be used as they are all converted
//	to strings
template <typename _Ty,class _Alloc = std::allocator<_Ty> >  class setEx : public std::set<_Ty>
{
	typedef std::set<_Ty> _set;
//	void add()	{};		//Needed for when the add(...) method runs out of arguments
public:
	setEx(){};
	//	Allows the initialisation of a set using an aggregate, e.g. setEx<string> a{{"hello","world"}};  
	//	Also supports = operator, e.g. a = {{"good morning","world"}};
	setEx(const _set & _V) :_set(_V) {};

	setEx & operator += (setEx<_Ty> a) {
		_set::insert(a.begin(), a.end());
		return *this;
	};


/*
//	Add one or more things
	template<typename... Rest> void add(const _Ty & value,const Rest & ... rest)	{
		this->insert(value);
		add(rest...);
	};

//	Where some of the things can be setExs
	template<typename... Rest> void add(const setEx & value,const Rest & ... rest)	{
		this->insert(value.begin(),value.end());
		add(rest...);
	};
*/

/*	template<typename... Values> setEx(const _Ty & value,const Values & ... values) {
		this->insert(value);
		add(values...);
	};
*/
	bool contains(const _Ty & value) const
	{
		return (_set::find(value) != _set::end());
	}

};

template <typename _Ty, class _Alloc = std::allocator<_Ty> >  class multisetEx : public std::multiset<_Ty>
{
};

//	Use a full definition so that the template parser class function that handles vectors will also handle vectorExs


template <typename _Ty,class _Alloc = std::allocator<_Ty> >  class vectorEx : public std::vector<_Ty>
{
public:
	typedef std::vector<_Ty> _vec;

	vectorEx() {};
	vectorEx(const size_t size) : _vec(size) {};
	vectorEx(const size_t size,const _Ty val) : _vec(size,val) {};
	//	This allows a vectorEx to be constructed from a vector, and so also allows it
	//	to be constructed from an aggregate, e.g. vectorEx<int> a{0,1,3}; or a = {{1,2,5}};
	vectorEx(const _vec & a) : _vec(a) {};
	vectorEx(const std::initializer_list<_Ty> a) : _vec(a) {};


//	Adds a vectorEx to the end of another vectorEx.  This can be done using an initialser list
//	vectorEx<int> a({1,2,3}),b({4,5,6});
//	a += b;
//	a += {7,8,9};
	vectorEx<_Ty> & operator += (const _vec & a) {
		_vec::insert(_vec::end(), a.begin(), a.end());
		return *this;
	};

	vectorEx<_Ty> & operator = (const vectorEx<_Ty> & _A)
	{
		_vec::operator = (_A);
		return *this;
	}
	vectorEx<_Ty> operator / (const vectorEx<_Ty> & _A) const
	{
		size_t S = std::min(_vec::size(),_A.size());
		vectorEx<_Ty> _V(S);
		for (size_t i = 0;i < S;i++)
			_V[i] = _vec::at(i)/_A[i];
		return _V;
	}

	_Ty mean() const
	{
		return accumulate(_vec::begin(), _vec::end(), 0.0) / _vec::size();
	}
	size_t find (const _Ty & value) const 
	{
		for (const auto & i : *this)
			if (i == value) return i;
		return std::string::npos;
	}
	size_t contains (const _Ty & value) const 
	{
		for (const auto & i : *this)
			if (i == value)	return true;
		return false;
	}
	_Ty max() const
	{
		return *max_element(_vec::begin(), _vec::end());
	}


//	operator vector<_Ty> & () {return This;};

/*	bool contains(const _Ty & value)
	{
		return (set<_Ty>::vector(value) != vector<_Ty>::end());
	}
*/
};
/*
namespace parserInternal
{
template <typename _Ty>
	inline void parseval(const char * start,vectorEx<_Ty> & value,size_t & len)
	{

	};
}
*/


/*********************************************
	Allows C++ 2011 range based revers interation through a map

	Combining this with the Pair definition allows:

	for (<className>::Pair revIt : reverse_iterate(<classInstance>))
	{

	}

*/


template <typename T>
struct reverse_range
{
private:
  T& x_;
 
public:
  reverse_range (T& x): x_ (x) {}
 
  auto begin () const -> decltype (x_.rbegin ())
  {
    return x_.rbegin ();
  }
 
  auto end () const -> decltype (x_.rend ())
  {
    return x_.rend ();
  }
};
 
template <typename T>
reverse_range<T> reverse_iterate (T& x)
{
  return reverse_range<T> (x);
};

#endif
