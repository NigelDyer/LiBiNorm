// ***************************************************************************
// nelderMeadOptimiser.h (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// An implementation of the NelderMead optimser.  This code is based on code 
//  originally given to me in 1981 by John Cook at the British Telecom 
//  Research Labs, Martlesham Heath.  It was originally in HP basic and was
//	subsequently converted to Pascal and then C++.  It has been used for various
//	jobs over the years including the design of the line card on the System X
//  telephone exchange
// ***************************************************************************

#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H
#include <vector>
#include "stringEx.h"
#include "dataVec.h"

//	Uncomment to add an adhoc tweak that was found to help a problem some years ago but 
//	appears to degrade the performance of LiBiNorm 
//#define NOISY_COMBINE


inline VEC_DATA_TYPE sqr(VEC_DATA_TYPE a) { return a * a; };

//	Holds the information relating to an item to be optimised as a pair
//	The 'first' item of the pair is the value being iptimised, and the second 
//	indicates whether optimisation is currently enabled
class optiItem
{
public:
	VEC_DATA_TYPE value;
	bool optimise;
	//	Constructor.  Pass in the initial value and whether it is beinh optimised
	optiItem(VEC_DATA_TYPE v, bool o) :value(v), optimise(o) {};
	optiItem() :value(0), optimise(false) {};

	void readFloat(const char * buffer) {
		value = (VEC_DATA_TYPE)atof(buffer);
		const char * p = strchr(buffer, ',');
		if (p)
			optimise = (*(p + 1) == 'Y');
	};

};

//	A vector of values being optimised
struct optiVector : public std::vector<optiItem>
{
	optiVector(size_t i = 0, optiItem oi = optiItem()) :std::vector<optiItem>(i, oi) {};
	//	Return true if any of the entries are currently being optimised
	bool optimising()
	{
		for (size_t i = 0; i < size(); i++)
			if (at(i).optimise)
				return true;
		return false;
	};
	void setOptimising(bool mode = true);
};

//	A set of such vectors
typedef std::vector<optiVector *> optiDataType;

//	A class to contain both the error value representing the basic difference
//	between the target and the model and the value that includes a prior, e.g.
//	additional weightings to prevent negative coefficients etc.
struct ErrorPair {
	VEC_DATA_TYPE WithPrior, NoWeightings;

	ErrorPair(VEC_DATA_TYPE err = 0.0) :WithPrior(err), NoWeightings(err) {};
	void operator = (const ErrorPair & val) {
		WithPrior = val.WithPrior;
		NoWeightings = val.NoWeightings;
	};
	//	Cast the duplet to the value with prior by default as this is the value
	//	used during optimisation
	operator VEC_DATA_TYPE() const { return WithPrior; };
};


//	Internal class.  A consolidated set of the items being optimised
class optiVec : public dataVec
{
public:
	optiVec(size_t size = 0) :dataVec(size), m_error(0.0) {	};

	//	Dont allow a copy constructor
	optiVec(const optiVec & data) = delete;
	//	But do allow an r-value constructor
	optiVec(optiVec && data) { swap(data); };

	//	For setting and returning the error value associated with this set of parameters
	//	This returns the value with prior
	ErrorPair Error() { return m_error; };
	void SetError(ErrorPair error) { m_error = error; };

	//	Some useful operators
	optiVec & operator = (const optiVec & data);
	optiVec & operator = (optiVec && data);
	optiVec & operator += (const optiVec & data);
	optiVec & operator -= (const optiVec & data);

#ifdef NOISY_COMBINE
	// A noisy combination of two parameter vectors had previously appeared to be a useful addition
	//	to the algorithm, but this is no longer being used
	optiVec & noisyCombine(const optiVec & A, const optiVec & B, const optiVec & C);
#endif

private:
	ErrorPair m_error;

};

inline std::ostream& operator<< (std::ostream &out, const optiVec & d)
{
	for (size_t i = 0; i < d.size(); i++)
	{
		out << ", " << d[i];
	}
	return out;
}


class optiArray : public std::vector<optiVec>
{
public:
	optiArray(size_t size = 0);

	//	void resize(size_t size);
	size_t initialise(const optiDataType & data);
};


class nelderMeadOptimiser
{
public:
	nelderMeadOptimiser(int & iterations,const char * rootDir = "");

	void InitData(VEC_DATA_TYPE delta);
	void CreateAndScoreOrthogonalVectors(VEC_DATA_TYPE delta);

	VEC_DATA_TYPE optimise(optiDataType & data, int maxIter, int outputEvery,
		VEC_DATA_TYPE delta,stringEx id);

	void loadExternalData(optiVec & data);

	void GetError(optiVec & data);

	//	To be implemented by the class representing the data to be optimised.
	//	This generates the error 
	virtual ErrorPair ErrorFunc() = 0;
	//	Can be used to print results to screen/file during the process.
	virtual void SaveResults(bool toFile) = 0;

	template <typename... Params>
	void SaveProgress(const std::string &,const Params & ...params);

protected:
	std::string m_rootDir;

	//	The lowest error so far and the associated parameters
	ErrorPair m_ErrLowest;
	optiVec m_DataLowest;

	//	The last error
	ErrorPair m_ErrLast;

public:
	optiDataType * m_externalData;

	int & m_totalIterations;
	int m_recalcCount;
	int m_outputCount;

	size_t m_size;
	optiArray m_data;

	int	m_maxIter;
	int m_outputEvery;

	VEC_DATA_TYPE 	Error_highest, Error_second_highest, Eps, Err;
	size_t   Err_pointer_highest, Err_pointer_sec_highest, Err_pointer_lowest;

};

#endif // !NELDER_MEDA_H
