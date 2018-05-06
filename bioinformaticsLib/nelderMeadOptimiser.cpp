// ***************************************************************************
// nelderMeadOptimiser.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 11 August 2017
// ---------------------------------------------------------------------------
// An implementation of the NelderMead optimser.  This code is based on code 
//  originally given to me in 1981 by John Cook at the British Telecom 
//  Research Labs, Martlesham Heath.  It was originally in HP basic and was
//	subsequently converted to Pascal and then C++.  It has been used for various
//	jobs over the years including the design of the line card on the System X
//  telephone exchange
// ***************************************************************************

#include <algorithm>
#include "libCommon.h"
#include "nelderMeadOptimiser.h"

//	Uncomment to create a progress file
//#define SAVE_PROGRESS


#define  Al  1
#define  Ga  2
#define  Be  0.5		//Original


using namespace std;

//	Indicate that all of the paramateres are being optimised
void optiVector::setOptimising(bool mode)
{
	for (size_t i = 0; i < size(); i++)
		at(i).optimise = mode;
}

VEC_DATA_TYPE cumulateDifference(const dataVec & data1, const dataVec & data2)
{
	VEC_DATA_TYPE difference = 0;
	for (size_t i = 0; i < min<size_t>(data1.size(), data2.size()); i++)
	{
		difference += fabs(data1[i] - data2[i]);
	}
	return difference;
}


optiVec & optiVec::operator = (const optiVec & data)
{
	dataVec::operator =(data);
	m_error = data.m_error;
	return (*this);
}
optiVec & optiVec::operator = (optiVec && data)
{
	swap(data);
	return (*this);
}

optiVec & optiVec::operator += (const optiVec & data)
{
	for (size_t i = 0; i < size(); i++)
		at(i) += data[i];
	m_error = 0;
	return (*this);
}

optiVec & optiVec::operator -= (const optiVec & data)
{
	for (size_t i = 0; i < size(); i++)
		at(i) -= data[i];
	m_error = 0;
	return (*this);
}

optiVec operator * (const optiVec & data, VEC_DATA_TYPE multiplier)
{
	optiVec dv(data.size());
	for (size_t i = 0; i < data.size(); i++)
		dv[i] = data[i] * multiplier;
	return dv;
}

optiVec operator / (const optiVec & data, VEC_DATA_TYPE divider)
{
	optiVec dv(data.size());
	for (size_t i = 0; i < data.size(); i++)
		dv[i] = data[i] / divider;
	return dv;
}

optiVec operator + (const optiVec & data1, const optiVec & data2)
{
	optiVec dv(data1.size());
	for (size_t i = 0; i < data1.size(); i++)
		dv[i] = data1[i] + data2[i];
	return dv;
}

//	A set of rValue operations for efficiency
optiVec operator - (optiVec && data1, const optiVec & data2)
{
	for (size_t i = 0; i < data1.size(); i++)
		data1[i] -= data2[i];
	return move(data1);
}

optiVec operator / (optiVec && data, VEC_DATA_TYPE divider)
{
	for (auto & i : data)
		i /= divider;
	return move(data);
}

optiVec operator + (optiVec && data1, const optiVec & data2)
{
	for (size_t i = 0; i < data1.size(); i++)
		data1[i] += data2[i];
	return move(data1);
}

optiVec operator - (const optiVec & data1, const optiVec & data2)
{
	optiVec dv(data1.size());

	for (size_t i = 0; i < data1.size(); i++)
		dv[i] = data1[i] - data2[i];
	return dv;
}

#ifdef NOISY_COMBINE
optiVec & optiVec::noisyCombine(const optiVec & A, const optiVec & B, const optiVec & C)
{
	//	Make the choice as to which one to use based on bits returned from
	//	the randomise function
	int j = 1;
	int R = rand();
	for (size_t i = 0; i < size(); i++)
	{
		if (R & j)
			at(i) = C[i];
		else
			at(i) = A[i];
		j *= 2;
		if (j > 65536)
		{
			j = 1;
			R = rand();
		}
	}
	return (*this);
}
#endif

optiArray::optiArray(size_t size) : vector<optiVec>(size + 1)
{
}

size_t optiArray::initialise(const optiDataType & data)
{
	//	A square array, plus one extra row for the reference 'M' row
	resize(1);
	at(0).clear();

	//	Setup the first vector in the array with the entries in the
	//	input data that are to be optimised
	//	The data forms a two dimensional array, so go through all the entries picking out the ones
	//	that are set o be optimised
	for (auto & i : data)
	{
		for (auto j : *i)
			if (j.optimise)
				at(0).push_back(j.value);
	}

	//And then resize the master array so that there are N+1 entries
	resize(at(0).size() + 1);

	//	And return the number of entries
	return at(0).size();
}



nelderMeadOptimiser::nelderMeadOptimiser(int & iterations,const char * rootDir) :
	m_rootDir(rootDir),
m_externalData(nullptr),
m_totalIterations(iterations), m_recalcCount(0), m_outputCount(0), Eps(0)
{
}

void nelderMeadOptimiser::loadExternalData(optiVec & data)
{
	//	Transfers from the internal arrays holding the current values of the parameters being
	//	optimised to the external data values that are used within the model
	size_t p = 0;
	/*	for (size_t i = 0; i < (*m_externalData).size(); i++)
		{
			for (size_t j = 0; j < m_externalData ->at(i)->size(); j++)
				if (m_externalData->at(i)->at(j).optimise)
					m_externalData->at(i)->at(j).value = data[p++];
		}
		*/
	for (auto & i : *m_externalData)
	{
		for (auto & j : *i)
			if (j.optimise)
				j.value = data[p++];
	}
}

//
//	The starting or restarting point for the nelder mead optimisation is N+1 vectors 
//	which are an vector corresponding to the starting point for the optimisation
//	and a set of vectors where each parameter in turn is changed by delta 
void nelderMeadOptimiser::CreateAndScoreOrthogonalVectors(VEC_DATA_TYPE delta)
{
	size_t i;

	//	We have N+1 copies of the orginal data
	for (i = 1; i <= m_size; i++)
		m_data[i] = m_data[0];

//	Add the increments to create the N orthogonal vectors, with the O'th vector
//	being the original vector
	for (i = 0; i < m_size; i++)
	{
		if (m_data[i + 1][i] > delta)
			m_data[i + 1][i] -= delta;
		//		else if (m_data[i+1][i] == 0)
		//			m_data[i+1][i] = dx[i];
		else
			m_data[i + 1][i] += delta;
	}


	//	And now find the errors for each vector, 
	for (i = 0; i <= m_size; i++)
	{
		GetError(m_data[i]);
		SaveResults(m_data[i].Error() == m_ErrLowest);
	}
}


void nelderMeadOptimiser::InitData(VEC_DATA_TYPE delta)
{
	size_t i;
	m_size = m_data.initialise(*m_externalData);

	VEC_DATA_TYPE rmsTot = 0;
	for (i = 0; i < m_size; i++)
		rmsTot += sqr(m_data[0][i]);

	rmsTot = sqrt(rmsTot / m_size);

	m_recalcCount = m_totalIterations = m_outputCount = -(int)m_size;
	m_ErrLowest = MAX_DOUBLE;

	CreateAndScoreOrthogonalVectors(delta);

}

void nelderMeadOptimiser::GetError(optiVec & data)
{
	m_totalIterations++;
	m_recalcCount++;
	m_outputCount++;

	loadExternalData(data);
	m_ErrLast = ErrorFunc();
	data.SetError(m_ErrLast);

	if (m_ErrLast < m_ErrLowest)
	{
		m_ErrLowest = m_ErrLast;
		m_DataLowest = data;
	}

}

#define SLOPE 0.05
#define THRESHOLD 0.01

//	The Nelder mead engine itself
VEC_DATA_TYPE nelderMeadOptimiser::optimise(optiDataType & data, int maxIter, int outputEvery, 
	VEC_DATA_TYPE delta, stringEx id)
{
	int average_count = 0;
	VEC_DATA_TYPE last_err = MAX_DOUBLE;

	Err_pointer_lowest = 0;
	Err_pointer_highest = 0;
	Error_second_highest = 0;

	//	Pointers to the locations of the external data being optimised
	m_externalData = &data;

	InitData(delta);

	bool Optiexit = false;
	bool averaging = true;

	optiVec aveVector(m_size);
	optiVec leadVector(m_size);
	optiVec trialVector(m_size);

	m_DataLowest.resize(m_size);

	do
	{
		Error_highest = 0;
		Error_second_highest = 0;
		m_ErrLowest = MAX_DOUBLE;

		size_t I;
		for (I = 0; I <= m_size; I++)
		{
			if (m_data[I].Error() > Error_highest)
			{
				Error_second_highest = Error_highest;
				Err_pointer_sec_highest = Err_pointer_highest;
				Error_highest = m_data[I].Error();
				Err_pointer_highest = I;
			}
			else if (m_data[I].Error() > Error_second_highest)
			{
				Error_second_highest = m_data[I].Error();
				Err_pointer_sec_highest = I;
			};
			if (m_data[I].Error() < m_ErrLowest)
			{
				m_ErrLowest = m_data[I].Error();
				Err_pointer_lowest = I;
			}
		};

		// Check for exit.}
		if (m_totalIterations >= (int)maxIter)
			Optiexit = true;

		//	Pointers to the original data.  
		m_DataLowest = m_data[Err_pointer_lowest];
		const optiVec & worstVector = m_data[Err_pointer_highest];

		Err = m_ErrLowest;

		if (Eps == 0)
			Eps = m_ErrLowest * 2;
		else
			Eps = (Eps * (1-SLOPE)) + m_ErrLowest * SLOPE;

		// Check for exit.}
		if ((m_ErrLowest + THRESHOLD) > Eps)
		{
			SaveProgress("Finishing early after ", m_totalIterations, " iterations");
			Optiexit = true;
		}
		else
		{
			//	Xo is the average of all the vectors except the worst
			aveVector = m_data[0];
			for (I = 1; I <= m_size; I++)
				aveVector += m_data[I];
			aveVector = (aveVector - worstVector) / m_size;

			//	aveVector is the average, plus a bit pointing away from the worst result
			//	ie downhill is in the opposite direction from uphill
			//  
			//	Default Al = 2
			leadVector = (aveVector * (1 + Al)) - (worstVector * Al);
			GetError(leadVector);

			//	GetError will have set the error to the lowest if it is a new low
			if (leadVector.Error() == m_ErrLowest)
			{
				SaveResults(true);
				//	New vector looks better, try a second that is the same more in the same direction
				trialVector = (leadVector * Ga) + (aveVector * (1 - Ga));	//Ga = 2
				GetError(trialVector);

				//	Pick the one that gave the best results, and replace the worst with it
				if (trialVector.Error() == m_ErrLowest)
				{
					SaveResults(true);

//					if (debugPrint) cerr << ".";

#ifdef NOISY_COMBINE
					optiVec noisyTrialVector(m_size);
					noisyTrialVector.noisyCombine(aveVector, leadVector, trialVector);
					GetError(noisyTrialVector);
					if (noisyTrialVector.Error() == m_ErrLowest)
					{
						//	Noisy Double better than the lowest
						m_data[Err_pointer_highest] = noisyTrialVector;
						SaveResults(true);
					}
					else
#endif
					{
						//	Double better than the lowest
						m_data[Err_pointer_highest] = trialVector;
					}
				}
				else
				{
					//	Better than the lowest
					m_data[Err_pointer_highest] = leadVector;
				}
			}
			else
			{
				//	The new one we have found is better than the second worst, so use it to replace
				//	the worst, (which it must also be better than, by definition
				if (leadVector.Error() <= Error_second_highest)
				{
					//	Better than the second highest
					m_data[Err_pointer_highest] = leadVector;
				}
				else
				{
					//	The new one is only better than the worst, so replace it
					if (leadVector.Error() < Error_highest)
						m_data[Err_pointer_highest] = leadVector;
	
					//	Trial vector which is half the worst plus half the average
					trialVector = (worstVector * Be) + (aveVector * (1 - Be));
					GetError(trialVector);

					//	If the average is worse than the highest then take more drastic action
					if (trialVector.Error() > Error_highest)
					{
						if ((last_err / m_ErrLowest) < 1.000001)
							average_count++;
						else
							average_count = 0;
						debugMessage(id,": Err = ", Err, ", Average count = ",average_count);
						last_err = m_ErrLowest;
						if (average_count > 50)
							Optiexit = true;

						Err = m_ErrLowest;

						//{ The second term does not make much sense, but was there in the original}
						//	Anyway, Eps is zero, so it had no effect anyway
//						if ((Err>Eps)) //|| (Err<Mx) )
						{
							//	Swap between two strategies for restarting if we get stuck.
							if (averaging)
							{
								//	Change all the vectors to be midway between the current value 
								//	and the lowest.  For the row that is the lowest, the values will
								//	end up being (essentially) unchanged
								debugMessage(id,": Err = ", Err, ", Averaging");
								for (I = 0; I <= m_size; I++)
								{	
									m_data[I] = (m_data[I] * 0.5) + (m_DataLowest * 0.5);
									GetError(m_data[I]);
									if (m_data[I].Error() == m_ErrLowest)
										SaveResults(true);

								};
							}
							else
							{
								debugMessage(id,": Err = ", Err, ", Reorthogonalising, Delta = ",delta);
								m_data[0] = m_DataLowest;
								CreateAndScoreOrthogonalVectors(delta);
							}
							averaging = !averaging;
						}
//						else
//							Optiexit = true;
					}
					else
					{
						m_data[Err_pointer_highest] = trialVector;
					};
				};
			};
		};
		if (m_outputCount > outputEvery)
		{
			m_outputCount = 0;
			loadExternalData(m_DataLowest);
			SaveProgress(id, m_DataLowest);
//			SaveResults(false);
		}
	} while (!Optiexit);

	loadExternalData(m_DataLowest);
	SaveResults(false);

	return m_ErrLowest.NoWeightings;
}
