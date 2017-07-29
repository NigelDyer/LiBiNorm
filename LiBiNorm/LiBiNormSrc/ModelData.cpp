// ***************************************************************************
// ModelData.cpp (c) 2017 Nigel Dyer
// School of Life Sciences, University of Warwick
// ---------------------------------------------------------------------------
// Last modified: 24 July 2017
// ---------------------------------------------------------------------------
// Code for each of the five models
// ***************************************************************************

#include "rand.h"
#include "stringEx.h"
#include "containerEx.h"
#include "ModelData.h"
#include "Options.h"
using namespace std;

//	Sets the methods to be used for the liklihood (sum of squares) and prior functions
void setSSfun(optionsType & options,modelType m)
{
	switch (m)
	{
	case noModelSpecified: break;
	case ModelA:
		options.ssfun = &FLL_ModelA;
		options.priorfun = &priorFunc;
		break;
	case ModelB:
		options.ssfun = &FLL_ModelB;
		options.priorfun = &priorFunc;
		break;
	case ModelC:
		options.ssfun = &FLL_ModelC;
		options.priorfun = &priorFunc;
		break;
	case ModelD:
		options.ssfun = &FLL_ModelD;
		options.priorfun = &priorFunc;
		break;
	case ModelE:
		options.ssfun = &FLL_ModelE;
		options.priorfun = &priorFuncE;
		break;
	case ModelBD:
		options.ssfun = &FLL_ModelBD;
		options.priorfun = &priorFunc;
		break;
	case findBestModel:
			break;
	}
}



//	Returns a list of all the models, which is used to iterate through the list
const std::vector<modelType> & allModels()
{
	static std::vector<modelType> list{ ModelA ,ModelB ,ModelC,ModelD,ModelE,ModelBD };
	return list;
};

//	Converts a modelType to a string
string conv(const modelType m, bool removeGaps)
{
	stringEx retVal;
	switch (m)
	{
	case noModelSpecified: retVal = "No model specified"; break;
	case ModelA: retVal = "Model A"; break;
	case ModelB: retVal = "Model B"; break;
	case ModelC: retVal = "Model C"; break;
	case ModelD: retVal = "Model D"; break;
	case ModelE: retVal = "Model E"; break;
	case ModelBD: retVal = "Model BD"; break;
	case findBestModel: retVal = "Best"; break;
	case none: retVal = "None"; break;
	}
	if (removeGaps)
		retVal.replace(" ", "");
	return retVal;
}

//	The following are needed by the stringEx class, but are useful in other contexts
//	Used by microsoft stringEx
string std::to_string(const modelType & m)
{
	return conv(m);
}

//	Used by gcc stringEx
ostream& operator<< (ostream &out, const modelType & m)
{
	out << conv(m);
	return out;
}

//  Selects a model based on a string, exits if the string is not valid.  Used, for example, when 
//	parsing the command line
modelType modelFromString(const string & desc)
{
	static map<string, modelType> mappings{
		{"A",ModelA },{"B",ModelB },{"C",ModelC },{"D",ModelD },{"E",ModelE },{"BD",ModelBD }, 
		{"a",ModelA },{"b",ModelB },{"c",ModelC },{"d",ModelD },{"e",ModelE },{"bd",ModelBD }, 
		{ "PolyA",ModelD },{"polya",ModelD },{ "POLYA",ModelD },{ "polyA",ModelD },
		{"random",ModelE},{ "RANDOM",ModelE },
		{ "smart",ModelBD },{ "SMART",ModelBD },{ "Smart",ModelBD },
		{"best",findBestModel },{"BEST",findBestModel },{ "Best",findBestModel },
		{"none",none},{"NONE",none},{ "None",none }
	};
	auto iter = mappings.find(desc);
	if (iter == mappings.end())
		exitFail("Unknown model description:", desc);
	return (*iter).second;
}
//	Used by the printTSV class to print modelTypes to files
bool printVal(outputDataFile * f, modelType m)
{
	fputs(conv(m).c_str(), f->fout);
	return true;
};

//	udentifier, min max and initial value for each of the parameters
#define PARAM_D { "d", -1, 2, -0.5 }
#define PARAM_H { "h", 0, 3, 1.5 }
#define PARAM_T1 { "t1", -5, -1, -3 } 
#define PARAM_T2 { "t2", -5, -1, -3 } 
#define PARAM_A { "a", 0, 1, 0.5 }


//	Loads the paremeter information associated with each of the models and sets the value to a random value within the allowed
//	range for the parameter. 
paramDescriptionSet GetModelParams(modelType model,dataVec * defaults, VEC_DATA_TYPE offset)
{
	paramDescriptionSet params;

	switch (model)
	{
	case noModelSpecified:
		break;
	case ModelB: case ModelD: case ModelE:
		params = { PARAM_D   // average length of fragments
			, PARAM_H   // the minimum length of fragmenation
			, PARAM_T1   // theta1
			, PARAM_T2 // theta2
		};
		break;
	case ModelC:
		params = { PARAM_D    // average length of fragments
			, PARAM_H  // the minimum length of fragmenation
			, PARAM_T2 // theta2
		};
		break;
	case ModelA:
		params = { PARAM_D    // average length of fragments
			, PARAM_H   // the minimum length of fragmenation
		};
		break;
	case ModelBD:
		params = { PARAM_D    // average length of fragments
			, PARAM_H   // the minimum length of fragmenation
			, PARAM_T1  // theta1
			, PARAM_T2 // theta2
			, PARAM_A // alpha strength of model B
		};
		break;
	case findBestModel:
	case none:
		break;
	};

#ifdef SIMULATE_MATLAB_BUG
	params[0].value = rand(3);
#endif

	if ((defaults) && (defaults->size()))
	{
		params.setValues(defaults->addNoise(offset));
//		params.setValues(*defaults);
	}
	else
	{
		dataVec r = randn(params.size());
		for (size_t i = 0; i < params.size(); i++)
			params[i].value = params[i].initial + r[i] * offset;
	}

	return params;
}

//	Returns the list of the headers using the data returned from GetModelParams
headerType getHeaders()
{
	headerType _retVal;
	for (modelType m : allModels())
	{
		paramDescriptionSet params = GetModelParams(m);
		for (auto i : params)
			_retVal[m].push_back(i.name);
	}
	return _retVal;
}

/* 
	Common Priors for the models
	This is used for all models and adds quadratic cost function for approaching too
	close to the extremes of the parameter range which is set as a prior

	The code also includes an option of setting a terget for the d and h parameters. This was found not 
	to be necessary once the Nelder Mead initialisation was performed in two stages.  The first being to find the values
	for t1 ,t2 and a, keeping d and h fixed, and then to look for d and h. The code has been retained in case
	it turns out to be of use in the future
*/

double priorFunc(const dataVec & data, const paramDescriptionSet & params)
{
	VEC_DATA_TYPE _retVal = 0;
	for (size_t i = 0; i < params.size(); i++)
	{
		VEC_DATA_TYPE diff = data[i] - params[i].max + 0.2;
		if (diff > 0)
			_retVal += (diff * diff) * EDGE_PENALTY_MULTIPLIER;
		else
		{
			diff = params[i].min - data[i] + 0.2;
			if (diff > 0)
				_retVal += (diff * diff) * EDGE_PENALTY_MULTIPLIER;
		}
	}

	//	Priors for d and h paremeters
#ifdef D_PRIOR_MULTIPLIER
	static double d_target = log10(D_PRIOR_TARGET);
	VEC_DATA_TYPE diff = abs(data[0] - d_target);
	//	Quadratic until the difference is 1 and then linear with matching slope
	if (diff < 1)
		_retVal += (diff * diff * D_PRIOR_MULTIPLIER);
	else
		_retVal += (2 * diff - 1) * D_PRIOR_MULTIPLIER;
#endif
#ifdef H_PRIOR_MULTIPLIER
	static double h_target = log10(H_PRIOR_TARGET);
	diff = abs(data[1] - h_target);
	if (diff < 1)
		_retVal += (diff * diff * H_PRIOR_MULTIPLIER);
	else
		_retVal += (2 * diff - 1) * H_PRIOR_MULTIPLIER;
#endif

	return _retVal;
};

/*
	This adds a target prior just for the h parameter in model E.  The cost function is weighted using a 
	quadratic for differences in the log value up to 1 away from the target (ie a factor of 10 in absolute terms)
	and then linearly beyond that
*/

double priorFuncE(const dataVec & data, const paramDescriptionSet & params)
{
	VEC_DATA_TYPE _retVal = priorFunc(data, params);
#ifdef E_H_PRIOR_MULTIPLIER
	{
		static double h_target = log10(E_H_PRIOR_TARGET);
		VEC_DATA_TYPE diff = abs(data[1] - h_target);
		if (diff < 1)
			_retVal += (diff * diff * E_H_PRIOR_MULTIPLIER);
		else
			_retVal += (2 * diff - 1) * E_H_PRIOR_MULTIPLIER;
		//	_retVal += data[1] * E_H_PRIOR_MULTIPLIER;
	}
#endif
#ifdef E_D_PRIOR_MULTIPLIER
	{
		static double d_target = log10(E_D_PRIOR_TARGET);
		VEC_DATA_TYPE diff = abs(data[1] - d_target);
		if (diff < 1)
			_retVal += (diff * diff * E_D_PRIOR_MULTIPLIER);
		else
			_retVal += (2 * diff - 1) * E_D_PRIOR_MULTIPLIER;
	}
#endif
	return _retVal;

}

/*

The following code implements the six different models for the bias within an RNA transcript.

In each case the parameters (between 2 and 5 ) are passed in with param and the data,
ie the information about the reads and the genes are passed in in data.

data.fragData contains the information about the individual reads, one entry per read
data.geneIndex indeicates which gene a read is associated with, one entry per read
data.geneData gives inforation about each gene
data.geneLengths is the lengths of the genes
data.geneFrequencies are the frequencies

One modification from the original matlab code arises from the fact that the normalisation
values are the same for all reads in the gene so only need to be calculated per gene.

The original MATLAB code is shown in comments
*/


//	The log liklyhood calculations for each of the models
double FLL_ModelA(const dataVec & param, const mcmcGeneData & data)
{
	double d = pow(10,param[0]);
	double h = pow(10,param[1]);

	const vector<int> & geneIndex = data.geneIndex;
	double LogL = 0;

	//	f_frag = (x> h).*(x < l-h) + 1/d;
	//	norm =  (2*h<l).*(l-2*h) + l/d;
	//LogL = -2*sum(log(f_frag./norm)./freq_l);
#ifdef VECTOR_MATHS
	const dataVec & x = data.fragData;
	const dataVec l = data.geneLengths.expand(geneIndex);
	const dataVec freq_l = data.geneFrequencies.expand(geneIndex);
	
	dataVec f_frag = (x> h)*(x < l-h) + 1/d;
	dataVec norm =  (2*h<l)*(l-2*h) + l/d;
	LogL = sum(log(f_frag/norm)/freq_l);
#else
	double last_l = 0;
	double norm=0;
	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const VEC_DATA_TYPE & x = data.fragData[i];
		const VEC_DATA_TYPE & l = data.geneLengths[geneIndex[i]];
		const VEC_DATA_TYPE & freq_l = data.geneFrequencies[geneIndex[i]];

		double f_frag = 1.0/d;
		if ((x> h) && (x < (l-h))) f_frag += 1.0;

		if (f_frag == 0) return 1E20;

		if (l != last_l)
		{
			norm =  l/d;
			if (2*h<l)
				norm +=  (l-(2*h));
			last_l = l;
		}
		LogL += log(f_frag/norm)/freq_l;
	}
#endif
	LogL = -2*LogL;
	return LogL;
}


double FLL_ModelB(const dataVec & param, const mcmcGeneData & data)
{

	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t1 = pow(10,param[2]);
	double t2 = pow(10,param[3]);

	const vector<int> & geneIndex = data.geneIndex;
	double LogL = 0;


//	f_frag = (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) + ...
//         1./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d;
//	 norm =  (2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l*(t1+t2)))/(t1 + t2)^2 + ...
//        (exp(-l*(t1+t2)).*(l.*t2^2+l.*t2*t1-t1)+t1)/(t1 + t2)^2/d;
#ifdef VECTOR_MATHS
	const dataVec & x = data.fragData;
	const dataVec l = data.geneLengths.expand(geneIndex);
	const dataVec freq_l = data.geneFrequencies.expand(geneIndex);
	dataVec f_frag = (x> h)*(x < l-h)/(t1+ t2) * (t1*exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) + 
	       1./(t1+ t2) * (t1*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d;
	dataVec norm =  (2*h<l)*(t1*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2)*(l-2*h)*exp(-l*(t1+t2)))/((t1 + t2) * (t1 + t2)) +
	        (exp(-l*(t1+t2))*((l*t2) *(l*t2)+l*t2*t1-t1)+t1)/(t1 + t2)^2/d;
	LogL = sum(log(f_frag / norm) / freq_l);
#else
	double last_l = 0;
	double norm = 0;
	double t1_p_t2 = t1+t2;
	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const VEC_DATA_TYPE & x = data.fragData[i];
		const VEC_DATA_TYPE & l = data.geneLengths[geneIndex[i]];
		const VEC_DATA_TYPE & freq_l = data.geneFrequencies[geneIndex[i]];

		double exp_l_t1_t2 = exp(-l * t1_p_t2);
		double t2_exp_l_t1_t2_full = t2*exp_l_t1_t2;
		double t1_exp_c_l_t1_t2 = t1*exp((x-l)*t1_p_t2);

		double f_frag = (t1_exp_c_l_t1_t2 + t2_exp_l_t1_t2_full)/d;
		if ((x> h) && (x < (l-h)))
			f_frag += (t1_exp_c_l_t1_t2/exp(h*t1_p_t2) + t2_exp_l_t1_t2_full);

		f_frag /= t1_p_t2;

		if (f_frag == 0)
			return 1E20;

		if (l != last_l)
		{
			norm = ((exp_l_t1_t2*(l*t2*t1_p_t2-t1)+t1)/d);
			if ((2*h)<l)
				norm += (t1*(exp(-2*h*t1_p_t2)-exp_l_t1_t2)+t2*t1_p_t2*(l-(2*h))*exp_l_t1_t2);
			norm /= (t1_p_t2*t1_p_t2);

			last_l = l;
		}

		LogL += log(f_frag/norm)/freq_l;
	}

#endif
	LogL = -2 * LogL;
	return LogL;
}

double FLL_ModelC(const dataVec & param, const mcmcGeneData & data)
{
	//	function [LogL] = FLL_Deng(param, data)

	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t2 = pow(10,param[2]);

	const vector<int> & geneIndex = data.geneIndex;

	double last_l = 0;
	double norm=0;
	double LogL = 0;

//    f_frag = (x> h).*(x < l-h).*exp(-t2*(x+h)) + ...
//        exp(-t2*(x))/d;
//    norm =  (2*h<l).*(exp(-2*h*t2) - exp(-l*( t2)))/t2 + ...
//        (1-exp(-l*( t2)))/t2/d;
	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const VEC_DATA_TYPE & x = data.fragData[i];
		const VEC_DATA_TYPE & l = data.geneLengths[geneIndex[i]];
		const VEC_DATA_TYPE & freq_l = data.geneFrequencies[geneIndex[i]];

		double f_frag = exp(-t2*(x))/d;
		if ((x> h) && (x < l-h))
			f_frag += exp(-t2*(x+h));

		if (f_frag == 0)
			return 1E20;


		if (l != last_l)
		{
			double exp_ml_t2 = exp(-l*( t2));

			norm = (1-exp_ml_t2)/t2/d;
			if (2*h<l)
				norm += (exp(-2*h*t2) - exp_ml_t2)/t2;
			last_l = l;
		}

//		LogL = -2*sum(log(f_frag/norm(geneIndex))/freq_l(geneIndex));
		LogL += log(f_frag/norm)/freq_l;
	}

	LogL = -2 * LogL;
	return LogL;
}


double FLL_ModelD(const dataVec & param, const mcmcGeneData & data)
{
	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t1 = pow(10,param[2]);
	double t2 = pow(10,param[3]);
	double LogL = 0;

	const vector<int> & geneIndex = data.geneIndex;

	double last_l = 0;
	double norm=0;
	double temp_l = 0;

//f_frag = (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
//         1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d;
	double t1_p_t2 = t1+t2;

	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const VEC_DATA_TYPE & x = data.fragData[i];
		const VEC_DATA_TYPE & l = data.geneLengths[geneIndex[i]];
		const VEC_DATA_TYPE & freq_l = data.geneFrequencies[geneIndex[i]];

		double f_frag = 1/t1_p_t2 * (t1*exp(-t1*(l-x)) + t2*exp(-t1*l-t2*(x)))/d;

		if ((x > h) && (x < l-h))
			f_frag += (t1*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2*exp(-t1*l-t2*(x+h)))/t1_p_t2;

		if (f_frag == 0)
			return 1E20;

		if (l != last_l)
		{

//    norm =  (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + ...
//    (1-exp(-l*(t1 + t2)))/(t1 + t2)/d;
			double exp_ml_t1_p_t2 = exp(-l*t1_p_t2);
			norm =  (2*h<l)*(exp(-2*h*t1_p_t2) - exp_ml_t1_p_t2)/t1_p_t2 + 
				(1-exp_ml_t1_p_t2)/t1_p_t2/d;

			last_l = l;


//			LogL = -2*sum(log(f_frag/norm(geneIndex))/Freq_l(geneIndex));
		}
		temp_l += log(f_frag/norm)/freq_l;
	}
	LogL = -2*temp_l;
	return LogL;
}

double FLL_ModelE(const dataVec & param, const mcmcGeneData & data)
{
	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t1 = pow(10,param[2]);
	double t2 = pow(10,param[3]);

	const vector<int> & geneIndex = data.geneIndex;
	double LogL = 0;

	//	f_frag = (x> h).*(x < l-h)./t1/(t1+ t2).*(exp(-2*h*(t1+ t2))-exp(-(h +x)*(t1+t2)) - exp(-t1*h-2*h*t2-(l-x)*t1) + exp(-h*t2-l*t1-x*t2)) + ...
	//     1./t1/(t1+ t2) .*(1 - exp(-x*(t1+t2)) - exp(-(l-x)*t1) + exp(-l*t1-x*t2))/d;
	//	norm =  (2*h<l).*(exp(-l*t1 - 2*h*t2)*(t1 + t2)^2 - exp(-l*(t1 + t2))*t1^2 + t1*t2*exp(-2*h*(t1 + t2))*(l*t2 -2*h*t1 -2*h*t2+l*t1 - t2/t1 - 2))/(t1 + t2)^2/t1^2/t2 + ...
	//     (l-1/(t1 + t2) - 1/t1 - t1/t2/(t1+t2)*exp(-l*(t1 + t2))+(t1 + t2)/t1/t2*exp(-l*t1))/(t1 + t2)/t1/d;
#ifdef VECTOR_MATHS
		const dataVec & x = data.fragData;
		const dataVec l = data.geneLengths.expand(geneIndex);
		const dataVec freq_l = data.geneFrequencies.expand(geneIndex);
		dataVec f_frag = (x > h)*(x < l - h) / t1 / (t1 + t2)*(exp(-2 * h*(t1 + t2)) - exp(-(x + h)*(t1 + t2)) - exp(-t1*h - 2 * h*t2 - (l - x)*t1) + exp(-h*t2 - l*t1 - x*t2)) +
			1 / t1 / (t1 + t2) *(1 - exp(-x*(t1 + t2)) - exp(-(l - x)*t1) + exp(-l*t1 - x*t2)) / d;
		dataVec	norm = (2 * h < l)*(exp(-l*t1 - 2 * h*t2)*(t1 + t2)*(t1 + t2) - exp(-l*(t1 + t2))*t1*t1 + t1*t2*exp(-2 * h*(t1 + t2))*(l*t2 - 2 * h*t1 - 2 * h*t2 + l*t1 - t2 / t1 - 2)) / ((t1 + t2)*(t1 + t2)) / (t1*t1) / t2 +
			(l - 1 / (t1 + t2) - 1 / t1 - t1 / t2 / (t1 + t2)*exp(-l*(t1 + t2)) + (t1 + t2) / t1 / t2*exp(-l*t1)) / (t1 + t2) / t1 / d;
	LogL = sum(log(f_frag / norm) / freq_l);
#else

		double last_l = 0;
		double norm = 0;

		double t1_p_t2 = t1 + t2;
		double t1_p_t2_sq = t1_p_t2*t1_p_t2;
		double exp_m2_h_t1_p_t2 = exp(-2 * h*t1_p_t2);

		for (size_t i = 0; i < data.fragData.size(); i++)
		{
			const VEC_DATA_TYPE & x = data.fragData[i];
			const VEC_DATA_TYPE & l = data.geneLengths[geneIndex[i]];
			const VEC_DATA_TYPE & freq_l = data.geneFrequencies[geneIndex[i]];

			double	f_frag = 1 / t1 / t1_p_t2 *(1 - exp(-x*t1_p_t2) - exp(-(l - x)*t1) + exp(-l*t1 - x*t2)) / d;

			if ((x > h) && (x < l - h))
				f_frag += (exp_m2_h_t1_p_t2 - exp(-(x + h)*t1_p_t2) - exp(-t1*h - 2 * h*t2 - (l - x)*t1) + exp(-h*t2 - l*t1 - x*t2)) / t1 / t1_p_t2;

			if (f_frag == 0)
				return 1E20;

			if (l != last_l)
			{
				norm = (l - 1 / t1_p_t2 - 1 / t1 - t1 / t2 / t1_p_t2*exp(-l*t1_p_t2) + t1_p_t2 / t1 / t2*exp(-l*t1)) / t1_p_t2 / t1 / d;

				if (2 * h < l)
					norm += (exp(-l*t1 - 2 * h*t2)*t1_p_t2_sq - exp(-l*t1_p_t2)*t1*t1 + t1*t2*exp(-2 * h*t1_p_t2)*(l*t2 - 2 * h*t1 - 2 * h*t2 + l*t1 - t2 / t1 - 2)) / (t1_p_t2*t1_p_t2*t1*t1*t2);

				last_l = l;
			}
		LogL += log(f_frag/norm)/freq_l;
		}
#endif
	LogL = -2 * LogL;
	return LogL;
	}

double FLL_ModelBD(const dataVec & param, const mcmcGeneData & data)
{

	double d = pow(10,param[0]);
	double h = pow(10,param[1]);
	double t1 = pow(10,param[2]);
	double t2 = pow(10,param[3]);
	double a = param[4];

	const vector<int> & geneIndex = data.geneIndex;

	double LogL = 0;
//	Original MATLAB code for reference
//	f_frag = a*( (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l-h+x)) + t2*exp(-l*(t1+t2))) + ...
//         1./(t1+ t2) .* (t1.*exp(-2*l*(t1+t2)+(t1+t2)*(l+x)) + t2*exp(-l*(t1+t2)))/d) + ...
//         (1-a)*( (x> h).*(x < l-h)./(t1+ t2) .* (t1.*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2.*exp(-t1*l-t2*(x+h))) + ...
//         1./(t1+ t2) .* (t1.*exp(-t1*(l-x)) + t2.*exp(-t1*l-t2*(x)))/d);
//	 norm = a*( (2*h<l).*(t1.*(exp(-2*h*(t1+t2))-exp(-l*(t1+t2)))+t2*(t1+t2).*(l-2*h).*exp(-l*(t1+t2)))/(t1 + t2)^2 + ...
//        (exp(-l.*(t1+t2)).*(l.*t2^2+l.*t2*t1-t1)+t1)/(t1 + t2)^2/d) + ...
//        (1-a)*( (2*h<l).*(exp(-2*h*(t1 + t2)) - exp(-l*(t1 + t2)))/(t1 + t2) + ...
//        (1-exp(-l.*(t1 + t2)))/(t1 + t2)/d);

#ifdef VECTOR_MATHS
//	Example code using full dataVec arithmatic.  Slower but closer to the original MATLAB code
	const dataVec & x = data.fragData;
	const dataVec l = data.geneLengths.expand(geneIndex);
	const dataVec freq_l = data.geneFrequencies.expand(geneIndex);

	dataVec f_frag = a*((x > h)*(x < l - h) / (t1 + t2) * (t1*exp(-2 * l*(t1 + t2) + (t1 + t2)*(l - h + x)) + t2*exp(-l*(t1 + t2))) +
		1 / (t1 + t2) * (t1*exp(-2 * l*(t1 + t2) + (t1 + t2)*(l + x)) + t2*exp(-l*(t1 + t2))) / d) +
		(1 - a)*((x > h)*(x < l - h) / (t1 + t2) * (t1*exp(-t1*(l - x) - 2 * t2*h - t1*h) + t2*exp(-t1*l - t2*(x + h))) +
			1 / (t1 + t2) * (t1*exp(-t1*(l - x)) + t2*exp(-t1*l - t2*(x))) / d);

	dataVec norm = a*((2 * h < l)*(t1*(exp(-2 * h*(t1 + t2)) - exp(-l*(t1 + t2))) + t2*(t1 + t2)*(l - 2 * h)*exp(-l*(t1 + t2))) / ((t1 + t2) * (t1 + t2)) +
		(exp(-l*(t1 + t2))*(l*(t2 *t2) + l*t2*t1 - t1) + t1) / ((t1 + t2) * (t1 + t2)) / d) +
		(1 - a)*((2 * h < l)*(exp(-2 * h*(t1 + t2)) - exp(-l*(t1 + t2))) / (t1 + t2) +
		(1 - exp(-l*(t1 + t2))) / (t1 + t2) / d);

	LogL = sum(log(f_frag / norm) / freq_l);
#else
	double last_l = 0;
	double norm = 0;

	double t1_p_t2 = t1+t2;
	double t1_p_t2_sq = t1_p_t2*t1_p_t2;

	for (size_t i = 0;i < data.fragData.size(); i++)
	{
		const VEC_DATA_TYPE & x = data.fragData[i];
		const VEC_DATA_TYPE & l = data.geneLengths[geneIndex[i]];
		const VEC_DATA_TYPE & freq_l = data.geneFrequencies[geneIndex[i]];

		double exp_ml_t1_p_t2 = exp(-l*t1_p_t2);
		double f_frag = 1/t1_p_t2 * (t1*exp(-2*l*t1_p_t2+t1_p_t2*(l+x)) + t2*exp_ml_t1_p_t2)/d;
		double f_frag_pt2 = (t1*exp(-t1*(l-x)) + t2*exp(-t1*l-t2*(x)))/t1_p_t2/d;

		if ((x > h) && (x < (l-h)))
		{
			f_frag += (t1 * exp(-2*l*t1_p_t2+t1_p_t2*(l-h+x)) + t2*exp_ml_t1_p_t2)/t1_p_t2;
			f_frag_pt2 += (t1*exp(-t1*(l-x) - 2*t2*h - t1*h) + t2*exp(-t1*l-t2*(x+h)))/t1_p_t2;
		}

		f_frag = a * f_frag + (1-a)*f_frag_pt2;

		if (f_frag == 0)
			return 1E20;

		if (l != last_l)
		{
			norm = (exp_ml_t1_p_t2*(l*t2*t2+l*t2*t1-t1)+t1)/t1_p_t2_sq/d;
			double norm_pt2 = (1-exp_ml_t1_p_t2)/t1_p_t2/d;

			if (2*h<l)
			{
				double exp_m2_h_t1_p_t2 = exp(-2*h*t1_p_t2);
				norm += (t1*(exp_m2_h_t1_p_t2-exp_ml_t1_p_t2)+t2*t1_p_t2*(l-2*h)*exp_ml_t1_p_t2)/t1_p_t2_sq;
				norm_pt2 += (exp_m2_h_t1_p_t2 - exp_ml_t1_p_t2)/t1_p_t2;
			}
			norm = a*norm +(1-a) * norm_pt2;
			last_l = l;
		}
		LogL += log(f_frag/norm)/freq_l;
	}

#endif
	LogL = -2 * LogL;
	return LogL;
}


void getBias(modelType m,dataVec & params,const dataVec & l, dataVec & bias)
{
	double d = pow(10, params[0]);
	double h = pow(10, params[1]);
	double t1, t2, a;
	switch (m)
	{
	case noModelSpecified:	//This should not happen
	case ModelA:
		break;
	case ModelB:
	case ModelD:
	case ModelE:
		t1 = pow(10, params[2]);
		t2 = pow(10, params[3]);
		break;
	case ModelC:
		t1 = 0;
		t2 = pow(10, params[2]);
		break;
	case ModelBD:
		t1 = pow(10, params[2]);
		t2 = pow(10, params[3]);
		a = params[4];
		break;
	case findBestModel:
		break;
	}

	bias.resize(l.size());

	switch (m)
	{
	case noModelSpecified:
		break;
	case ModelA:
		bias = (2 * h < l)*(l - 2 * h) + l / d;
		break;
	case ModelB:
		bias = ((2 * h < l)*(t1*(exp(-2 * h*(t1 + t2)) - exp(-l*(t1 + t2))) + t2*(t1 + t2)*(l - 2 * h)*exp(-l*(t1 + t2))) / ((t1 + t2)*(t1 + t2)) +
			exp(-l*(t1 + t2))*(l*t2*t2 + t1*(exp(l*(t1 + t2)) + l*t2 - 1)) / ((t1 + t2)*(t1 + t2)) / d);
		break;
	case ModelC:
	{
		dataVec exp_ml_t2 = exp(-l*(t2));
		bias = (2 * h < l)*(exp(-2 * h*t2) - exp_ml_t2) / t2 + (1 - exp_ml_t2) / t2 / d;
	}

	break;
	case ModelD:
		for (size_t i = 0; i < l.size(); i++)
		{
			if (2 * h < l[i])
				bias[i] = (exp(-2 * h*(t1 + t2)) - exp(-l[i] * (t1 + t2))) / (t1 + t2) + (1 - exp(-l[i] * (t1 + t2))) / (t1 + t2) / d;
			else
				bias[i] = (1 - exp(-l[i] * (t1 + t2))) / (t1 + t2) / d;
		}
		break;
	case ModelE:
		bias = (2 * h < l)*(exp(-l*t1 - 2 * h*t2)*(t1 + t2)*(t1 + t2) - exp(-l*(t1 + t2))*t1*t1 + t1*t2*exp(-2 * h*(t1 + t2))*(l*t2 - 2 * h*t1 - 2 * h*t2 + l*t1 - t2 / t1 - 2)) / ((t1 + t2) * (t1 + t2)) / (t1 * t1) / t2 +
			(l - 1 / (t1 + t2) - 1 / t1 - t1 / t2 / (t1 + t2)*exp(-l*(t1 + t2)) + (t1 + t2) / t1 / t2*exp(-l*t1)) / (t1 + t2) / t1 / d;

		bias /= t1;
		break;
	case ModelBD:
	{
		bias = a*((2 * h < l)*(t1*(exp(-2 * h*(t1 + t2)) - exp(-l*(t1 + t2))) + t2*(t1 + t2)*(l - 2 * h)*exp(-l*(t1 + t2))) / ((t1 + t2) * (t1 + t2)) +
			(exp(-l*(t1 + t2))*(l*t2 *t2 + l*t2*t1 - t1) + t1) / ((t1 + t2) *(t1 + t2)) / d) +
			(1 - a)*((2 * h < l)*(exp(-2 * h*(t1 + t2)) - exp(-l*(t1 + t2))) / (t1 + t2) +
			(1 - exp(-l*(t1 + t2))) / (t1 + t2) / d);
		break;
	}
	case findBestModel:
		break;
	}
	bias = bias * l[0] / bias[0];
	bias /= l;
}
