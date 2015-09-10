#ifndef OPTIMIZER_DATA_STRUCTURES_HPP
#define OPTIMIZER_DATA_STRUCTURES_HPP

#include "Map.hpp"

namespace optimizercore	{

	enum class StopCriterionType : int { OptimalPoint = 0, Precision = 1 };
	enum class LocalTuningMode : int { None = 0, Maximum = 1, Adaptive = 2 };

	struct OptimizerTrialPoint
	{
		double x, val;
		int v;

		OptimizerTrialPoint()
		{}
		OptimizerTrialPoint(double _x, double _val, int _index)
		{
			x = _x;
			val = _val;
			v = _index;
		}
	};
	
	inline bool operator<(const OptimizerTrialPoint& t1,
		const OptimizerTrialPoint& t2){
		return t1.x < t2.x;
	}
	

	struct OptimaizerInterval
	{
		OptimizerTrialPoint left, right;
		double rank, localM;

		OptimaizerInterval()
		{}

		OptimaizerInterval(OptimizerTrialPoint _left,
			OptimizerTrialPoint _right, double _rank, double _localM)
		{
			left = _left;
			right = _right;
			rank = _rank;
			localM = _localM;
		}
	};

	struct OptimizerParameters
	{
		int maxIterationsNumber;
		int numberOfThreads;
		double eps;
		double* r;
		double* reserves;
		unsigned algDimention;

		int localExponent;
		int localMixParameter;
		int localAlgStartIterationNumber;
		bool localVerification;

		LocalTuningMode localTuningMode;
		MapType mapType;
		int mapTightness;
		OptimizerParameters()
		{
			r = reserves = nullptr;
		}
		OptimizerParameters(int _maxIterationsNumber, int _numberOfThreads, double _eps,
			double* _r, double* _reserves, unsigned _algDimention, int _localExponent,
			int _localMixParameter, int _localAlgStartIterationNumber, MapType _mapType, int _mapTightness,
			bool _localVerification)
		{
			maxIterationsNumber = _maxIterationsNumber;
			numberOfThreads = _numberOfThreads;
			eps = _eps;
			r = _r;
			reserves = _reserves;
			algDimention = _algDimention;
			localExponent = _localExponent;
			localMixParameter = _localMixParameter;
			localAlgStartIterationNumber = _localAlgStartIterationNumber;
			mapType = _mapType;
			mapTightness = _mapTightness;
			localVerification = _localVerification;
		}
		OptimizerParameters(int _maxIterationsNumber, int _numberOfThreads, double _eps,
			double* _r, double* _reserves, unsigned _algDimention, int _localExponent,
			int _localMixParameter, int _localAlgStartIterationNumber, MapType _mapType, int _mapTightness,
			bool _localVerification, LocalTuningMode _localTuningMode)
		{
			maxIterationsNumber = _maxIterationsNumber;
			numberOfThreads = _numberOfThreads;
			eps = _eps;
			r = _r;
			reserves = _reserves;
			algDimention = _algDimention;
			localExponent = _localExponent;
			localMixParameter = _localMixParameter;
			localAlgStartIterationNumber = _localAlgStartIterationNumber;
			mapType = _mapType;
			mapTightness = _mapTightness;
			localVerification = _localVerification;
			localTuningMode = _localTuningMode;
		}
	};

	class IndxSet
	{

	private:
		int mCurrentSetSize;
		int mMemReallocationStep;
		int mMemSize;
		OptimizerTrialPoint* mMem;
		double mCurrentZMin;

	public:

		IndxSet();
		IndxSet(int memSize, int memReallocStep);

		void Add(const OptimizerTrialPoint& trial);
		OptimizerTrialPoint Get(int i) const;

		int GetSize() const;
		double GetMinimumValue() const;
		OptimizerTrialPoint GetMinimumPoint() const;

		void Reset();
		OptimizerTrialPoint& operator[](int i);

		~IndxSet();

	};

}

#endif