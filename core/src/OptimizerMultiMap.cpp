#include "OptimizerMultiMap.hpp"
#include "Map.hpp"

#include <cassert>
#include <cmath>
#include <algorithm>

optimizercore::OptimizerMultiMap::OptimizerMultiMap()
{
	mIsInitialized = false;
}

optimizercore::OptimizerMultiMap::OptimizerMultiMap(MultimapType mapType, int n, int m, int l)
{
	assert(n > 1);
	assert(m > 1);
	assert(l > 0);

	mDimension = n;
	mTightness = m;
	mNumberOfMaps = l;
	mMapType = mapType;
	mRotationPlanes = nullptr;

	p2 = new double[mDimension];
	mTmpVector = SharedVector(p2, utils::array_deleter<double>());

	switch (mMapType)
	{
	case MultimapType::Set:
		mCurrentMap = &OptimizerMultiMap::GetSetMapImage;
		mCurrentInvMap = &OptimizerMultiMap::GetAllSetMapPreimages;
		break;
	case MultimapType::Rotated:
		assert(l <= n*(n - 1) + 1);
		InitRotatedMap();
		mCurrentMap = &OptimizerMultiMap::GetRotatedMapImage;
		mCurrentInvMap = &OptimizerMultiMap::GetAllRotatedMapPreimages;
	}

	mIsInitialized = true;
}

optimizercore::OptimizerMultiMap::~OptimizerMultiMap()
{
//	if (mMapType == MultimapType::Rotated && mIsInitialized)
//		utils::DeleteMatrix(mRotationPlanes, mRotatedMapPlanesCount);
}

int optimizercore::OptimizerMultiMap::GetNumberOfMaps() const
{
	return mNumberOfMaps;
}

void optimizercore::OptimizerMultiMap::GetRotatedMapImage(double x, double y[])
{
	int intx = (int)x;//Номер интервала
	x = x - intx;//дробная часть x
	mapd(x, mTightness, y, mDimension);//получаем точку y[] в исходных координатах
	if (intx == 0 || mNumberOfMaps == 1)
		return;//Если начальный интервал или одна развертка - далее ничего не делаем
	
	int PlaneIndex = intx - 1;//Теперь PlaneNumber - номер перестановки
	PlaneIndex = PlaneIndex % mRotatedMapPlanesCount;
	//Преобразование координат
	double tmp = y[mRotationPlanes[PlaneIndex][1]];
	y[mRotationPlanes[PlaneIndex][1]] = y[mRotationPlanes[PlaneIndex][0]];
	y[mRotationPlanes[PlaneIndex][0]] = -tmp;

	if (intx > mRotatedMapPlanesCount)//Меняем знак преобразования
	{
		y[mRotationPlanes[PlaneIndex][0]] = -y[mRotationPlanes[PlaneIndex][0]];
		y[mRotationPlanes[PlaneIndex][1]] = -y[mRotationPlanes[PlaneIndex][1]];
	}
}
void optimizercore::OptimizerMultiMap::GetSetMapImage(double x, double y[])
{
	double del;
	int i, intx = (int)x;
	x = x - intx;

	if (intx == 0)
		del = 0.0;
	else
		for (i = 1, del = 1; i < intx + 1; del /= 2, i++);

	mapd(x, mTightness + 1, y, mDimension);

	for (i = 0; i < mDimension; i++)
		y[i] = 2 * y[i] + 0.5 - del;
}

void optimizercore::OptimizerMultiMap::GetImage(double x, double y[])
{
	(this->*mCurrentMap)(x, y);
}

int optimizercore::OptimizerMultiMap::GetAllSetMapPreimages(double * p, double xp[])
{
	double xx;
	double del = 0.5;
	for (int i = 1; i < mNumberOfMaps; i++)
	{
		for (int j = 0; j < mDimension; j++)
			p2[j] = (p[j] + del - 0.5) * 0.5;

		xyd(&xx, mTightness + 1, p2, mDimension);
		xp[i] = xx + i;
		del *= 0.5;
	}
	del = 0.0;
	for (int j = 0; j < mDimension; j++)
		p2[j] = (p[j] + del - 0.5) * 0.5;
	xyd(&xx, mTightness + 1, p2, mDimension);
	xp[0] = xx;

	return mNumberOfMaps;
}

int optimizercore::OptimizerMultiMap::GetAllRotatedMapPreimages(double * p, double xp[])
{
	double xx;
	std::copy_n(p, mDimension, p2);
	xyd(&xx, mTightness, p2, mDimension);
	xp[0] = xx;
	//Если одна развертка - далее ничего не делаем
	if (mNumberOfMaps == 1)return 1;

	for (int i = 1; i < mNumberOfMaps; i++)
	{
		std::copy_n(p, mDimension, p2);
		//Обратное преобразование координат
		int PlaneIndex = (i - 1) % mRotatedMapPlanesCount;

		double tmp = p[mRotationPlanes[PlaneIndex][1]];
		p2[mRotationPlanes[PlaneIndex][1]] = -p[mRotationPlanes[PlaneIndex][0]];
		p2[mRotationPlanes[PlaneIndex][0]] = tmp;

		if (i > mRotatedMapPlanesCount)//Меняем знак преобразования
		{
			p2[mRotationPlanes[PlaneIndex][0]] = -p2[mRotationPlanes[PlaneIndex][0]];
			p2[mRotationPlanes[PlaneIndex][1]] = -p2[mRotationPlanes[PlaneIndex][1]];
		}

		xyd(&xx, mTightness, p2, mDimension);
		xp[i] = xx + (i);
	}

	return mNumberOfMaps;
}

int optimizercore::OptimizerMultiMap::GetAllPreimages(double * p, double xp[])
{
	return (this->*mCurrentInvMap)(p, xp);
}

void optimizercore::OptimizerMultiMap::InitRotatedMap()
{
	mRotatedMapPlanesCount = mDimension*(mDimension - 1) / 2;
	mRotationPlanes = utils::AllocateMatrix<int>(mRotatedMapPlanesCount, 2);//Номера осей плоскостей, вокруг которых будут совершаться повороты
	mRotationPlanesHolder = std::shared_ptr<MatrixMemHolder>(
		new MatrixMemHolder(mRotationPlanes, mRotatedMapPlanesCount));

	const int k = 2;//Подмножества из двух элементов
	int plane[k];//Два номера под элементы

	for (int i = 0; i < k; i++)
		plane[i] = i;

	if (mDimension <= k)
	{
		for (int i = 0; i < k; i++)
			mRotationPlanes[0][i] = plane[i];
		return;
	}
	int p = k - 1;
	int counter = 0;//счетчик числа перестановок
	while (p >= 0)
	{
		for (int i = 0; i < k; i++)
			mRotationPlanes[counter][i] = plane[i];
		counter++;

		if (plane[k - 1] == mDimension - 1)
			p--;
		else
			p = k - 1;

		if (p >= 0)
			for (int i = k - 1; i >= p; i--)
				plane[i] = plane[p] + i - p + 1;
	}
}
