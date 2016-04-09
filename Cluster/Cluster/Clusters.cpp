#include "StdAfx.h"
#include "Clusters.h"

#define AtomIndexCheck(index) assert(index >= 0 && index < _N)
#define AtomNumberCheck(Number) assert(Number == _N)

#pragma region 构造函数

Clusters::Clusters(int N)
	:_N(N) { this->HandlesInit(); }

Clusters::Clusters(int N,ATOM_TYPE type)
	:_N(N)
{
	this->HandlesInit();
	_alloy.push_back(type);
	_alloyNum.push_back(N);
}

Clusters::Clusters(int N,ATOM_TYPE type1,ATOM_TYPE type2,int numberOfType1)
	:_N(N)
{
	this->HandlesInit();
	_alloy.push_back( type1 );
	_alloyNum.push_back( numberOfType1 );
	_alloy.push_back(type2);
	_alloyNum.push_back(N-numberOfType1);
}

Clusters::Clusters(int N, Alloy &types, AlloyNum &numbers)
	:_N(N)
{
	this->HandlesInit();
	_alloy = types;
	_alloyNum = numbers;
}

Clusters::~Clusters(void) { }

#pragma endregion 构造函数

#pragma region 私有方法

void Clusters::HandlesInit()
{
	_isChanged = true;
	for (int i = 0; i < _N; i ++)
		_atoms.push_back(Atom(ATOM_None,AtomPos(0,0,0)));
}


#pragma endregion 私有方法

#pragma region 通用

int Clusters::GetAtomsNumber() { return _N; }

double* Clusters::GetDistancePointer()
{
	if (_R.size() != _N * _N)	_R.resize(_N * _N,0);

	if ( _isChanged )
	{
		_isChanged = false;

		for (int i = 0; i < _N; i++ )
			_R[ i*_N + i ] = 0;

		for (int i = 0; i < _N-1; i++)
			for (int j = i+1; j < _N; j++)
			{
				AtomPos posAti = this->_atoms[i].GetPos();
				AtomPos posAtj = this->_atoms[j].GetPos();

				_R[i*_N + j] = sqrt((posAti.x-posAtj.x) * (posAti.x-posAtj.x)  + (posAti.y-posAtj.y) * (posAti.y-posAtj.y) + (posAti.z-posAtj.z) * (posAti.z-posAtj.z));
				if ( fabs(_R[i*_N + j]) < 0.0001 )
				{
					double dr = 0.01;
					posAtj.x += ( RANDI - 0.5 ) * dr;
					posAtj.y += ( RANDI - 0.5 ) * dr;
					posAtj.z += ( RANDI - 0.5 ) * dr;
					this->_atoms[j].SetPos( posAtj );
					_R[i*_N + j] = sqrt((posAti.x-posAtj.x) * (posAti.x-posAtj.x)  + (posAti.y-posAtj.y) * (posAti.y-posAtj.y) + (posAti.z-posAtj.z) * (posAti.z-posAtj.z));
				}
				_R[j*_N + i] = _R[i*_N + j];
			}
	}

	return &_R[0];
}

#pragma endregion 通用

#pragma region 原子类别

Alloy Clusters::GetAlloy() { return _alloy; }

AlloyNum Clusters::GetAlloyNum() { return _alloyNum; }

#pragma endregion 原子类别

#pragma region 原子

Atom Clusters::GetAtomAtIndex(int index)
{
	assert(index >= 0 && index < _N);
	return _atoms[index];
}

void Clusters::SetPosOfAtomIndex(int index, AtomPos pos) { _isChanged = true; _atoms[index].SetPos(pos); }

void Clusters::SetNoteOfAtomIndex(int index, int note) { _isChanged = true; _atoms[index].SetNote(note); }

void Clusters::SetAtomOfAtomIndex(int index, Atom atom) { _isChanged = true; _atoms[index] = atom; }

#pragma endregion 原子

#pragma region 能量

void Clusters::SetEnergy(double energy) { _E = energy; }

double Clusters::GetEnergy() { return _E; }

void Clusters::SetEnergyOfAtomIndex(int index, double energy)
{
	if ( _atomsEnergy.size() != _N ) _atomsEnergy.resize(_N,0);
	
	_atomsEnergy[index] = energy;
}

double Clusters::GetEnergyOfAtomIndex(int index) { return _atomsEnergy[index]; }

void Clusters::SetEnergyVectorOfAtoms(vector<double>& energyVector) { _atomsEnergy = energyVector; }

double* Clusters::GetEnergyPointer() { 
	if ( _atomsEnergy.size() != _N ) _atomsEnergy.resize(_N,0);
	return &_atomsEnergy[0]; 
}

#pragma endregion 能量

#pragma region 受力

void Clusters::SetForceOfCluster(double force) { _F = force; }

double Clusters::GetForceOfCluster() { return _F; }

void Clusters::SetForceXYZOfAtomIndex(int index,double forceX,double forceY,double forceZ) 
{
	AtomIndexCheck(index);

	if ( _FX.size() != _N ) _FX.resize(_N,0);
	if ( _FY.size() != _N ) _FY.resize(_N,0);
	if ( _FZ.size() != _N ) _FZ.resize(_N,0);

	_FX[index] = forceX;
	_FY[index] = forceY;
	_FZ[index] = forceZ;
}

double Clusters::GetForceXOfAtomIndex(int index) { AtomIndexCheck(index); return _FX[index]; }

double Clusters::GetForceYOfAtomIndex(int index) { AtomIndexCheck(index); return _FY[index]; }

double Clusters::GetForceZOfAtomIndex(int index) { AtomIndexCheck(index); return _FZ[index]; }

void Clusters::SetForceXYZVectorOfAtoms(vector<double>& forceXVector,vector<double>& forceYVector, vector<double>& forceZVector)
{
	_FX = forceXVector;
	_FY = forceYVector;
	_FZ = forceZVector;
}

double* Clusters::GetForceXPointer() 
{
	if ( _FX.size() != _N ) _FX.resize(_N);

	return &_FX[0];
}

double* Clusters::GetForceYPointer() 
{
	if ( _FY.size() != _N ) _FY.resize(_N);

	return &_FY[0];
}

double* Clusters::GetForceZPointer() 
{
	if ( _FZ.size() != _N ) _FZ.resize(_N);

	return &_FZ[0];
}

#pragma endregion 受力

#pragma region 序列改变

void Clusters::AtomsOrderByZ()
{	
	_isChanged = true;

	for (int i = _N - 1; i > 1; i--)
	{
		for (int j = 0; j < i; j++)
		{
			if (_atoms[j].GetPos().z < _atoms[j+1].GetPos().z)
				swap(_atoms[j],_atoms[j+1]);
		}
	}
}

#pragma endregion 序列改变

#pragma region 原子的位移，旋转，顺序变化

void Clusters::MoveCenterToOrigin()
{
	double middleOfX = 0, middleOfY = 0, middleOfZ = 0;

	for (int i = 0; i < _N; i ++)
	{
		AtomPos pos = _atoms[i].GetPos();
		middleOfX += pos.x;
		middleOfY += pos.y;
		middleOfZ += pos.z;
	}

	middleOfX /= _N; middleOfY /= _N; middleOfZ /= _N;

	this->MoveCenter(-middleOfX,-middleOfY,-middleOfZ);
}

void Clusters::MoveCenter( double dx, double dy, double dz )
{
	//_isChanged = true;

	for (int i= 0; i < _N; i++)
	{
		AtomPos pos = _atoms[i].GetPos();
		pos.x += dx;
		pos.y += dy;
		pos.z += dz;
		_atoms[i].SetPos(pos);
	}
}

//http://blog.csdn.net/changbaolong/article/details/8307052
void Clusters::RotateOriginVector(double theta, double phi, double alpha)
{
	double ax,ay,az;
	ax = sin(theta) * cos(phi);
	ay = sin(theta) * sin(phi);
	az = cos(theta);
	this->RotateWithVectorAndAngle(ax,ay,az,alpha);
}

//void Clusters::RotateWithVectorAndAngle(int index, double ax, double ay, double az, double angle )
//{
//	vector<bool> flag(_N,false);
//	flag[index] = true;
//	this->RotateWithVectorAndAngle(flag,ax,ay,az,angle);
//}

void Clusters::RotateWithVectorAndAngle(double ax,double ay, double az, double angle)
{
	vector<bool> flag(_N,true);
	this->RotateWithVectorAndAngle(flag,ax,ay,az,angle);
}

void Clusters::RotateWithVectorAndAngle(vector<bool>& flag, double ax, double ay, double az, double angle )
{
	_isChanged = true;

	double A[3][3],B[3][3],M[3][3],I[3][3];

	A[0][0] = ax * ax; A[0][1] = ax * ay; A[0][2] = ax * az;
	A[1][0] = ay * ax; A[1][1] = ay * ay; A[1][2] = ay * az;
	A[2][0] = az * ax; A[2][1] = az * ay; A[2][2] = az * az;

	B[0][0] =   0; B[0][1] = -az; B[0][2] =  ay;
	B[1][0] =  az; B[1][1] =   0; B[1][2] = -ax;
	B[2][0] = -ay; B[2][1] =  ax; B[2][2] =   0;

	I[0][0] = 1; I[0][1] = 0; I[0][2] = 0;
	I[1][0] = 0; I[1][1] = 1; I[1][2] = 0;
	I[2][0] = 0; I[2][1] = 0; I[2][2] = 1;

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			M[i][j] = A[i][j] + cos(angle) * (I[i][j] - A[i][j]) + sin(angle) * B[i][j];
		}
	}

	for (int i = 0; i < _N; i++)
	{
		if ( flag[i] )
		{
			double coodAfterRotate[3];
			for (int j = 0; j < 3; j++)
			{
				AtomPos pos = _atoms[i].GetPos();
				coodAfterRotate[j] = pos.x * M[j][0] + pos.y * M[j][1] + pos.z * M[j][2];
			}
			_atoms[i].SetPos(AtomPos(coodAfterRotate[0],coodAfterRotate[1],coodAfterRotate[2]));
		}
	}
}

#pragma endregion 原子的位移，旋转，顺序变化

bool SortClusterByEnergy( Clusters& c1, Clusters& c2)
{
	return c1.GetEnergy() > c2.GetEnergy();
}






