#pragma once
#include "Atom.h"

class Clusters
{
public:
	//无团簇类型构造函数
	Clusters(int N);
	//单晶团簇构造函数
	Clusters(int N,ATOM_TYPE type);
	//二合金团簇构造函数
	Clusters(int N,ATOM_TYPE type1,ATOM_TYPE type2,int numberOftype1);
	//多合金构造函数
	Clusters(int N, Alloy &types, AlloyNum &numbers);

	virtual ~Clusters(void);

	//通用
	int GetAtomsNumber();
	double* GetDistancePointer();

	//原子类别
	Alloy GetAlloy();
	AlloyNum GetAlloyNum();
	
	//原子
	Atom GetAtomAtIndex(int index);
	void SetPosOfAtomIndex(int index,AtomPos pos);
	void SetNoteOfAtomIndex(int index,int note);
	void SetAtomOfAtomIndex(int index, Atom atom);

	//能量
	void SetEnergy(double energy);
	double GetEnergy();
	double* GetEnergyPointer();
	void SetEnergyOfAtomIndex(int index, double energy);
	double GetEnergyOfAtomIndex(int index);
	void SetEnergyVectorOfAtoms(vector<double>& energyVector);
	
	//受力
	void SetForceOfCluster(double force);
	double GetForceOfCluster();
	void SetForceXYZOfAtomIndex(int index,double forceX,double forceY,double forceZ);
	double GetForceXOfAtomIndex(int index);
	double GetForceYOfAtomIndex(int index);
	double GetForceZOfAtomIndex(int index);
	void SetForceXYZVectorOfAtoms(vector<double>& forceXVector,vector<double>& forceYVector, vector<double>& forceZVector);
	double* GetForceXPointer();
	double* GetForceYPointer();
	double* GetForceZPointer();

	//序列改变
	void AtomsOrderByZ();

	//位置改变，包括位移，旋转
	void MoveCenterToOrigin();
	void MoveCenter( double dx, double dy, double dz );
	void RotateOriginVector(double theta, double phi, double alpha);
	//void RotateWithVectorAndAngle(int index, double ax, double ay, double az, double angle );
	void RotateWithVectorAndAngle(double ax,double ay, double az, double angle);
	void RotateWithVectorAndAngle(vector<bool>& flag, double ax, double ay, double az, double angle );

private:
	int _N;									//原子总数
	vector<Atom> _atoms;					//原子实例

	Alloy _alloy;							//原子类型
	AlloyNum _alloyNum;						//各类型个数

	vector<double> _R;						//距离矩阵
	double _E;								//势能值
	vector<double> _atomsEnergy;			//每个原子的能量
	double _F;								//团簇力的值
	vector<double> _FX,_FY,_FZ;				//各个原子X,Y,Z方向的受力值

	bool _isChanged;

	void HandlesInit();
};

bool SortClusterByEnergy( Clusters& c1, Clusters& c2);

