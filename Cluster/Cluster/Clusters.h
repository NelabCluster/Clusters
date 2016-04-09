#pragma once
#include "Atom.h"

class Clusters
{
public:
	//���Ŵ����͹��캯��
	Clusters(int N);
	//�����Ŵع��캯��
	Clusters(int N,ATOM_TYPE type);
	//���Ͻ��Ŵع��캯��
	Clusters(int N,ATOM_TYPE type1,ATOM_TYPE type2,int numberOftype1);
	//��Ͻ��캯��
	Clusters(int N, Alloy &types, AlloyNum &numbers);

	virtual ~Clusters(void);

	//ͨ��
	int GetAtomsNumber();
	double* GetDistancePointer();

	//ԭ�����
	Alloy GetAlloy();
	AlloyNum GetAlloyNum();
	
	//ԭ��
	Atom GetAtomAtIndex(int index);
	void SetPosOfAtomIndex(int index,AtomPos pos);
	void SetNoteOfAtomIndex(int index,int note);
	void SetAtomOfAtomIndex(int index, Atom atom);

	//����
	void SetEnergy(double energy);
	double GetEnergy();
	double* GetEnergyPointer();
	void SetEnergyOfAtomIndex(int index, double energy);
	double GetEnergyOfAtomIndex(int index);
	void SetEnergyVectorOfAtoms(vector<double>& energyVector);
	
	//����
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

	//���иı�
	void AtomsOrderByZ();

	//λ�øı䣬����λ�ƣ���ת
	void MoveCenterToOrigin();
	void MoveCenter( double dx, double dy, double dz );
	void RotateOriginVector(double theta, double phi, double alpha);
	//void RotateWithVectorAndAngle(int index, double ax, double ay, double az, double angle );
	void RotateWithVectorAndAngle(double ax,double ay, double az, double angle);
	void RotateWithVectorAndAngle(vector<bool>& flag, double ax, double ay, double az, double angle );

private:
	int _N;									//ԭ������
	vector<Atom> _atoms;					//ԭ��ʵ��

	Alloy _alloy;							//ԭ������
	AlloyNum _alloyNum;						//�����͸���

	vector<double> _R;						//�������
	double _E;								//����ֵ
	vector<double> _atomsEnergy;			//ÿ��ԭ�ӵ�����
	double _F;								//�Ŵ�����ֵ
	vector<double> _FX,_FY,_FZ;				//����ԭ��X,Y,Z���������ֵ

	bool _isChanged;

	void HandlesInit();
};

bool SortClusterByEnergy( Clusters& c1, Clusters& c2);

