#pragma once;

enum ATOM_TYPE
{
	ATOM_None,
	ATOM_Fe = 1,
	ATOM_Cu,
	ATOM_Pt,
	ATOM_Pd,
	ATOM_Au,
	ATOM_Ag,
	ATOM_Co,
	ATOM_Zn,
	ATOM_Cd
};

typedef vector<ATOM_TYPE> Alloy;
typedef vector<int> AlloyNum;

struct AtomPos 
{
	double x;
	double y;
	double z;
	AtomPos(double a,double b,double c):x(a),y(b),z(c){}
};

class Atom
{
public:

	Atom(int note,AtomPos pos):_note(note),_pos(pos){};
	
	void SetPos(AtomPos pos);
	AtomPos GetPos();

	void SetNote(int note);
	int GetNote();

protected:
	int _note;
	AtomPos _pos;
};

enum CRYSTAL_STRUCTURE
{
	BCC,FCC,HCP
};

struct AtomPara
{
	CRYSTAL_STRUCTURE stru;
	double a;
	double c;
	double r0;

	AtomPara(double m):r0(m) {};
	AtomPara(CRYSTAL_STRUCTURE s, double m, double n):stru(s),a(m),c(n) {};
	static AtomPara AtomParaWithType( ATOM_TYPE type );
	static double NearestNeighborSeparation( ATOM_TYPE type );
	static double NearestNeighborSeparation( Alloy alloy );
};

string AtomTypeName( ATOM_TYPE type ); 