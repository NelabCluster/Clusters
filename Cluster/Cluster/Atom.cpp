#include "StdAfx.h"
#include "Atom.h"


void Atom::SetPos(AtomPos pos) { _pos = pos; }
AtomPos Atom::GetPos() { return _pos; }

void Atom::SetNote(int note) { _note = note; }
int Atom::GetNote() { return _note; }

AtomPara AtomPara::AtomParaWithType(ATOM_TYPE type)
{
	switch (type)
	{
		case ATOM_Cu : return AtomPara(2.556);
		case ATOM_Au : return AtomPara(2.884);
		break;
	}

	return 2.6;
}

double AtomPara::NearestNeighborSeparation(ATOM_TYPE type)
{
	return AtomParaWithType(type).r0;

	AtomPara para = AtomPara::AtomParaWithType(type);
	switch (para.stru)
	{
	case FCC:
		return para.a / sqrt(2.0);
	case BCC:
		return para.a * sqrt(3.0) / 2;
	case HCP:
		return (para.a > para.c)?para.a:para.c;
	break;
	}
}

double AtomPara::NearestNeighborSeparation( Alloy alloy )
{
	int alloyNum = alloy.size();
	double aveR0 = 0;
	
	for ( int i = 0; i < alloyNum; i++)
		aveR0 += NearestNeighborSeparation( alloy.at(i) );

	aveR0 /= alloyNum;

	return aveR0;
}

string AtomTypeName( ATOM_TYPE type )
{
	switch ( type )
	{
		case ATOM_None: return "LJ";
		case ATOM_Fe: return "Fe";
		case ATOM_Cu: return "Cu";
		case ATOM_Pt: return "Pt";
		case ATOM_Pd: return "Pd";
		case ATOM_Au: return "Au";
		case ATOM_Ag: return "Ag";
		case ATOM_Co: return "Co";
		case ATOM_Zn: return "Zn";
		case ATOM_Cd: return "Cd";
	}
}