// added by Jieling, 05/02/2018
// to do delaunay triangulation

#ifndef DELAUNAY_CPP
#define DELAUNAY_CPP

#include <iostream>
#include <vector>
#include <cmath>
#include <memory.h>
#include <algorithm>
#include "delaunay.h"

using namespace std;

// type conversion implemented for INTEL architecture
// long = 4 bytes
// float = 4 bytes
// double = 8 bytes

void CExact::Dump()
{
	int i;
	printf("%c s:%d e:%d\t", Plus ? '+' : '-', Size, Exp);
	for (i = 0;i<Size;i++) printf(" %02x", S[i]);
	printf("\t%#.18lg\n", (double)(*this));
}

CExact::CExact(const CExact&E)
{
	Plus = E.Plus; Exp = E.Exp; Size = E.Size; S = 0;
	if (Size) { S = new unsigned char[Size]; memcpy(S, E.S, Size); }
}

const CExact& CExact::operator=(const CExact&E)
{
	if (this == &E) return *this;
	Plus = E.Plus; Exp = E.Exp; Size = E.Size; delete[] S;
	if (Size)
	{
		S = new unsigned char[Size]; memcpy(S, E.S, Size);
	}
	else S = 0;
	return *this;
}

void CExact::Init(long n)
{
	int i;
	Exp = 0;
	if (n<0) { n = -n; Plus = false; }
	else Plus = true;
	if (!n) { S = 0; Size = 0; return; }
	S = new unsigned char[sizeof(long)];
	for (i = 0;i<sizeof(long);i++)
	{
		if (n & 0xFF) break;
		Exp++; n >>= 8;
	}
	*((long*)S) = n;
	for (Size = sizeof(long);Size;Size--) if (S[Size - 1]) break;
}

void CExact::Init(double d)
{
	int i, j, k, m;
	unsigned char *a = (unsigned char *)&d;
	unsigned char b = 1;

	if (!d) { S = 0; Size = 0; Exp = 0; Plus = true; return; }
	k = ((a[7] & 0x7f) << 4) + (a[6] >> 4);
	if (k >= 2047)	//// Infinity - return 0
	{
		S = 0; Size = 0; Exp = 0; Plus = true; return;
	}

	S = new unsigned char[8]; memset(S, 0, 8);
	Plus = !(a[7] & 0x80);
	if (!k) { b = 0; k = 1; } ////	SUBNORMAL!
	Exp = (k + 1) / 8 - 128 - 7;

	m = (k + 1) % 8; S[7] = b << m;
	if (m <= 4)
	{
		m = 4 - m; S[7] += (a[6] & 0xf) >> m; j = 6;
	}
	else
	{
		m -= 4; S[7] += (a[6] & 0xf) << m;
		m = 8 - m; S[7] += a[5] >> m; j = 5;
	}

	for (i = 6;i >= 0;i--)
	{
		S[i] = a[j] << (8 - m); j--;
		if (j<0) break;
		S[i] += a[j] >> m;
	}
	for (Size = 8;Size;Size--) if (S[Size - 1]) break;
	for (j = 0;j<Size;j++) if (S[j]) break;

	if (j)
	{
		Size -= j; Exp += j;
		for (i = 0;i<Size;i++) S[i] = S[i + j];
	}
}

CExact::operator double() const
{
	int i, j, k, m, sn;
	unsigned char b, c;
	double d;
	unsigned char *a = (unsigned char *)&d;
	if (!Size) return 0;
	memset(a, 0, sizeof(double));
	k = Size + Exp - 1;
	b = 0x80; c = S[Size - 1];
	for (m = 7;m;m--) { if (b&c) break; b >>= 1; }
	// need check for too small or too big k
	k = k * 8 + m + 1023;
	if (k >= 2047)	// Infinity
	{
		a[7] = Plus ? 0x7f : 0xff; a[6] = 0xf0; return d;
	}
	if (k<-111)	// Zero
	{
		a[7] = Plus ? 0 : 0x80; return d;
	}

	if (k<1)		// Subnormal
	{
		sn = k - 1; k = 1;
	}
	else sn = 0;

	a[7] = k >> 4; a[6] = (k & 0xf) << 4;
	if (!Plus) a[sizeof(double) - 1] |= 0x80;
	j = Size - 1;
	if (m >= 4) { m -= 4; a[6] |= (S[j] >> m) & 0xf; }
	else
	{
		a[6] |= (S[j] << (4 - m)) & 0xf; m += 4; j--;
		if (j >= 0) a[6] |= S[j] >> m;
	}
	if (j >= 0) for (i = 5;i >= 0;i--)
	{
		a[i] = S[j] << (8 - m); j--;
		if (j<0) break;
		a[i] |= S[j] >> m;
	}
	if (sn) d *= pow(2.0, sn);
	return d;
}

char CExact::CompAbs(const CExact&E) const
{
	int a, b;
	a = Size + Exp; b = E.Size + E.Exp;
	if (a>b) return 1;
	if (a<b) return -1;
	a = Size; b = E.Size;
	while (a&&b)
	{
		a--; b--;
		if (S[a]>E.S[b]) return 1;
		if (S[a]<E.S[b]) return -1;
	}
	if (a) return 1;
	if (b) return -1;
	return 0;
}

char CExact::Compare(const CExact&E) const
{
	if (Plus)
	{
		if (E.Plus) return CompAbs(E);
		else return 1;
	}
	else
	{
		if (E.Plus) return -1;
		else return -CompAbs(E);
	}
}

bool CExact::operator==(const CExact&E) const { return !Compare(E); }
bool CExact::operator!=(const CExact&E) const { return Compare(E) != 0; }
bool CExact::operator>(const CExact&E) const { return Compare(E) == 1; }
bool CExact::operator<(const CExact&E) const { return Compare(E) == -1; }
bool CExact::operator>=(const CExact&E) const { return Compare(E) != -1; }
bool CExact::operator<=(const CExact&E) const { return Compare(E) != 1; }

void CExact::Add(const CExact&A, const CExact&B, int&Sz, int&Ex, unsignedcharptr&St) const
{
	int i, j, k, n;
	int a, b;

	Ex = A.Exp<B.Exp ? A.Exp : B.Exp;
	a = A.Exp + A.Size; b = B.Exp + B.Size;
	n = a>b ? a : b; n++; Sz = n - Ex;
	St = new unsigned char[Sz];

	b = 0; j = Ex - A.Exp; k = Ex - B.Exp;
	for (i = 0;i<Sz;i++)
	{
		a = b;
		if (j >= 0 && j<A.Size) a += A.S[j];
		if (k >= 0 && k<B.Size) a += B.S[k];
		St[i] = a & 0xff;
		j++; k++; b = a >> 8;
	}
	for (;Sz;Sz--) if (St[Sz - 1]) break;
	for (j = 0;j<Sz;j++) if (St[j]) break;
	if (j)
	{
		Sz -= j; Ex += j;
		for (i = 0;i<Sz;i++) St[i] = St[i + j];
	}
}

//// Assumes A>B>=0
void CExact::Sub(const CExact&A, const CExact&B, int&Sz, int&Ex, unsignedcharptr&St) const
{
	int i, j, k, n;
	int a, b;

	Ex = A.Exp<B.Exp ? A.Exp : B.Exp;
	n = A.Exp + A.Size; //b=B.Exp+B.Size;
						//n=a>b?a:b; 
	Sz = n - Ex;
	St = new unsigned char[Sz];

	b = 0; j = Ex - A.Exp; k = Ex - B.Exp;
	for (i = 0;i<Sz;i++)
	{
		a = 0x100 - b;
		if (j >= 0 && j<A.Size) a += A.S[j];
		if (k >= 0 && k<B.Size) a -= B.S[k];
		St[i] = a & 0xff;
		j++; k++; b = a<0x100 ? 1 : 0;
	}
	for (;Sz;Sz--) if (St[Sz - 1]) break;
	for (j = 0;j<Sz;j++) if (St[j]) break;
	if (j)
	{
		Sz -= j; Ex += j;
		for (i = 0;i<Sz;i++) St[i] = St[i + j];
	}
}

const CExact CExact::operator-() const
{
	CExact A = (*this); if (Size) A.Plus = !Plus; return A;
}

const CExact CExact::operator+(const CExact&E) const
{
	int Sz, Ex;
	unsigned char *St;
	CExact A;
	if (Plus == E.Plus)
	{
		Add(*this, E, Sz, Ex, St);
		A.Size = Sz; A.Exp = Ex; A.S = St; A.Plus = Plus;
	}
	else
	{
		switch (CompAbs(E))
		{
		case 0:  return A;
		case 1:  Sub(*this, E, Sz, Ex, St); A.Plus = Plus; break;
		default: Sub(E, *this, Sz, Ex, St); A.Plus = E.Plus;
		}
		A.Size = Sz; A.Exp = Ex; A.S = St;
	}
	return A;
}

const CExact& CExact::operator+=(const CExact&E)
{
	int Sz, Ex;
	unsigned char *St;
	if (Plus == E.Plus)
	{
		Add(*this, E, Sz, Ex, St);
		delete[] S; Size = Sz; Exp = Ex; S = St;
	}
	else
	{
		switch (CompAbs(E))
		{
		case 0:  delete[] S; Size = 0; Exp = 0; Plus = true; return *this;
		case 1:  Sub(*this, E, Sz, Ex, St); break;
		default: Sub(E, *this, Sz, Ex, St); Plus = !Plus;
		}
		delete[]S; Size = Sz; Exp = Ex; S = St;
	}
	return *this;
}

const CExact CExact::operator-(const CExact&E) const
{
	int Sz, Ex;
	unsigned char *St;
	CExact A;
	if (Plus != E.Plus)
	{
		Add(*this, E, Sz, Ex, St);
		A.Size = Sz; A.Exp = Ex; A.S = St; A.Plus = Plus;
	}
	else
	{
		switch (CompAbs(E))
		{
		case 0:  return A;
		case 1:  Sub(*this, E, Sz, Ex, St); A.Plus = Plus; break;
		default: Sub(E, *this, Sz, Ex, St); A.Plus = !Plus;
		}
		A.Size = Sz; A.Exp = Ex; A.S = St;
	}
	return A;
}

const CExact& CExact::operator-=(const CExact&E)
{
	int Sz, Ex;
	unsigned char *St;
	if (Plus != E.Plus)
	{
		Add(*this, E, Sz, Ex, St);
		delete[] S; Size = Sz; Exp = Ex; S = St;
	}
	else
	{
		switch (CompAbs(E))
		{
		case 0:  delete[] S; Size = 0; Exp = 0; Plus = true; return *this;
		case 1:  Sub(*this, E, Sz, Ex, St); break;
		default: Sub(E, *this, Sz, Ex, St); Plus = !Plus;
		}
		delete[]S; Size = Sz; Exp = Ex; S = St;
	}
	return *this;
}

//// Assumes A,B!=0
void CExact::Mul(const CExact&A, const CExact&B, int&Sz, int&Ex, unsignedcharptr&St) const
{
	int i, j;
	int a, b;

	Ex = A.Exp + B.Exp;
	Sz = A.Size + B.Size;
	St = new unsigned char[Sz];
	memset(St, 0, Sz);

	for (i = 0;i<A.Size;i++)
	{
		b = 0;
		for (j = 0;j<B.Size;j++)
		{
			a = A.S[i] * B.S[j] + St[i + j] + b;
			St[i + j] = a & 0xff;
			b = a >> 8;
		}
		St[i + j] = b;
	}

	for (;Sz;Sz--) if (St[Sz - 1]) break;
	for (j = 0;j<Sz;j++) if (St[j]) break;
	if (j)
	{
		Sz -= j; Ex += j;
		for (i = 0;i<Sz;i++) St[i] = St[i + j];
	}
}

const CExact CExact::operator*(const CExact&E) const
{
	int Sz, Ex;
	unsigned char *St;
	CExact A;
	if (!Size || !E.Size) return A;
	Mul(*this, E, Sz, Ex, St);
	A.Size = Sz; A.Exp = Ex; A.S = St;
	A.Plus = E.Plus == Plus;
	return A;
}

const CExact& CExact::operator*=(const CExact&E)
{
	int Sz, Ex;
	unsigned char *St;
	if (!Size) return *this;
	if (!E.Size) { delete[] S; Size = 0; Exp = 0; Plus = true; return *this; }
	Mul(*this, E, Sz, Ex, St);
	Size = Sz; Exp = Ex; Plus = E.Plus == Plus;
	delete[] S; S = St;
	return *this;
}

inline CInexact::operator bool() const { return fabs(D)>E; }
inline bool CInexact::is_zero() const { return fabs(D) <= E; }
inline const CInexact CInexact::operator-() const { CInexact A = *this; A.D = -D; return A; }
inline const CInexact CInexact::operator+(const CInexact&I) const { CInexact A = *this; return (A += I); }
inline const CInexact CInexact::operator-(const CInexact&I) const { CInexact A = *this; return (A -= I); }
inline const CInexact CInexact::operator*(const CInexact&I) const { CInexact A = *this; return (A *= I); }

inline const CInexact& CInexact::operator+=(const CInexact&I)
{
	D += I.D; E += I.E;
	double e = Err(); if (E<e) E = e;
	return *this;
}

inline const CInexact& CInexact::operator-=(const CInexact&I)
{
	D -= I.D; E += I.E;
	double e = Err(); if (E<e) E = e;
	return *this;
}

inline const CInexact& CInexact::operator*=(const CInexact&I)
{
	E = (fabs(D) + E)*I.E + fabs(I.D)*E; D *= I.D;
	double e = Err(); if (E<e) E = e;
	return *this;
}

void CInexact::Dump()
{
	printf("%#.18lg\t+-%#.18lg\n", D, E);
}

double CInexact::Err()
{
	long e;
	unsigned char *a = (unsigned char *)&D;
	e = ((a[7] & 0x7f) << 4) + (a[6] >> 4);
	if (!e) e = 1;
	e -= 1023 + 52 - 1;
	return pow(2.0, e);
}

void CInexact::Init(float d)
{
	long e;
	unsigned char *a = (unsigned char *)&d;
	D = d;
	e = (a[3] & 0x7f) << 1;
	if (a[2] & 0x80) e++;
	if (!e) e = 1;
	e -= 127 + 23 - 1;
	E = pow(2.0, e);
}

void CInexact::Init(double d) { D = d; E = Err(); }

CXYZW::~CXYZW()
{
	/*S0 = NULL;
	p = NULL;*/
}

CEdge::~CEdge()
{
	/*A = NULL;
	B = NULL;
	Cell = NULL;
	Rev = NULL;
	Enx = NULL;
	Fnx = NULL;
	S1 = NULL;
	S2 = NULL;*/
}

CCell::~CCell()
{
	/*VV[0] = NULL;
	VV[1] = NULL;
	VV[2] = NULL;
	ab = NULL;
	bc = NULL;
	ca = NULL;
	CC[0] = NULL;
	CC[1] = NULL;
	CC[2] = NULL;*/
}

CD3DW::CD3DW(double H)
{
	CXYZW *a, *b, *c, *d;
	L = H; 
	V0.N = -1; V1.N = -2; V2.N = -3; V3.N = -4;
	V0.set(0, 0, L, 0); 
	V1.set(0, L,-L, 0); 
	V2.set(-L,-L,-L,0); 
	V3.set(L, -L,-L,0);
	a = &V0; b = &V1; c = &V2; d = &V3;
	
	EdgePtr abc, bca, cab, bac, cba, acb;
	new6(abc, bca, cab, bac, cba, acb, a, b, c);
	EdgePtr abd, bda, dab, bad, dba, adb;
	new6(abd, bda, dab, bad, dba, adb, a, b, d);
	EdgePtr adc, dca, cad, dac, cda, acd;
	new6(adc, dca, cad, dac, cda, acd, a, d, c);
	EdgePtr dbc, bcd, cdb, bdc, cbd, dcb;
	new6(dbc, bcd, cdb, bdc, cbd, dcb, d, b, c);

	abc->Fnx = abd; abd->Fnx = abc;
	bac->Fnx = bad; bad->Fnx = bac;

	bca->Fnx = bcd; bcd->Fnx = bca;
	cba->Fnx = cbd; cbd->Fnx = cba;

	cab->Fnx = cad; cad->Fnx = cab;
	acb->Fnx = acd; acd->Fnx = acb;

	adb->Fnx = adc; adc->Fnx = adb;
	dab->Fnx = dac; dac->Fnx = dab;

	bda->Fnx = bdc; bdc->Fnx = bda;
	dba->Fnx = dbc; dbc->Fnx = dba;

	cda->Fnx = cdb; cdb->Fnx = cda;
	dca->Fnx = dcb; dcb->Fnx = dca;

	CCell*C = new CCell; Cs.push_back(C);
	C->VV[0] = a; C->VV[1] = b; C->VV[2] = c; C->VV[3] = d;
	C->abc = abc; C->adb = adb; C->acd = acd; C->bdc = bdc;

	abc->Cell = C; bca->Cell = C; cab->Cell = C;
	bac->Cell = C; cba->Cell = C; acb->Cell = C;
	abd->Cell = C; bda->Cell = C; dab->Cell = C;
	bad->Cell = C; dba->Cell = C; adb->Cell = C;
	adc->Cell = C; dca->Cell = C; cad->Cell = C;
	dac->Cell = C; cda->Cell = C; acd->Cell = C;
	dbc->Cell = C; bcd->Cell = C; cdb->Cell = C;
	bdc->Cell = C; cbd->Cell = C; dcb->Cell = C;

	abc->Flag = 2; bca->Flag = 2; cab->Flag = 2;
	bac->Flag = 2; cba->Flag = 2; acb->Flag = 2;
	abd->Flag = 2; bda->Flag = 2; dab->Flag = 2;
	bad->Flag = 2; dba->Flag = 2; adb->Flag = 2;
	adc->Flag = 2; dca->Flag = 2; cad->Flag = 2;
	dac->Flag = 2; cda->Flag = 2; acd->Flag = 2;
	dbc->Flag = 2; bcd->Flag = 2; cdb->Flag = 2;
	bdc->Flag = 2; cbd->Flag = 2; dcb->Flag = 2;

	N23 = 0; N32 = 0; N41 = 0; NLoc = 0; Nep = 0; Nei = 0;
	Ntp = 0; Nrp = 0; Ndd = 0; Nti = 0;
}

CD3DW::~CD3DW()
{
	int i;
	for (i=0;i<(int)Vs.size();i++) {delete Vs[i];Vs[i] = NULL;}
	for (i=0;i<(int)Es.size();i++) {delete Es[i];Es[i] = NULL;}
	for (i=0;i<(int)Cs.size();i++) {delete Cs[i];Cs[i] = NULL;}
	Vs.clear();
	Es.clear();
	Cs.clear();
	//delete V0;delete V1;delete V2;delete V3;
	//V0 = NULL;V1 = NULL;V2 = NULL;
}

void CD3DW::add(double x, double y, double z, double w, int n, void *l)
{
	CXYZW *a, *b, *c, *d, *e;
	CCell *abcd, *abce, *acde, *adbe, *bdce;

	e = new CXYZW; e->Label = l;
	e->set(x, y, z, w); 
	e->N = n;
	Vs.push_back(e);
	double tx = (abs(x) + abs(y) + abs(z) + abs(w)) * 1000;
	if (L < tx) fprintf(stderr, "coordinates too large!");

	EdgePtr abc, bca, cab;
	abc = locate(e); bca = abc->Enx; cab = bca->Enx;

	EdgePtr adb, dba, bad;
	bad = abc->Fnx->Rev; adb = bad->Enx; dba = adb->Enx;
	EdgePtr acd, cda, dac;
	acd = cab->Fnx->Rev; cda = acd->Enx; dac = cda->Enx;
	EdgePtr bdc, dcb, cbd;
	cbd = bca->Fnx->Rev; bdc = cbd->Enx; dcb = bdc->Enx;
	a = abc->A; b = abc->B; c = abc->C; d = abc->Fnx->C;

	if ((*e) == (*a))
	{	
		N41++; return;
	}
	if ((*e) == (*b))
	{	
		N41++; return;
	}
	if ((*e) == (*c))
	{	
		N41++; return;
	}
	if ((*e) == (*d))
	{	
		N41++; return;
	}

	if (!insphere(abc, e))
	{	
		N41++; return;
	}

	EdgePtr abe, bea, eab, bae, eba, aeb;
	new6(abe, bea, eab, bae, eba, aeb, a, b, e);
	EdgePtr bce, ceb, ebc, cbe, ecb, bec;
	new6(bce, ceb, ebc, cbe, ecb, bec, b, c, e);
	EdgePtr cae, aec, eca, ace, eac, cea;
	new6(cae, aec, eca, ace, eac, cea, c, a, e);
	EdgePtr ade, dea, ead, dae, eda, aed;
	new6(ade, dea, ead, dae, eda, aed, a, d, e);
	EdgePtr bde, deb, ebd, dbe, edb, bed;
	new6(bde, deb, ebd, dbe, edb, bed, b, d, e);
	EdgePtr cde, dec, ecd, dce, edc, ced;
	new6(cde, dec, ecd, dce, edc, ced, c, d, e);

	abe->Fnx = abc->Fnx; abc->Fnx = abe;
	bae->Fnx = bad->Fnx; bad->Fnx = bae;

	ace->Fnx = acd->Fnx; acd->Fnx = ace;
	cae->Fnx = cab->Fnx; cab->Fnx = cae;

	ade->Fnx = adb->Fnx; adb->Fnx = ade;
	dae->Fnx = dac->Fnx; dac->Fnx = dae;

	bce->Fnx = bca->Fnx; bca->Fnx = bce;
	cbe->Fnx = cbd->Fnx; cbd->Fnx = cbe;

	bde->Fnx = bdc->Fnx; bdc->Fnx = bde;
	dbe->Fnx = dba->Fnx; dba->Fnx = dbe;

	cde->Fnx = cda->Fnx; cda->Fnx = cde;
	dce->Fnx = dcb->Fnx; dcb->Fnx = dce;

	ead->Fnx = eac; eac->Fnx = eab; eab->Fnx = ead;
	aed->Fnx = aeb; aeb->Fnx = aec; aec->Fnx = aed;

	eba->Fnx = ebc; ebc->Fnx = ebd; ebd->Fnx = eba;
	bea->Fnx = bed; bed->Fnx = bec; bec->Fnx = bea;

	eca->Fnx = ecd; ecd->Fnx = ecb; ecb->Fnx = eca;
	cea->Fnx = ceb; ceb->Fnx = ced; ced->Fnx = cea;

	eda->Fnx = edb; edb->Fnx = edc; edc->Fnx = eda;
	dea->Fnx = dec; dec->Fnx = deb; deb->Fnx = dea;

	abcd = abc->Cell;
	abce = new CCell; Cs.push_back(abce);
	abce->abc = abc; abce->adb = aeb; abce->acd = ace; abce->bdc = bec;
	abce->VV[0] = a; abce->VV[1] = b; abce->VV[2] = c; abce->VV[3] = e;

	acde = new CCell; Cs.push_back(acde);
	acde->abc = acd; acde->adb = aec; acde->acd = ade; acde->bdc = ced;
	acde->VV[0] = a; acde->VV[1] = c; acde->VV[2] = d; acde->VV[3] = e;

	adbe = new CCell; Cs.push_back(adbe);
	adbe->abc = adb; adbe->adb = aed; adbe->acd = abe; adbe->bdc = deb;
	adbe->VV[0] = a; adbe->VV[1] = d; adbe->VV[2] = b; adbe->VV[3] = e;

	bdce = new CCell; Cs.push_back(bdce);
	bdce->abc = bdc; bdce->adb = bed; bdce->acd = bce; bdce->bdc = dec;
	bdce->VV[0] = b; bdce->VV[1] = d; bdce->VV[2] = c; bdce->VV[3] = e;

	abcd->CC[0] = abce; abcd->CC[1] = acde; abcd->CC[2] = adbe; abcd->CC[3] = bdce;

	abc->Cell = abce; bca->Cell = abce; cab->Cell = abce;
	ace->Cell = abce; cea->Cell = abce; eac->Cell = abce;
	aeb->Cell = abce; eba->Cell = abce; bae->Cell = abce;
	bec->Cell = abce; ecb->Cell = abce; cbe->Cell = abce;

	acd->Cell = acde; cda->Cell = acde; dac->Cell = acde;
	ade->Cell = acde; dea->Cell = acde; ead->Cell = acde;
	aec->Cell = acde; eca->Cell = acde; cae->Cell = acde;
	dce->Cell = acde; ced->Cell = acde; edc->Cell = acde;

	adb->Cell = adbe; dba->Cell = adbe; bad->Cell = adbe;
	abe->Cell = adbe; bea->Cell = adbe; eab->Cell = adbe;
	aed->Cell = adbe; eda->Cell = adbe; dae->Cell = adbe;
	deb->Cell = adbe; ebd->Cell = adbe; bde->Cell = adbe;

	bdc->Cell = bdce; dcb->Cell = bdce; cbd->Cell = bdce;
	bce->Cell = bdce; ceb->Cell = bdce; ebc->Cell = bdce;
	bed->Cell = bdce; edb->Cell = bdce; dbe->Cell = bdce;
	cde->Cell = bdce; dec->Cell = bdce; ecd->Cell = bdce;

	St.push_back(abc); abc->Flag = 1;
	St.push_back(acd); acd->Flag = 1;
	St.push_back(adb); adb->Flag = 1;
	St.push_back(bdc); bdc->Flag = 1;

	bool b0, b1, b2;
	while (St.size())
	{
		abc = St[St.size() - 1]; St.pop_back(); abc->Flag &= 2;
		if (abc->Flag) { continue; }
		a = abc->A; b = abc->B; c = abc->C; d = abc->Rev->Fnx->C;
		if (!insphere(abc, d)) continue;
		b0 = positive(abc->Fnx->Rev, d);
		b1 = positive(abc->Enx->Fnx->Rev, d);
		b2 = positive(abc->Enx->Enx->Fnx->Rev, d);
		if (b0&&b1&&b2) { flip23(abc); continue; }

		if (b1&&b2 && !b0)
		{
			if (abc->Fnx->Fnx->C == d) flip32(abc);
			continue;
		}
		if (b0&&b2 && !b1)
		{
			abc = abc->Enx;
			if (abc->Fnx->Fnx->C == d) flip32(abc);
			continue;
		}
		if (b0&&b1 && !b2)
		{
			abc = abc->Enx->Enx;
			if (abc->Fnx->Fnx->C == d) flip32(abc);
			continue;
		}
		if (b0 && !b1 && !b2)
		{
			if (abc->Enx->Fnx->Fnx->C == d && abc->Enx->Enx->Fnx->Fnx->C == d) flip41(abc);
			continue;
		}
		if (b1 && !b0 && !b2)
		{
			abc = abc->Enx;
			if (abc->Enx->Fnx->Fnx->C == d && abc->Enx->Enx->Fnx->Fnx->C == d) flip41(abc);
			continue;
		}
		if (b2 && !b0 && !b1)
		{
			abc = abc->Enx->Enx;
			if (abc->Enx->Fnx->Fnx->C == d && abc->Enx->Enx->Fnx->Fnx->C == d) flip41(abc);
			continue;
		}
	}
}

inline void CD3DW::flip32(CEdge *abc)
{
	N32++;
	CXYZW *a, *b, *c, *d, *e;
	EdgePtr bca, cab, bac, acb, cba;
	bca = abc->Enx; cab = bca->Enx;
	bac = abc->Rev; acb = bac->Enx; cba = acb->Enx;

	EdgePtr abe, bea, eab, bae, aeb, eba;
	abe = abc->Fnx; bea = abe->Enx; eab = bea->Enx;
	bae = abe->Rev; aeb = bae->Enx; eba = aeb->Enx;
	EdgePtr bce, ceb, ebc, cbe, bec, ecb;
	bce = bca->Fnx; ceb = bce->Enx; ebc = ceb->Enx;
	cbe = bce->Rev; bec = cbe->Enx; ecb = bec->Enx;
	EdgePtr cae, aec, eca, ace, cea, eac;
	cae = cab->Fnx; aec = cae->Enx; eca = aec->Enx;
	ace = cae->Rev; cea = ace->Enx; eac = cea->Enx;

	EdgePtr bad, adb, dba, abd, bda, dab;
	bad = bac->Fnx; adb = bad->Enx; dba = adb->Enx;
	abd = bad->Rev; bda = abd->Enx; dab = bda->Enx;
	EdgePtr acd, cda, dac, cad, adc, dca;
	acd = acb->Fnx; cda = acd->Enx; dac = cda->Enx;
	cad = acd->Rev; adc = cad->Enx; dca = adc->Enx;
	EdgePtr cbd, bdc, dcb, bcd, cdb, dbc;
	cbd = cba->Fnx; bdc = cbd->Enx; dcb = bdc->Enx;
	bcd = cbd->Rev; cdb = bcd->Enx; dbc = cdb->Enx;

	EdgePtr  ead, ade, dea, aed, eda, dae;
	ead = eab->Fnx; ade = ead->Enx; dea = ade->Enx;
	aed = ead->Rev; eda = aed->Enx; dae = eda->Enx;
	EdgePtr dbe, bed, edb, bde, deb, ebd;
	dbe = dba->Fnx; bed = dbe->Enx; edb = bed->Enx;
	bde = dbe->Rev; deb = bde->Enx; ebd = deb->Enx;

	a = abc->A; b = abc->B; c = abc->C; d = dba->A; e = eba->A;

	if (deb->Fnx != dea) { return; }

	EdgePtr cde, dec, ecd, dce, edc, ced;
	new6(cde, dec, ecd, dce, edc, ced, c, d, e);

	cde->Fnx = cda; cdb->Fnx = cde;
	dce->Fnx = dcb; dca->Fnx = dce;

	dec->Fnx = dea; deb->Fnx = dec;
	edc->Fnx = edb; eda->Fnx = edc;

	ecd->Fnx = eca; ecb->Fnx = ecd;
	ced->Fnx = ceb; cea->Fnx = ced;

	cad->Fnx = cae; dae->Fnx = dac; eac->Fnx = ead;
	ace->Fnx = acd; adc->Fnx = ade; aed->Fnx = aec;

	cbe->Fnx = cbd; dbc->Fnx = dbe; ebd->Fnx = ebc;
	bcd->Fnx = bce; bde->Fnx = bdc; bec->Fnx = bed;

	CCell *abce, *acbd, *abed, *cdea, *cedb;

	abce = abc->Cell; acbd = cba->Cell; abed = abe->Cell;
	cdea = new CCell; Cs.push_back(cdea);
	cdea->abc = cde; cdea->adb = cad; cdea->acd = cea; cdea->bdc = dae;
	cdea->VV[0] = c; cdea->VV[1] = d; cdea->VV[2] = e; cdea->VV[3] = a;

	cedb = new CCell; Cs.push_back(cedb);
	cedb->abc = ced; cedb->adb = cbe; cedb->acd = cdb; cedb->bdc = ebd;
	cedb->VV[0] = c; cedb->VV[1] = e; cedb->VV[2] = d; cedb->VV[3] = b;

	abce->CC[0] = cdea; abce->CC[1] = cedb;
	acbd->CC[0] = cdea; acbd->CC[1] = cedb;
	abed->CC[0] = cdea; abed->CC[1] = cedb;

	cde->Cell = cdea; dec->Cell = cdea; ecd->Cell = cdea;
	ace->Cell = cdea; cea->Cell = cdea; eac->Cell = cdea;
	adc->Cell = cdea; dca->Cell = cdea; cad->Cell = cdea;
	aed->Cell = cdea; eda->Cell = cdea; dae->Cell = cdea;

	ced->Cell = cedb; edc->Cell = cedb; dce->Cell = cedb;
	bcd->Cell = cedb; cdb->Cell = cedb; dbc->Cell = cedb;
	bde->Cell = cedb; deb->Cell = cedb; ebd->Cell = cedb;
	bec->Cell = cedb; ecb->Cell = cedb; cbe->Cell = cedb;

	abc->Flag |= 2;
	bca->Flag |= 2;
	cab->Flag |= 2;
	bac->Flag |= 2;
	acb->Flag |= 2;
	cba->Flag |= 2;

	abd->Flag |= 2;
	bda->Flag |= 2;
	dab->Flag |= 2;
	bad->Flag |= 2;
	adb->Flag |= 2;
	dba->Flag |= 2;

	abe->Flag |= 2;
	bea->Flag |= 2;
	eab->Flag |= 2;
	bae->Flag |= 2;
	aeb->Flag |= 2;
	eba->Flag |= 2;

	St.push_back(adc); adc->Flag = 1;
	St.push_back(bcd); bcd->Flag = 1;
}

inline void CD3DW::flip23(CEdge *abc)
{
	N23++;
	CXYZW *a, *b, *c, *d, *e;
	EdgePtr bca, cab, bac, acb, cba;
	bca = abc->Enx; cab = bca->Enx;
	bac = abc->Rev; acb = bac->Enx; cba = acb->Enx;

	EdgePtr abe, bea, eab, bae, aeb, eba;
	abe = abc->Fnx; bea = abe->Enx; eab = bea->Enx;
	bae = abe->Rev; aeb = bae->Enx; eba = aeb->Enx;
	EdgePtr bce, ceb, ebc, cbe, bec, ecb;
	bce = bca->Fnx; ceb = bce->Enx; ebc = ceb->Enx;
	cbe = bce->Rev; bec = cbe->Enx; ecb = bec->Enx;
	EdgePtr cae, aec, eca, ace, cea, eac;
	cae = cab->Fnx; aec = cae->Enx; eca = aec->Enx;
	ace = cae->Rev; cea = ace->Enx; eac = cea->Enx;

	EdgePtr bad, adb, dba, abd, bda, dab;
	bad = bac->Fnx; adb = bad->Enx; dba = adb->Enx;
	abd = bad->Rev; bda = abd->Enx; dab = bda->Enx;
	EdgePtr acd, cda, dac, cad, adc, dca;
	acd = acb->Fnx; cda = acd->Enx; dac = cda->Enx;
	cad = acd->Rev; adc = cad->Enx; dca = adc->Enx;
	EdgePtr cbd, bdc, dcb, bcd, cdb, dbc;
	cbd = cba->Fnx; bdc = cbd->Enx; dcb = bdc->Enx;
	bcd = cbd->Rev; cdb = bcd->Enx; dbc = cdb->Enx;

	a = abc->A; b = abc->B; c = abc->C; d = dba->A; e = eba->A;

	EdgePtr dea, ead, ade, eda, aed, dae;
	new6(dea, ead, ade, eda, aed, dae, d, e, a);
	EdgePtr deb, ebd, bde, edb, bed, dbe;
	new6(deb, ebd, bde, edb, bed, dbe, d, e, b);
	EdgePtr dec, ecd, cde, edc, ced, dce;
	new6(dec, ecd, cde, edc, ced, dce, d, e, c);

	dea->Fnx = deb; deb->Fnx = dec; dec->Fnx = dea;
	eda->Fnx = edc; edc->Fnx = edb; edb->Fnx = eda;

	eac->Fnx = ead; ead->Fnx = eab;
	aeb->Fnx = aed; aed->Fnx = aec;

	eba->Fnx = ebd; ebd->Fnx = ebc;
	bec->Fnx = bed; bed->Fnx = bea;

	ecb->Fnx = ecd; ecd->Fnx = eca;
	cea->Fnx = ced; ced->Fnx = ceb;


	adc->Fnx = ade; ade->Fnx = adb;
	dab->Fnx = dae; dae->Fnx = dac;

	bda->Fnx = bde; bde->Fnx = bdc;
	dbc->Fnx = dbe; dbe->Fnx = dba;

	cdb->Fnx = cde; cde->Fnx = cda;
	dca->Fnx = dce; dce->Fnx = dcb;

	abd->Fnx = abe; bae->Fnx = bad;
	bcd->Fnx = bce; cbe->Fnx = cbd;
	cad->Fnx = cae; ace->Fnx = acd;

	CCell *abce, *acbd, *deab, *debc, *deca;
	abce = abc->Cell; acbd = cba->Cell;

	deab = new CCell; Cs.push_back(deab);
	deab->abc = dea; deab->adb = dbe; deab->acd = dab; deab->bdc = eba;
	deab->VV[0] = d; deab->VV[1] = e; deab->VV[2] = a; deab->VV[3] = b;

	debc = new CCell; Cs.push_back(debc);
	debc->abc = deb; debc->adb = dce; debc->acd = dbc; debc->bdc = ecb;
	debc->VV[0] = d; debc->VV[1] = e; debc->VV[2] = b; debc->VV[3] = c;

	deca = new CCell; Cs.push_back(deca);
	deca->abc = dec; deca->adb = dae; deca->acd = dca; deca->bdc = eac;
	deca->VV[0] = d; deca->VV[1] = e; deca->VV[2] = c; deca->VV[3] = a;

	abce->CC[0] = deab; abce->CC[1] = debc; abce->CC[2] = deca;
	acbd->CC[0] = deab; acbd->CC[1] = debc; acbd->CC[2] = deca;

	abd->Cell = deab; bda->Cell = deab; dab->Cell = deab;
	bcd->Cell = debc; cdb->Cell = debc; dbc->Cell = debc;
	cad->Cell = deca; adc->Cell = deca; dca->Cell = deca;

	eba->Cell = deab; bae->Cell = deab; aeb->Cell = deab;
	ecb->Cell = debc; cbe->Cell = debc; bec->Cell = debc;
	eac->Cell = deca; ace->Cell = deca; cea->Cell = deca;

	dea->Cell = deab; ead->Cell = deab; ade->Cell = deab;
	deb->Cell = debc; ebd->Cell = debc; bde->Cell = debc;
	dec->Cell = deca; ecd->Cell = deca; cde->Cell = deca;

	eda->Cell = deca; dae->Cell = deca; aed->Cell = deca;
	edb->Cell = deab; dbe->Cell = deab; bed->Cell = deab;
	edc->Cell = debc; dce->Cell = debc; ced->Cell = debc;

	abc->Flag |= 2; bca->Flag |= 2; cab->Flag |= 2;
	cba->Flag |= 2; bac->Flag |= 2; acb->Flag |= 2;

	St.push_back(abd); abd->Flag = 1;
	St.push_back(bcd); bcd->Flag = 1;
	St.push_back(cad); cad->Flag = 1;
}

inline void CD3DW::flip41(CEdge *abc)
{
	N41++;
	CXYZW *a, *b, *d, *e;
	EdgePtr bca, cab, bac, acb, cba;
	bca = abc->Enx; cab = bca->Enx;
	bac = abc->Rev; acb = bac->Enx; cba = acb->Enx;

	EdgePtr abe, bea, eab, bae, aeb, eba;
	abe = abc->Fnx; bea = abe->Enx; eab = bea->Enx;
	bae = abe->Rev; aeb = bae->Enx; eba = aeb->Enx;
	EdgePtr bce, ceb, ebc, cbe, bec, ecb;
	bce = bca->Fnx; ceb = bce->Enx; ebc = ceb->Enx;
	cbe = bce->Rev; bec = cbe->Enx; ecb = bec->Enx;
	EdgePtr cae, aec, eca, ace, cea, eac;
	cae = cab->Fnx; aec = cae->Enx; eca = aec->Enx;
	ace = cae->Rev; cea = ace->Enx; eac = cea->Enx;

	EdgePtr bad, adb, dba, abd, bda, dab;
	bad = bac->Fnx; adb = bad->Enx; dba = adb->Enx;
	abd = bad->Rev; bda = abd->Enx; dab = bda->Enx;
	EdgePtr acd, cda, dac, cad, adc, dca;
	acd = acb->Fnx; cda = acd->Enx; dac = cda->Enx;
	cad = acd->Rev; adc = cad->Enx; dca = adc->Enx;
	EdgePtr cbd, bdc, dcb, bcd, cdb, dbc;
	cbd = cba->Fnx; bdc = cbd->Enx; dcb = bdc->Enx;
	bcd = cbd->Rev; cdb = bcd->Enx; dbc = cdb->Enx;

	EdgePtr dae, aed, eda, ade, dea, ead;
	dae = dac->Fnx; aed = dae->Enx; eda = aed->Enx;
	ade = dae->Rev; dea = ade->Enx; ead = dea->Enx;
	EdgePtr bde, deb, ebd, dbe, bed, edb;
	bde = bdc->Fnx; deb = bde->Enx; ebd = deb->Enx;
	dbe = bde->Rev; bed = dbe->Enx; edb = bed->Enx;
	EdgePtr dec, ecd, cde, edc, dce, ced;
	dec = dea->Fnx; ecd = dec->Enx; cde = ecd->Enx;
	edc = dec->Rev; dce = edc->Enx; ced = dce->Enx;

	a = abc->A; b = abc->B; d = dba->A; e = eba->A;

	abd->Fnx = abe; bda->Fnx = bde; dab->Fnx = dae;
	aeb->Fnx = aed; eba->Fnx = ebd; bae->Fnx = bad;
	ade->Fnx = adb; dea->Fnx = deb; ead->Fnx = eab;
	bed->Fnx = bea; edb->Fnx = eda; dbe->Fnx = dba;

	CCell *abdc, *baec, *adec, *bedc, *abde;
	abde = new CCell; Cs.push_back(abde);
	abde->abc = abd; abde->adb = aeb; abde->acd = ade; abde->bdc = bed;
	abde->VV[0] = a; abde->VV[1] = b; abde->VV[2] = d; abde->VV[3] = e;

	abdc = abd->Cell; abdc->CC[0] = abde;
	baec = bae->Cell; baec->CC[0] = abde;
	adec = ade->Cell; adec->CC[0] = abde;
	bedc = bed->Cell; bedc->CC[0] = abde;

	abd->Cell = abde; bda->Cell = abde; dab->Cell = abde;
	bae->Cell = abde; aeb->Cell = abde; eba->Cell = abde;
	ade->Cell = abde; dea->Cell = abde; ead->Cell = abde;
	dbe->Cell = abde; bed->Cell = abde; edb->Cell = abde;

	abc->Flag |= 2;
	bca->Flag |= 2;
	cab->Flag |= 2;
	bac->Flag |= 2;
	acb->Flag |= 2;
	cba->Flag |= 2;

	adc->Flag |= 2;
	dca->Flag |= 2;
	cad->Flag |= 2;
	dac->Flag |= 2;
	acd->Flag |= 2;
	cda->Flag |= 2;

	bdc->Flag |= 2;
	dcb->Flag |= 2;
	cbd->Flag |= 2;
	dbc->Flag |= 2;
	bcd->Flag |= 2;
	cdb->Flag |= 2;

	ace->Flag |= 2;
	cea->Flag |= 2;
	eac->Flag |= 2;
	cae->Flag |= 2;
	aec->Flag |= 2;
	eca->Flag |= 2;

	bce->Flag |= 2;
	ceb->Flag |= 2;
	ebc->Flag |= 2;
	cbe->Flag |= 2;
	bec->Flag |= 2;
	ecb->Flag |= 2;

	dce->Flag |= 2;
	ced->Flag |= 2;
	edc->Flag |= 2;
	cde->Flag |= 2;
	dec->Flag |= 2;
	ecd->Flag |= 2;
}

inline void CD3DW::new6(EdgePtr&abc,EdgePtr&bca,EdgePtr&cab,
	EdgePtr&bac,EdgePtr&cba,EdgePtr&acb,CXYZW*a,CXYZW*b,CXYZW*c)
{
	abc = new CEdge(a, b, c); Es.push_back(abc);
	bca = new CEdge(b, c, a); Es.push_back(bca);
	cab = new CEdge(c, a, b); Es.push_back(cab);
	bac = new CEdge(b, a, c); Es.push_back(bac);
	cba = new CEdge(c, b, a); Es.push_back(cba);
	acb = new CEdge(a, c, b); Es.push_back(acb);
	abc->Rev = bac; bac->Rev = abc;
	bca->Rev = cba; cba->Rev = bca;
	cab->Rev = acb; acb->Rev = cab;
	abc->Enx = bca; bca->Enx = cab; cab->Enx = abc;
	bac->Enx = acb; acb->Enx = cba; cba->Enx = bac;
}

inline bool CD3DW::inside(CCell *C,CXYZW *V)
{
	return positive(C->abc, V)
		&& positive(C->adb, V)
		&& positive(C->acd, V)
		&& positive(C->bdc, V);
}

inline CEdge* CD3DW::locate(CXYZW *V)
{
	CCell *C = Cs[0];
	while (C->CC[0])
	{
		NLoc++;
		if (!C->CC[1]) { C = C->CC[0]; continue; }
		if (inside(C->CC[0], V)) { C = C->CC[0]; continue; }
		if (!C->CC[2]) { C = C->CC[1]; continue; }
		if (inside(C->CC[1], V)) { C = C->CC[1]; continue; }
		if (!C->CC[3]) { C = C->CC[2]; continue; }
		if (inside(C->CC[2], V)) { C = C->CC[2]; continue; }
		C = C->CC[3];
	}
	return C->abc;
}

inline void CEdge::check_minors()
{
	if (Minors) return;
	CInexact M01 = B->X - A->X;
	CInexact M02 = B->Y - A->Y;
	CInexact M03 = B->Z - A->Z;
	CInexact M12 = A->X * B->Y - B->X * A->Y;
	CInexact M13 = A->X * B->Z - B->X * A->Z;
	CInexact M23 = A->Y * B->Z - B->Y * A->Z;

	M012 = M12 - C->X * M02 + C->Y * M01;
	M013 = M13 - C->X * M03 + C->Z * M01;
	M023 = M23 - C->Y * M03 + C->Z * M02;
	M123 = C->X * M23 - C->Y * M13 + C->Z * M12;

	Enx->M012 = M012;
	Enx->M013 = M013;
	Enx->M023 = M023;
	Enx->M123 = M123;

	Enx->Enx->M012 = M012;
	Enx->Enx->M013 = M013;
	Enx->Enx->M023 = M023;
	Enx->Enx->M123 = M123;

	Rev->M012 = -M012;
	Rev->M013 = -M013;
	Rev->M023 = -M023;
	Rev->M123 = -M123;

	Rev->Enx->M012 = -M012;
	Rev->Enx->M013 = -M013;
	Rev->Enx->M023 = -M023;
	Rev->Enx->M123 = -M123;

	Rev->Enx->Enx->M012 = -M012;
	Rev->Enx->Enx->M013 = -M013;
	Rev->Enx->Enx->M023 = -M023;
	Rev->Enx->Enx->M123 = -M123;

	Minors = true;
	Enx->Minors = true;
	Enx->Enx->Minors = true;
	Rev->Minors = true;
	Rev->Enx->Minors = true;
	Rev->Enx->Enx->Minors = true;
}

inline bool CD3DW::positive(CEdge *E,CXYZW *V)
{
	bool valid = false;
	
	Ntp++;
	if (E->LastV == V) { Nrp++; return E->LastP; }

	if (E->Minors) Ndd++;
	E->check_minors();

	CInexact det = E->M123 - V->X * E->M023 + V->Y * E->M013 - V->Z * E->M012;

	if (det.is_zero()) { valid = e_positive(*(E->A),*(E->B),*(E->C),*V); }
	else { valid = (0<(double)det); }

	E->LastP = valid; E->LastV = V; E = E->Enx;
	E->LastP = valid; E->LastV = V; E = E->Enx;
	E->LastP = valid; E->LastV = V; E = E->Rev;
	E->LastP = !valid; E->LastV = V; E = E->Enx;
	E->LastP = !valid; E->LastV = V; E = E->Enx;
	E->LastP = !valid; E->LastV = V;

	return valid;
}

bool CD3DW::insphere(CEdge *E,CXYZW *V)
{
	Nti++;
	CXYZW *A = E->A;
	CXYZW *B = E->B;
	CXYZW *C = E->C;
	CXYZW *D = E->Fnx->C;

	E->check_minors();
	CInexact M4 = E->M023 * D->X - E->M013 * D->Y + E->M012 * D->Z - E->M123;
	CInexact M3 = E->M023 * V->X - E->M013 * V->Y + E->M012 * V->Z - E->M123;
	CEdge*H = E->Fnx; H->check_minors();	// abd
	CInexact M2 = H->M023 * V->X - H->M013 * V->Y + H->M012 * V->Z - H->M123;
	H = E->Enx->Enx->Fnx->Rev; H->check_minors();
	CInexact M1 = H->M023 * V->X - H->M013 * V->Y + H->M012 * V->Z - H->M123;
	H = E->Enx->Fnx; H->check_minors();	// bcd
	CInexact M0 = H->M023 * V->X - H->M013 * V->Y + H->M012 * V->Z - H->M123;

	CInexact det = A->U * M0 - B->U * M1 + C->U * M2 - D->U * M3 + V->U * M4;

	if (det.is_zero()) { return e_insphere(*A,*B,*C,*D,*V); }
	return 0 < (double)det;
}

inline bool CD3DW::e_positive(CXYZW&A, CXYZW&B, CXYZW&C, CXYZW&D)
{
	int i, j;
	bool b = true;
	CXYZW *p;
	CXYZW *P[4]; P[0] = &A; P[1] = &B; P[2] = &C; P[3] = &D;
	Nep++;
	for (i = 1;i < 4;i++)
	{
		for (j = i;j > 0;j--)
		{
			if (P[j]->N > P[j - 1]->N) break;
			p = P[j - 1]; P[j - 1] = P[j]; P[j] = p; b = !b;
		}
	}
	return sign4(P) == b;
}

inline bool CD3DW::e_insphere(CXYZW&A, CXYZW&B, CXYZW&C, CXYZW&D, CXYZW&E)
{
	int i, j;
	bool b = true;
	CXYZW *p;
	CXYZW *P[5]; P[0] = &A; P[1] = &B; P[2] = &C; P[3] = &D; P[4] = &E;
	Nei++;
	for (i = 1;i<5;i++) for (j = i;j>0;j--)
	{
		if (P[j]->N>P[j - 1]->N) break;
		p = P[j - 1]; P[j - 1] = P[j]; P[j] = p; b = !b;
	}
	return sign5(P) == b;
}

inline bool CD3DW::sign4(CXYZW**P)
{
	CExact m23x23 = P[2]->eZ - P[3]->eZ;
	CExact m13x23 = P[1]->eZ - P[3]->eZ;
	CExact m12x23 = P[1]->eZ - P[2]->eZ;
	CExact m123x123 = P[1]->eY*m23x23 - P[2]->eY*m13x23 + P[3]->eY*m12x23;
	CExact m03x23 = P[0]->eZ - P[3]->eZ;
	CExact m02x23 = P[0]->eZ - P[2]->eZ;
	CExact m023x123 = P[0]->eY*m23x23 - P[2]->eY*m03x23 + P[3]->eY*m02x23;
	CExact m01x23 = P[0]->eZ - P[1]->eZ;
	CExact m013x123 = P[0]->eY*m13x23 - P[1]->eY*m03x23 + P[3]->eY*m01x23;
	CExact m012x123 = P[0]->eY*m12x23 - P[1]->eY*m02x23 + P[2]->eY*m01x23;
	CExact m0123x0123 = P[0]->eX*m123x123 - P[1]->eX*m023x123 + P[2]->eX*m013x123 - P[3]->eX*m012x123;
	if ((bool)m0123x0123) return m0123x0123.Positive();

	CExact m23x13 = P[2]->eY - P[3]->eY;
	CExact m13x13 = P[1]->eY - P[3]->eY;
	CExact m12x13 = P[1]->eY - P[2]->eY;
	CExact m123x013 = P[1]->eX*m23x13 - P[2]->eX*m13x13 + P[3]->eX*m12x13;
	if ((bool)m123x013) return m123x013.Positive();

	CExact m123x023 = P[1]->eX*m23x23 - P[2]->eX*m13x23 + P[3]->eX*m12x23;
	if ((bool)m123x023) return !m123x023.Positive();

	if ((bool)m123x123) return m123x123.Positive();

	CExact m03x13 = P[0]->eY - P[3]->eY;
	CExact m02x13 = P[0]->eY - P[2]->eY;
	CExact m023x013 = P[0]->eX*m23x13 - P[2]->eX*m03x13 + P[3]->eX*m02x13;
	if ((bool)m023x013) return !m023x013.Positive();

	CExact m23x03 = P[2]->eX - P[3]->eX;
	if ((bool)m23x03) return m23x03.Positive();

	if ((bool)m23x13) return !m23x13.Positive();

	CExact m023x023 = P[0]->eX*m23x23 - P[2]->eX*m03x23 + P[3]->eX*m02x23;
	if ((bool)m023x023) return m023x023.Positive();

	if ((bool)m23x23) return m23x23.Positive();

	if ((bool)m023x123) return !m023x123.Positive();

	CExact m01x13 = P[0]->eY - P[1]->eY;
	CExact m013x013 = P[0]->eX*m13x13 - P[1]->eX*m03x13 + P[3]->eX*m01x13;
	if ((bool)m013x013) return m013x013.Positive();

	CExact m13x03 = P[1]->eX - P[3]->eX;
	if ((bool)m13x03) return !m13x03.Positive();

	if ((bool)m13x13) return m13x13.Positive();

	CExact m03x03 = P[0]->eX - P[3]->eX;
	if ((bool)m03x03) return m03x03.Positive();

	return true;
}

inline bool CD3DW::sign5(CXYZW**P)
{
	CExact m34x34 = P[3]->eU - P[4]->eU;
	CExact m24x34 = P[2]->eU - P[4]->eU;
	CExact m23x34 = P[2]->eU - P[3]->eU;
	CExact m234x234 = P[2]->eZ*m34x34 - P[3]->eZ*m24x34 + P[4]->eZ*m23x34;
	CExact m14x34 = P[1]->eU - P[4]->eU;
	CExact m13x34 = P[1]->eU - P[3]->eU;
	CExact m134x234 = P[1]->eZ*m34x34 - P[3]->eZ*m14x34 + P[4]->eZ*m13x34;
	CExact m12x34 = P[1]->eU - P[2]->eU;
	CExact m124x234 = P[1]->eZ*m24x34 - P[2]->eZ*m14x34 + P[4]->eZ*m12x34;
	CExact m123x234 = P[1]->eZ*m23x34 - P[2]->eZ*m13x34 + P[3]->eZ*m12x34;
	CExact m1234x1234 = P[1]->eY*m234x234 - P[2]->eY*m134x234 + P[3]->eY*m124x234 - P[4]->eY*m123x234;
	CExact m04x34 = P[0]->eU - P[4]->eU;
	CExact m03x34 = P[0]->eU - P[3]->eU;
	CExact m034x234 = P[0]->eZ*m34x34 - P[3]->eZ*m04x34 + P[4]->eZ*m03x34;
	CExact m02x34 = P[0]->eU - P[2]->eU;
	CExact m024x234 = P[0]->eZ*m24x34 - P[2]->eZ*m04x34 + P[4]->eZ*m02x34;
	CExact m023x234 = P[0]->eZ*m23x34 - P[2]->eZ*m03x34 + P[3]->eZ*m02x34;
	CExact m0234x1234 = P[0]->eY*m234x234 - P[2]->eY*m034x234 + P[3]->eY*m024x234 - P[4]->eY*m023x234;
	CExact m01x34 = P[0]->eU - P[1]->eU;
	CExact m014x234 = P[0]->eZ*m14x34 - P[1]->eZ*m04x34 + P[4]->eZ*m01x34;
	CExact m013x234 = P[0]->eZ*m13x34 - P[1]->eZ*m03x34 + P[3]->eZ*m01x34;
	CExact m0134x1234 = P[0]->eY*m134x234 - P[1]->eY*m034x234 + P[3]->eY*m014x234 - P[4]->eY*m013x234;
	CExact m012x234 = P[0]->eZ*m12x34 - P[1]->eZ*m02x34 + P[2]->eZ*m01x34;
	CExact m0124x1234 = P[0]->eY*m124x234 - P[1]->eY*m024x234 + P[2]->eY*m014x234 - P[4]->eY*m012x234;
	CExact m0123x1234 = P[0]->eY*m123x234 - P[1]->eY*m023x234 + P[2]->eY*m013x234 - P[3]->eY*m012x234;
	CExact m01234x01234 = P[0]->eX*m1234x1234 - P[1]->eX*m0234x1234 + P[2]->eX*m0134x1234 - P[3]->eX*m0124x1234 + P[4]->eX*m0123x1234;
	if ((bool)m01234x01234) return m01234x01234.Positive();

	CExact m34x24 = P[3]->eZ - P[4]->eZ;
	CExact m24x24 = P[2]->eZ - P[4]->eZ;
	CExact m23x24 = P[2]->eZ - P[3]->eZ;
	CExact m234x124 = P[2]->eY*m34x24 - P[3]->eY*m24x24 + P[4]->eY*m23x24;
	CExact m14x24 = P[1]->eZ - P[4]->eZ;
	CExact m13x24 = P[1]->eZ - P[3]->eZ;
	CExact m134x124 = P[1]->eY*m34x24 - P[3]->eY*m14x24 + P[4]->eY*m13x24;
	CExact m12x24 = P[1]->eZ - P[2]->eZ;
	CExact m124x124 = P[1]->eY*m24x24 - P[2]->eY*m14x24 + P[4]->eY*m12x24;
	CExact m123x124 = P[1]->eY*m23x24 - P[2]->eY*m13x24 + P[3]->eY*m12x24;
	CExact m1234x0124 = P[1]->eX*m234x124 - P[2]->eX*m134x124 + P[3]->eX*m124x124 - P[4]->eX*m123x124;
	if ((bool)m1234x0124) return !m1234x0124.Positive();

	CExact m234x134 = P[2]->eY*m34x34 - P[3]->eY*m24x34 + P[4]->eY*m23x34;
	CExact m134x134 = P[1]->eY*m34x34 - P[3]->eY*m14x34 + P[4]->eY*m13x34;
	CExact m124x134 = P[1]->eY*m24x34 - P[2]->eY*m14x34 + P[4]->eY*m12x34;
	CExact m123x134 = P[1]->eY*m23x34 - P[2]->eY*m13x34 + P[3]->eY*m12x34;
	CExact m1234x0134 = P[1]->eX*m234x134 - P[2]->eX*m134x134 + P[3]->eX*m124x134 - P[4]->eX*m123x134;
	if ((bool)m1234x0134) return m1234x0134.Positive();

	CExact m1234x0234 = P[1]->eX*m234x234 - P[2]->eX*m134x234 + P[3]->eX*m124x234 - P[4]->eX*m123x234;
	if ((bool)m1234x0234) return !m1234x0234.Positive();

	if ((bool)m1234x1234) return m1234x1234.Positive();

	CExact m04x24 = P[0]->eZ - P[4]->eZ;
	CExact m03x24 = P[0]->eZ - P[3]->eZ;
	CExact m034x124 = P[0]->eY*m34x24 - P[3]->eY*m04x24 + P[4]->eY*m03x24;
	CExact m02x24 = P[0]->eZ - P[2]->eZ;
	CExact m024x124 = P[0]->eY*m24x24 - P[2]->eY*m04x24 + P[4]->eY*m02x24;
	CExact m023x124 = P[0]->eY*m23x24 - P[2]->eY*m03x24 + P[3]->eY*m02x24;
	CExact m0234x0124 = P[0]->eX*m234x124 - P[2]->eX*m034x124 + P[3]->eX*m024x124 - P[4]->eX*m023x124;
	if ((bool)m0234x0124) return m0234x0124.Positive();

	CExact m34x14 = P[3]->eY - P[4]->eY;
	CExact m24x14 = P[2]->eY - P[4]->eY;
	CExact m23x14 = P[2]->eY - P[3]->eY;
	CExact m234x014 = P[2]->eX*m34x14 - P[3]->eX*m24x14 + P[4]->eX*m23x14;
	if ((bool)m234x014) return m234x014.Positive();

	CExact m234x024 = P[2]->eX*m34x24 - P[3]->eX*m24x24 + P[4]->eX*m23x24;
	if ((bool)m234x024) return !m234x024.Positive();

	if ((bool)m234x124) return m234x124.Positive();

	CExact m034x134 = P[0]->eY*m34x34 - P[3]->eY*m04x34 + P[4]->eY*m03x34;
	CExact m024x134 = P[0]->eY*m24x34 - P[2]->eY*m04x34 + P[4]->eY*m02x34;
	CExact m023x134 = P[0]->eY*m23x34 - P[2]->eY*m03x34 + P[3]->eY*m02x34;
	CExact m0234x0134 = P[0]->eX*m234x134 - P[2]->eX*m034x134 + P[3]->eX*m024x134 - P[4]->eX*m023x134;
	if ((bool)m0234x0134) return !m0234x0134.Positive();

	CExact m234x034 = P[2]->eX*m34x34 - P[3]->eX*m24x34 + P[4]->eX*m23x34;
	if ((bool)m234x034) return m234x034.Positive();

	if ((bool)m234x134) return !m234x134.Positive();

	CExact m0234x0234 = P[0]->eX*m234x234 - P[2]->eX*m034x234 + P[3]->eX*m024x234 - P[4]->eX*m023x234;
	if ((bool)m0234x0234) return m0234x0234.Positive();

	if ((bool)m234x234) return m234x234.Positive();

	if ((bool)m0234x1234) return !m0234x1234.Positive();

	CExact m01x24 = P[0]->eZ - P[1]->eZ;
	CExact m014x124 = P[0]->eY*m14x24 - P[1]->eY*m04x24 + P[4]->eY*m01x24;
	CExact m013x124 = P[0]->eY*m13x24 - P[1]->eY*m03x24 + P[3]->eY*m01x24;
	CExact m0134x0124 = P[0]->eX*m134x124 - P[1]->eX*m034x124 + P[3]->eX*m014x124 - P[4]->eX*m013x124;
	if ((bool)m0134x0124) return !m0134x0124.Positive();

	CExact m14x14 = P[1]->eY - P[4]->eY;
	CExact m13x14 = P[1]->eY - P[3]->eY;
	CExact m134x014 = P[1]->eX*m34x14 - P[3]->eX*m14x14 + P[4]->eX*m13x14;
	if ((bool)m134x014) return !m134x014.Positive();

	CExact m134x024 = P[1]->eX*m34x24 - P[3]->eX*m14x24 + P[4]->eX*m13x24;
	if ((bool)m134x024) return m134x024.Positive();

	if ((bool)m134x124) return !m134x124.Positive();

	CExact m04x14 = P[0]->eY - P[4]->eY;
	CExact m03x14 = P[0]->eY - P[3]->eY;
	CExact m034x014 = P[0]->eX*m34x14 - P[3]->eX*m04x14 + P[4]->eX*m03x14;
	if ((bool)m034x014) return m034x014.Positive();

	CExact m34x04 = P[3]->eX - P[4]->eX;
	if ((bool)m34x04) return !m34x04.Positive();

	if ((bool)m34x14) return m34x14.Positive();

	CExact m034x024 = P[0]->eX*m34x24 - P[3]->eX*m04x24 + P[4]->eX*m03x24;
	if ((bool)m034x024) return !m034x024.Positive();

	if ((bool)m34x24) return !m34x24.Positive();

	if ((bool)m034x124) return m034x124.Positive();

	CExact m014x134 = P[0]->eY*m14x34 - P[1]->eY*m04x34 + P[4]->eY*m01x34;
	CExact m013x134 = P[0]->eY*m13x34 - P[1]->eY*m03x34 + P[3]->eY*m01x34;
	CExact m0134x0134 = P[0]->eX*m134x134 - P[1]->eX*m034x134 + P[3]->eX*m014x134 - P[4]->eX*m013x134;
	if ((bool)m0134x0134) return m0134x0134.Positive();

	CExact m134x034 = P[1]->eX*m34x34 - P[3]->eX*m14x34 + P[4]->eX*m13x34;
	if ((bool)m134x034) return !m134x034.Positive();

	if ((bool)m134x134) return m134x134.Positive();

	CExact m034x034 = P[0]->eX*m34x34 - P[3]->eX*m04x34 + P[4]->eX*m03x34;
	if ((bool)m034x034) return m034x034.Positive();

	if ((bool)m34x34) return m34x34.Positive();

	if ((bool)m034x134) return !m034x134.Positive();

	CExact m0134x0234 = P[0]->eX*m134x234 - P[1]->eX*m034x234 + P[3]->eX*m014x234 - P[4]->eX*m013x234;
	if ((bool)m0134x0234) return !m0134x0234.Positive();

	if ((bool)m134x234) return !m134x234.Positive();

	if ((bool)m034x234) return m034x234.Positive();

	if ((bool)m0134x1234) return m0134x1234.Positive();

	CExact m012x124 = P[0]->eY*m12x24 - P[1]->eY*m02x24 + P[2]->eY*m01x24;
	CExact m0124x0124 = P[0]->eX*m124x124 - P[1]->eX*m024x124 + P[2]->eX*m014x124 - P[4]->eX*m012x124;
	if ((bool)m0124x0124) return m0124x0124.Positive();

	CExact m12x14 = P[1]->eY - P[2]->eY;
	CExact m124x014 = P[1]->eX*m24x14 - P[2]->eX*m14x14 + P[4]->eX*m12x14;
	if ((bool)m124x014) return m124x014.Positive();

	CExact m124x024 = P[1]->eX*m24x24 - P[2]->eX*m14x24 + P[4]->eX*m12x24;
	if ((bool)m124x024) return !m124x024.Positive();

	if ((bool)m124x124) return m124x124.Positive();

	CExact m02x14 = P[0]->eY - P[2]->eY;
	CExact m024x014 = P[0]->eX*m24x14 - P[2]->eX*m04x14 + P[4]->eX*m02x14;
	if ((bool)m024x014) return !m024x014.Positive();

	CExact m24x04 = P[2]->eX - P[4]->eX;
	if ((bool)m24x04) return m24x04.Positive();

	if ((bool)m24x14) return !m24x14.Positive();

	CExact m024x024 = P[0]->eX*m24x24 - P[2]->eX*m04x24 + P[4]->eX*m02x24;
	if ((bool)m024x024) return m024x024.Positive();

	if ((bool)m24x24) return m24x24.Positive();

	if ((bool)m024x124) return !m024x124.Positive();

	CExact m01x14 = P[0]->eY - P[1]->eY;
	CExact m014x014 = P[0]->eX*m14x14 - P[1]->eX*m04x14 + P[4]->eX*m01x14;
	if ((bool)m014x014) return m014x014.Positive();

	CExact m14x04 = P[1]->eX - P[4]->eX;
	if ((bool)m14x04) return !m14x04.Positive();

	if ((bool)m14x14) return m14x14.Positive();

	CExact m04x04 = P[0]->eX - P[4]->eX;
	if ((bool)m04x04) return m04x04.Positive();

	return true;
}

inline bool edge_less(EdgePtr A, EdgePtr B)
{
	if (A->A->N<B->A->N) return true;
	if (A->A->N>B->A->N) return false;
	return A->B->N<B->B->N;
}

void CD3DW::list_edges(std::vector<EdgePtr>*V)
{
	int i, j;
	V->resize(Es.size());
	for (i = 0;i<Es.size();i++) (*V)[i] = Es[i];
	std::stable_sort(V->begin(), V->end(), edge_less);
	j = 0;
	for (i = 0;i<V->size();i++)
	{
		if ((*V)[i]->Flag & 2) continue;
		if ((*V)[i]->A->N<0) continue;
		if ((*V)[i]->B->N<0) continue;
		if (j && (*V)[i]->A == (*V)[j - 1]->A
			&& (*V)[i]->B == (*V)[j - 1]->B) continue;
		(*V)[j] = (*V)[i]; j++;
	}
	V->resize(j);
}

CComplex::CComplex(CD3DW *D)
{
	CEdge *E,*F;
	C0Simplex *S0;
	C1Simplex *S1;
	C2Simplex *S2;
	C3Simplex *S3;
	for (int i=0;i<(int)D->Es.size();i++) 
	{
		D->Es[i]->S1 = 0; D->Es[i]->S2 = 0; D->Es[i]->S3 = 0;
	}
	for (int i=0;i<(int)D->Es.size();i++)
	{
		E = D->Es[i]; if (E->Flag & 2) continue;
		if (E->A->N<0 || E->B->N<0 || E->C->N<0 || E->Fnx->C->N<0) continue;
		if (E->S3) continue;
		while (E->A->N>E->B->N || E->A->N>E->C->N) E = E->Enx;
		if (E->B->N>E->C->N) continue;
		if (E->A->N<E->Fnx->C->N && E->C->N>E->Fnx->C->N) continue;
		S3 = new C3Simplex;
		S3->E = E; S3->A = E->A; S3->B = E->B; S3->C = E->C; S3->D = E->Fnx->C;
		E->S3 = S3; E->Enx->S3 = S3; E->Enx->Enx->S3 = S3;
		F = E->Fnx->Rev;
		F->S3 = S3; F->Enx->S3 = S3; F->Enx->Enx->S3 = S3;
		F = E->Enx->Fnx->Rev;
		F->S3 = S3; F->Enx->S3 = S3; F->Enx->Enx->S3 = S3;
		F = E->Enx->Enx->Fnx->Rev;
		F->S3 = S3; F->Enx->S3 = S3; F->Enx->Enx->S3 = S3;
		S3->geometry();
		P.push_back(S3);
	}
	for (int i=0;i<(int)D->Es.size();i++)
	{
		E = D->Es[i]; if (E->Flag & 2) continue;
		if (E->S2) continue;
		if (E->A->N<0 || E->B->N<0 || E->C->N<0) continue;
		while (E->A->N>E->B->N || E->A->N>E->C->N) E = E->Enx;
		if (E->B->N>E->C->N) E = E->Rev;
		S2 = new C2Simplex;
		S2->E = E; S2->A = E->A; S2->B = E->B; S2->C = E->C;
		E->S2 = S2; E->Enx->S2 = S2; E->Enx->Enx->S2 = S2;
		E->Rev->S2 = S2; E->Rev->Enx->S2 = S2; E->Rev->Enx->Enx->S2 = S2;
		S2->geometry();
		P.push_back(S2);
	}
	for (int i = 0;i < D->Es.size();i++)
	{
		E = D->Es[i]; if (E->Flag & 2) continue;
		if (E->S1) continue;
		if (E->A->N<0 || E->B->N<0 || E->B->N<E->A->N) continue;
		S1 = new C1Simplex;
		S1->E = E; S1->A = E->A; S1->B = E->B;
		E->S1 = S1; E->Rev->S1 = S1;
		for (F = E->Fnx;F != E;F = F->Fnx) { F->S1 = S1; F->Rev->S1 = S1; }
		S1->geometry();
		P.push_back(S1);

		S0 = S1->A->S0;
		if (!S0)
		{
			S0 = new C0Simplex;
			S0->A = S1->A; S1->A->S0 = S0;
			S0->Si = -S1->A->W; S0->Sii = 0;
			if (S1->At)
			{
				S0->Mu0 = S1->Mu0; S0->Mu0i = S1->Mu0i;
			}
			else
			{
				S0->Mu0 = S1->Si; S0->Mu0i = S1->Sii;
			}
			S0->Mu1 = S1->Si; S0->Mu1i = S1->Sii;
			P.push_back(S0);
		}
		if (S0->Si > S1->Si && !S1->Sii) S0->At = 1;
		if (!S1->At && !S1->Sii && (S0->Mu0i || S1->Si<S0->Mu0)) { S0->Mu0 = S1->Si; S0->Mu0i = 0; }
		if (!S1->Mu0i && (S0->Mu0i || S1->Mu0<S0->Mu0)) { S0->Mu0 = S1->Mu0; S0->Mu0i = 0; }
		if (S1->Sii || S1->Mu1i) S0->Mu1i = 1;
		else
		{
			if (S1->Si>S0->Mu1) S0->Mu1 = S1->Si;
			if (S1->Mu1>S0->Mu1) S0->Mu1 = F->S1->Mu1;
		}

		S0 = S1->B->S0;
		if (!S0)
		{
			S0 = new C0Simplex;
			S0->A = S1->B; S1->B->S0 = S0;
			S0->Si = -S1->B->W; S0->Sii = 0;
			if (S1->At)
			{
				S0->Mu0 = S1->Mu0; S0->Mu0i = S1->Mu0i;
			}
			else
			{
				S0->Mu0 = S1->Si; S0->Mu0i = S1->Sii;
			}
			S0->Mu1 = S1->Si; S0->Mu1i = S1->Sii;
			P.push_back(S0);
		}
		if (S0->Si > S1->Si && !S1->Sii) S0->At = 1;
		if (!S1->At && !S1->Sii && (S0->Mu0i || S1->Si<S0->Mu0)) { S0->Mu0 = S1->Si; S0->Mu0i = 0; }
		if (!S1->Mu0i && (S0->Mu0i || S1->Mu0<S0->Mu0)) { S0->Mu0 = S1->Mu0; S0->Mu0i = 0; }
		if (S1->Sii || S1->Mu1i) S0->Mu1i = 1;
		else
		{
			if (S1->Si>S0->Mu1) S0->Mu1 = S1->Si;
			if (S1->Mu1>S0->Mu1) S0->Mu1 = F->S1->Mu1;
		}
	}
}

CComplex::~CComplex()
{
	int i;
	for (i=0;i<(int)P.size();i++) {delete P[i];P[i] = NULL;}
	P.clear();
}

CXYZW* C0Simplex::operator [](int i)
{
	switch(i)
	{
		case 0: return A;
	}
	return 0;
}

CXYZW* C1Simplex::operator [](int i)
{
	switch(i)
	{
		case 0: return A;
		case 1: return B;
	}
	return 0;
}

CXYZW* C2Simplex::operator [](int i)
{
	switch(i)
	{
		case 0: return A;
		case 1: return B;
		case 2: return C;
	}
	return 0;
}

CXYZW* C3Simplex::operator[](int i)
{
	switch (i)
	{
		case 0: return A;
		case 1: return B;
		case 2: return C;
		case 3: return D;
	}
	return 0;
}

bool C2Simplex::on_surface(double a)
{
	if (At) { if (Mu0 > 1 || Mu0i) return false; }
	else { if (Si > a || Sii) return false; }
	if (Mu1 < a && !Mu1i) { return false; }
	return true;
}

bool C1Simplex::on_surface(double a)
{
	if (At) { if (Mu0 > a || Mu0i) return false; }
	else { if (Si > a || Sii) return false; }
	if (Mu1 < a && !Mu1i) { return false; }
	return true;
}

bool C0Simplex::on_surface(double a)
{
	if (At) { if (Mu0 > a || Mu0i) return false; }
	else { if (Si > a || Sii) return false; }
	if (Mu1 < a && !Mu1i) { return false; }
	return true;
}

void C1Simplex::geometry()
{
	double a, x, y, z, u;
	bool n;
	CEdge *F;
	Sii = 0; Mu0i = 0; Mu1i = 0; At = 0; Conv = 0;
	if ((double)A->X == (double)B->X && (double)A->Y == (double)B->Y && (double)A->Z == (double)B->Z)
	{
		if ((double)A->U != (double)B->U) Sii = 1;
		else { x = A->X; y = A->Y; z = A->Z; u = x*x + y*y + z*z; Si = 0; }
	}
	else
	{
		x = (double)A->X - (double)B->X;
		y = (double)A->Y - (double)B->Y;
		z = (double)A->Z - (double)B->Z;

		a = (double)B->U - (double)A->U + 2 * ((double)A->X*x + (double)A->Y*y + (double)A->Z*z);
		a /= 2 * (x*x + y*y + z*z);

		x = (double)A->X - a*x;
		y = (double)A->Y - a*y;
		z = (double)A->Z - a*z;
		u = 2 * (x*(double)A->X + y*(double)A->Y + z*(double)A->Z) - (double)A->U;
		Si = x*x + y*y + z*z - u;
	}
	F = E; n = 1;
	while (1)
	{
		if (F->C->N<0) Conv = 1;
		else
		{
			if (u + (double)F->C->U <= 2 * (x*(double)F->C->X + y*(double)F->C->Y + z*(double)F->C->Z)) At = 1;
			if (n)
			{
				if (F->S2->At)
				{
					Mu0 = F->S2->Mu0; Mu0i = F->S2->Mu0i;
				}
				else
				{
					Mu0 = F->S2->Si; Mu0i = F->S2->Sii;
				}
				Mu1 = F->S2->Si; Mu1i = F->S2->Sii;
				//				if(!Mu1i && (F->S2->Mu1i || F->S2->Mu1>Mu1))
				//				{	Mu1=F->S2->Mu1; Mu1i=F->S2->Mu1i;}
				n = 0;
			}
			if (!F->S2->At && !F->S2->Sii && (Mu0i || F->S2->Si<Mu0)) { Mu0 = F->S2->Si; Mu0i = 0; }
			if (!F->S2->Mu0i && (Mu0i || F->S2->Mu0<Mu0)) { Mu0 = F->S2->Mu0; Mu0i = 0; }
			if (F->S2->Sii || F->S2->Mu1i) Mu1i = 1;
			else
			{
				if (F->S2->Si>Mu1) Mu1 = F->S2->Si;
				if (F->S2->Mu1>Mu1) Mu1 = F->S2->Mu1;
			}
		}
		F = F->Fnx; if (F == E) break;
	}
}

void C2Simplex::geometry()
{
	double x, x1, x2, x3, y, y1, y2, y3, z, z1, z2, z3, u, u1, u2, u3, Lx, Ly, Lz, Lu;
	CXYZW*D;
	E->check_minors();
	CInexact L = E->M023*E->M023 + E->M013*E->M013 + E->M012*E->M012;
	if (L.is_zero())
	{
		fprintf(stderr, "3 - degeneracy!\n"); return;
	}
	x1 = (double)A->X; x2 = (double)B->X; x3 = (double)C->X;
	y1 = (double)A->Y; y2 = (double)B->Y; y3 = (double)C->Y;
	z1 = (double)A->Z; z2 = (double)B->Z; z3 = (double)C->Z;
	u1 = (double)A->U; u2 = (double)B->U; u3 = (double)C->U;

	double m014 = u1*(x3 - x2) - u2*(x3 - x1) + u3*(x2 - x1);
	double m024 = u1*(y3 - y2) - u2*(y3 - y1) + u3*(y2 - y1);
	double m034 = u1*(z3 - z2) - u2*(z3 - z1) + u3*(z2 - z1);
	double m124 = u1*(x2*y3 - x3*y2) - u2*(x1*y3 - x3*y1) + u3*(x1*y2 - x2*y1);
	double m134 = u1*(x2*z3 - x3*z2) - u2*(x1*z3 - x3*z1) + u3*(x1*z2 - x2*z1);
	double m234 = u1*(y2*z3 - y3*z2) - u2*(y1*z3 - y3*z1) + u3*(y1*z2 - y2*z1);

	Lu = 2 * (double)E->M123*(double)E->M123 - (double)E->M012*m124 - (double)E->M013*m134 - (double)E->M023*m234;	//*4
	Lx = 2 * (double)E->M123*(double)E->M023 - (double)E->M013*m034 - (double)E->M012*m024;		//*2
	Ly = -2 * (double)E->M123*(double)E->M013 - (double)E->M023*m034 + (double)E->M012*m014;		//*2
	Lz = 2 * (double)E->M123*(double)E->M012 + (double)E->M023*m024 + (double)E->M013*m014;		//*2

	u1 = (double)L;
	Si = ((Lx*Lx + Ly*Ly + Lz*Lz) / 4 / u1 - Lu) / u1;
	u = Lu / u1; u1 *= 2; x = Lx / u1; y = Ly / u1; z = Lz / u1;
	Sii = 0; Mu0i = 0; Mu1i = 0;
	At = 0; Conv = 0;

	D = E->Fnx->C;
	if (D->N<0) { Conv = 1; Mu1i = 1; }
	else
	{
		if (u + (double)D->U <= 2 * (x*(double)D->X + y*(double)D->Y + z*(double)D->Z)) At = 1;
		if (E->S3->Sii) { Mu0i = 1; Mu1i = 1; }
		else { Mu0 = E->S3->Si; Mu1 = E->S3->Si; }
	}
	D = E->Rev->Fnx->C;
	if (D->N<0) { Conv = 1; Mu1i = 1; }
	else
	{
		if (u + (double)D->U <= 2 * (x*(double)D->X + y*(double)D->Y + z*(double)D->Z)) At = 1;
		if (E->Rev->S3->Sii) { Mu1i = 1; if (Conv) Mu0i = 1; }
		else
		{
			if (Conv || E->Rev->S3->Si<Mu0) Mu0 = E->Rev->S3->Si;
			if (Conv || E->Rev->S3->Si>Mu1) Mu1 = E->Rev->S3->Si;
		}
	}
}

void C3Simplex::geometry()
{
	double U1, U2, U3, U4, Lx, Ly, Lz, Lu;
	CEdge *abc, *abd, *acd, *bcd;
	abc = E;
	abd = E->Fnx;
	acd = E->Enx->Enx->Fnx->Rev;
	bcd = E->Enx->Fnx;
	abc->check_minors(); abd->check_minors();
	acd->check_minors(); bcd->check_minors();
	CInexact L = bcd->M123 - acd->M123 + abd->M123 - abc->M123;
	if (L.is_zero())
	{
		fprintf(stderr, "4 - degeneracy!\n"); return;
	}
	U1 = A->U; U2 = B->U; U3 = C->U; U4 = D->U;
	Lu = U4*(double)abc->M123 - U3*(double)abd->M123 + U2*(double)acd->M123 - U1*(double)bcd->M123;	//*8
	Lx = U4*(double)abc->M023 - U3*(double)abd->M023 + U2*(double)acd->M023 - U1*(double)bcd->M023;	//*4
	Ly = U1*(double)bcd->M013 - U2*(double)acd->M013 + U3*(double)abd->M013 - U4*(double)abc->M013;	//*4
	Lz = U4*(double)abc->M012 - U3*(double)abd->M012 + U2*(double)acd->M012 - U1*(double)bcd->M012;	//*4
	U1 = (double)L;
	Si = ((Lx*Lx + Ly*Ly + Lz*Lz) / 4 / U1 - Lu) / U1;
	Sii = 0;
	Mu0i = Sii; Mu1i = Sii;
	At = 0; Conv = 0;
}

void C1Simplex::dump()
{
	printf("size: 2\t%s;\t%s;\t", At ? "Attached" : "Unattached", Conv ? "On Conv" : "Not on Conv");
	if (Sii) printf("Sigma: Infinity\t");
	else printf("Sigma: %lg\t", Si);
	if (Mu0i) printf("Mu0: Infinity\t");
	else printf("Mu0: %lg\t", Mu0);
	if (Mu1i) printf("Mu1: Infinity\t");
	else printf("Mu1: %lg\t", Mu1);
	printf("%d %d", A->N, B->N);
	printf("\n");
}

void C2Simplex::dump()
{
	printf("size: 3\t%s;\t%s;\t", At ? "Attached" : "Unattached", Conv ? "On Conv" : "Not on Conv");
	if (Sii) printf("Sigma: Infinity\t");
	else printf("Sigma: %lg\t", Si);
	if (Mu0i) printf("Mu0: Infinity\t");
	else printf("Mu0: %lg\t", Mu0);
	if (Mu1i) printf("Mu1: Infinity\t");
	else printf("Mu1: %lg\t", Mu1);
	printf("%d %d %d", A->N, B->N, C->N);
	printf("\n");
}

void C3Simplex::dump()
{
	printf("size: 4\tUnattached;\tNot on Conv;\t");
	if (Sii) printf("Sigma: Infinity\t");
	else printf("Sigma: %lg\t", Si);
	printf("%d %d %d %d", A->N, B->N, C->N, D->N);
	printf("\n");
}

#endif

