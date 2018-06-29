// added by Jieling, 05/02/2018
// to do delaunay triangulation

#ifndef DELAUNAY_H
#define DELAUNAY_H

#include <string>
#include <vector>
#include <math.h>
#include "../model/Elements/ModelElement.h"

using namespace std;

class CExact;
class CInexact;
typedef unsigned char* unsignedcharptr;
class CXYZW;
class CEdge;
class CCell;
class CD3DW;
typedef CEdge* EdgePtr;
class CSimplex;
class C0Simplex; // vertex
class C1Simplex; // edge
class C2Simplex; // triangle
class C3Simplex; // tetrahedra
class CComplex;

class CExact
{
	bool Plus;
	int Size, Exp;
	unsigned char *S;
private:
	void Init(long);
	void Init(double);
	void Add(const CExact&, const CExact&, int&Sz, int&Ex, unsignedcharptr&) const;
	void Sub(const CExact&, const CExact&, int&Sz, int&Ex, unsignedcharptr&) const;
	void Mul(const CExact&, const CExact&, int&Sz, int&Ex, unsignedcharptr&) const;
public:
	CExact() { S = 0; Size = 0; Exp = 0; Plus = true; }
	CExact(int n) { Init((long)n); }
	CExact(long n) { Init(n); }
	CExact(float f) { Init((double)f); }
	CExact(double f) { Init(f); }
	CExact(const CExact&);
	~CExact() { delete[] S; }
	const CExact& operator=(const CExact&);
	const CExact& operator=(int n) { return operator=((long)n); }
	const CExact& operator=(long n) { delete[]S; Init(n); return *this; }
	const CExact& operator=(float f) { return operator=((double)f); }
	const CExact& operator=(double f) { delete[]S; Init(f); return *this; }
	operator bool() const { return Size != 0; }
	operator double() const;
	operator float() const { return (float) operator double(); }
	operator long() const { return (long) operator double(); }
	operator int() const { return (int) operator double(); }
	char CompAbs(const CExact&) const;
	char Compare(const CExact&) const;
	bool operator==(const CExact&) const;
	bool operator!=(const CExact&) const;
	bool operator>(const CExact&) const;
	bool operator<(const CExact&) const;
	bool operator>=(const CExact&) const;
	bool operator<=(const CExact&) const;
	const CExact operator-() const;
	const CExact operator+(const CExact&) const;
	const CExact operator-(const CExact&) const;
	const CExact operator*(const CExact&) const;
	const CExact& operator+=(const CExact&);
	const CExact& operator-=(const CExact&);
	const CExact& operator*=(const CExact&);
	bool Positive() const { return Plus; }
	void Dump();
};

class CInexact
{
	double D;
	double E;
private:
	void Init(float);
	void Init(double);
	double Err();
public:
	CInexact() { D = 0; E = 0; }
	CInexact(int n) { D = n; E = 0; }
	CInexact(long n) { D = n; E = 0; }
	CInexact(float f) { Init(f); }
	CInexact(double f) { Init(f); }

	const CInexact& operator=(int n) { D = n; E = 0; }
	const CInexact& operator=(long n) { D = n; E = 0; }
	const CInexact& operator=(float f) { Init(f); return *this; }
	const CInexact& operator=(double f) { Init(f); return *this; }

	inline operator bool() const;
	operator double() const { return D; }
	inline bool is_zero() const;

	inline const CInexact operator-() const;
	inline const CInexact operator+(const CInexact&) const;
	inline const CInexact operator-(const CInexact&) const;
	inline const CInexact operator*(const CInexact&) const;
	inline const CInexact& operator+=(const CInexact&);
	inline const CInexact& operator-=(const CInexact&);
	inline const CInexact& operator*=(const CInexact&);
	void Dump();
};

class CXYZW
{
	protected:
		CInexact X, Y, Z, W, U;
		CExact eX, eY, eZ, eW, eU;
		int N;
		void *Label;
		C0Simplex *S0;
		ModelElement *p;

	public:
		CXYZW() { N = 0;S0 = NULL; };
		~CXYZW();
		int n() { return N; }
		double x() { return X; }
		double y() { return Y; }
		double z() { return Z; }
		double w() { return W; }
		ModelElement* P() { return p; } 
		void set_p(ModelElement* P) { p = P;  }
		void* label() const { return Label; }
		void set(double x, double y, double z, double w)
		{ 
			X = x;Y = y;Z = z;W = w;U = X*X + Y*Y + Z*Z - W; 
			eX= x;eY= y;eZ= z;eW= w;eU= eX*eX+eY*eY+eZ*eZ-eW;
		}
		void set(double x,double y,double z,double w,int n,void *l) { set(x,y,z,w);N=n;Label=l; }
		bool operator==(CXYZW &A) const { return (double)X == (double)A.X && (double)Y == (double)A.Y && (double)Z == (double)A.Z; }
	friend class CD3DW;
	friend class CEdge;
	friend class CCell;
	friend class CSimplex;
	friend class C0Simplex;
	friend class C1Simplex;
	friend class C2Simplex;
	friend class C3Simplex;
	friend class CComplex;
	friend bool edge_less(EdgePtr,EdgePtr);
};

class CEdge
{
	protected:
		CXYZW *A,*B,*C;
		CCell *Cell;
		CEdge *Rev,*Enx,*Fnx;
		CXYZW *LastV;
		bool LastP;
		char Flag; // 1:in stack;2:deleted
		int N;
		CInexact M123, M023, M013, M012;
		bool Minors;
		C1Simplex *S1;
		C2Simplex *S2;
		C3Simplex *S3;

	public:
		CEdge(CXYZW *a, CXYZW *b, CXYZW *c) { A = a;B = b;C = c;Flag = 0; LastV = NULL; Minors = false; }
		~CEdge();
		void setA(CXYZW* a) { A = a; }
		void setB(CXYZW* b) { B = b; }
		void setC(CXYZW* c) { C = c; }
		const CXYZW* a() const { return A; }
		const CXYZW* b() const { return B; }
		const CXYZW* c() const { return C; }
		inline void check_minors();
	friend class CD3DW;
	friend class CCell;
	friend class CSimplex;
	friend class C0Simplex;
	friend class C1Simplex;
	friend class C2Simplex;
	friend class C3Simplex;
	friend class CComplex;
	friend bool edge_less(EdgePtr,EdgePtr);
};

class CCell
{
	protected:
		CXYZW *VV[4];
		CEdge *abc,*adb,*acd,*bdc;
		CCell *CC[4];

	public:
		CCell() { CC[0] = 0;CC[1] = 0;CC[2] = 0;CC[3] = 0; }
		~CCell();
		CXYZW* get_VV(int i) {return VV[i];}
		CEdge* getabc() { return abc; }
		CEdge* getadb() { return adb; }
		CEdge* getacd() { return acd; }
		CEdge* getbdc() { return bdc; }
	friend class CD3DW;
	friend class CEdge;
	friend class CSimplex;
	friend class C0Simplex;
	friend class C1Simplex;
	friend class C2Simplex;
	friend class C3Simplex;
	friend class CComplex;
};

class CD3DW
{
	protected:
		std::vector<CEdge*> Es; // edge(triangle) list
		std::vector<CXYZW*> Vs; // vertex list
		std::vector<CCell*> Cs; // tetrahedra list
		std::vector<CEdge*> St; // edge(triangle) list for flip process
		CXYZW V0,V1,V2,V3;
		int N23, N32, N41, NLoc;
		int Ntp, Nrp, Ndd, Nep;
		int Nti, Nei;
		double L;
		
	public:
		CD3DW(double H=1.0e19); // the boundary vertices, removed after triangulation is done
		~CD3DW();
		inline CEdge* locate(CXYZW *V);                         // locate the triangle which contains the inserted point
		inline bool inside(CCell *C, CXYZW *V);                  // test if the point is inside one triangle
		inline void flip23(CEdge *abc);                         // flip the edge
		inline void flip32(CEdge *abc);
		inline void flip41(CEdge *abc);
		int size() { return (int)Vs.size(); }
		int nedges() { return (int)Es.size(); }
		int ncells() { return (int)Cs.size(); }
		int n23() { return N23; }
		int n32() { return N32; }
		int n41() { return N41; }
		int nloc() { return NLoc; }
		int nep() { return Nep; }
		int nei() { return Nei; }
		int ntp() { return Ntp; }
		int nti() { return Nti; }
		int nrp() { return Nrp; }
		int ndd() { return Ndd; }
		CXYZW* operator[](int n) {return Vs[n];}
		vector<CCell*> *get_Cs() {return &Cs;}
		vector<CXYZW*> *get_Vs() {return &Vs;}
		void add(double x, double y, double z, double w, int n, void *l);	// each time, randomly toss one point into the box
		inline bool positive(CEdge *E, CXYZW *V);
		inline bool insphere(CEdge *E,CXYZW *V);                // test if the point is inside the circumcircle of the triangle
		inline bool e_positive(CXYZW&A,CXYZW&B,CXYZW&C,CXYZW&D);
		inline bool e_insphere(CXYZW&A, CXYZW&B, CXYZW&C, CXYZW&D,CXYZW&E);
		inline bool sign4(CXYZW**P);
		inline bool sign5(CXYZW**P);
		inline void new6(EdgePtr&abc,EdgePtr&bca,EdgePtr&cab,
			EdgePtr&bac,EdgePtr&cba,EdgePtr&acb,CXYZW*v1,CXYZW*v2,CXYZW*v3);
		void list_edges(std::vector<EdgePtr>*);
		
		friend class CComplex;
};

class CSimplex
{
	protected:
		bool At; // attached
		double Si, Mu0, Mu1;
		double Sii, Mu0i, Mu1i;
		CEdge *E;
		virtual void geometry() = 0;
	public:
		~CSimplex() {};
		bool Conv; // on convex hull
		virtual CXYZW* operator[](int)=0;
		virtual int size() = 0;
		virtual bool attached() { return At; }
		virtual bool on_surface(double a = 0) { return false; }
		double x(int i) { return operator[](i)->X; }
		double y(int i) { return operator[](i)->Y; }
		double z(int i) { return operator[](i)->Z; }
		double w(int i) { return operator[](i)->W; }
		int n(int i) { return operator[](i)->N; }
		void* label(int i) {return operator[](i)->Label;}

	friend class CComplex;
	friend class C0Simplex;
	friend class C1Simplex;
	friend class C2Simplex;
};

class C0Simplex: public CSimplex
{
	public:
		virtual ~C0Simplex() {};
		CXYZW *A;
		int size() {return 1;}
		bool on_surface(double a);
		void geometry() {};
		CXYZW* operator[](int);
	friend class CComplex;
};

class C1Simplex: public CSimplex
{
	public:
		virtual ~C1Simplex() {};
		CXYZW *A,*B;
		int size() {return 2;}
		CXYZW* operator[](int);
		bool on_surface(double a);
		void dump();
		void geometry();       // alphavalue calculation
	friend class CComplex;
};

class C2Simplex: public CSimplex
{
	public:
		virtual ~C2Simplex() {};
		CXYZW *A,*B,*C;
		int size() {return 3;}
		CXYZW* operator[](int);
		bool on_surface(double a);
		void dump();
		void geometry();       // alphavalue calculation
	friend class CComplex;
};

class C3Simplex : public CSimplex
{
public:
	virtual ~C3Simplex() {};
	CXYZW *A, *B, *C, *D;
	int size() { return 4; }
	CXYZW* operator[](int);
	void dump();
	void geometry();
	friend class CComplex;
	friend class C2Simplex;
};

class CComplex
{
	public:
		CComplex(CD3DW*D);
		~CComplex();
		std::vector<CSimplex*> P;
};

#endif

