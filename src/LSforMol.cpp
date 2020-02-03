# include <iostream>
# include <string>
# include <cstring>
# include <fstream>
# include <stdlib.h>
# include <vector>
# include <cmath>
# include <set>
# include <sstream>
# include <utility>
# include <algorithm>

# include "LSPF.h"
# include "Rtlib.h"
# include <Eigen/Dense>

using namespace std;
using namespace Eigen;

/*
class DistFunc {
	public:
		DistFunc() {}
		DistFunc(int j, string type) : j(j), type(type) {}
		DistFunc(double r, int j, string type) : r(r), j(j), type(type) {}
		// DistFunc(const DistFunc& df) : j(df.j), type(df.type) {}
		DistFunc(const DistFunc& df) : r(df.r), j(df.j), type(df.type) {}
		void Set(int j, string type) { this->j = j; this->type = type; }
		void Set(double r, int j, string type) { this->r = r; this->j = j; this->type = type; }

		double func(double r) {
			if( strstr(type.c_str(),"inverse") ) return pow(r,-j);
			else if( strstr(type.c_str(),"exp") ) return pow( (1 - exp(-0.5 * r) ),j);
			else return -1;
		}
		double func() {
			if( strstr(type.c_str(),"inverse") ) return pow(r,-j);
			else if( strstr(type.c_str(),"exp") ) return pow( (1 - exp(-0.5 * r) ),j);
			else return -1;
		}

	private:
		double r; // distance
		int j; // degree
		string type; // function type
};


class BaseFunc {
	public:
		BaseFunc() {} // constructor
		// BaseFunc(const BaseFunc& bf) : vec_df(bf.vec_df), type(bf.type) {}
		BaseFunc(const BaseFunc& bf) : vec_df(bf.vec_df)  {}
		// void Set(vector<DistFunc> vec_df, string type) { this->vec_df = vec_df; this->type = type; }
		void Set(vector<DistFunc> vec_df) { this->vec_df = vec_df; } // set base function

		double func() {
			double v = 1;
			for(int i = 0;i < (int)vec_df.size();i++) v *= vec_df[i].func();
			return v;
		}

	private:
		vector<DistFunc> vec_df;
		// string type;
};


class FittingFunc {
	public:
		FittingFunc() {}
		FittingFunc(VectorXf a) : a(a) { cout << "Initialized with Coefficients VectorXf a :\n" << a << endl; }
		FittingFunc(const FittingFunc& ff) : a(ff.a), vec_bf(ff.vec_bf) {}
		void SetCoeff(VectorXf a) { this->a = a; }
		void SetBasissets(vector<BaseFunc> vec_bf) { this->vec_bf = vec_bf; }

		double func() {
			double v = 0; 
			for(int i = 0;i < (int)vec_bf.size();i++) v += a(i) * vec_bf[i].func();
			return v;
		} // get value

	private:
		VectorXf a; // Coefficients
		vector<BaseFunc> vec_bf; // basis sets
};


class LSFitting {
	public:
		LSFitting() {}
		LSFitting(int nrow, int nclm) : g(nrow,nclm), F(nrow) {}
		LSFitting(const LSFitting& lsf) : g(lsf.g), F(lsf.F) {}

		void SetValf(vector<double> vals) {
			for(int i = 0;i < (int)vals.size();i++) F(i) = vals[i]; 
			cout << "VectorXf F :\n" << F << endl;
		}

		void SetMatf(vector< vector<double> > mat) { 
			for(int i = 0;i < (int)mat.size();i++) { for(int j = 0;j < (int)mat[i].size();j++) g(i,j) = mat[i][j]; } 
			cout << "MatrixXf g :\n" << g << endl;
		}

		VectorXf fit() { return g.bdcSvd(ComputeThinU | ComputeThinV).solve(F); } // fitting here

	private:
		MatrixXf g;
		VectorXf F;
};


vector<BaseFunc> setBaseFuncvec(vector<double> distlist, string type, int rst) { // set vector of BaseFunc
	int npair = (int)distlist.size(); // distlist contains number of combination (number of atom pair)
	vector<BaseFunc> vec_bf;

	for(int i = 0;i < rst;i++) { // for each i ( < n) ; n means sumation restriction
		vector< vector<int> > mat_i;
		vector<int> vec_i(npair,0);

		vec_i[0] = i;
		do { mat_i.push_back(vec_i); } while( next_permutation( vec_i.begin(), vec_i.end() ) );

		int ncmb = (int)mat_i.size(); // number of combination (permutation)
		for(int j = 0;j < ncmb;j++) { // for each permunation (combination)
			BaseFunc bf;
			vector<DistFunc> vec_df; // set of DistFunc

			for(int k = 0;k < npair;k++) { // for each atom pair
				if(mat_i[j][k] > 0) { // if the degree is not zero
					DistFunc df(distlist[j],mat_i[j][k],type);
					vec_df.push_back(df); // add the base func
				} 
			}

			bf.Set(vec_df);
			vec_bf.push_back(bf);
		}

		for(int j = 0;j < ncmb;j++) {
			cout << "(";
			for(int k = 0;k < npair;k++) {
				cout << " " << mat_i[j][k];
			}
			cout << " )\n";
		}

	}

	cout << "set BaseFunc ok\n";
	return vec_bf;
}


void inPESdata(ifstream& ifs, vector< vector<double> >& mat_f, vector<double>& vec_f) {
	// vec_f[i]  \t  mat_f[i][0]  \t  mat_f[i][1]  \t  ...  \n

	string s;
	for(int i = 0;getline(ifs, s);i++) { cout << s << endl;
		stringstream ssline(s);
		string line;
		vector<double> vals; // vector of distance of each atom pair

		for(int j = 0;getline(ssline,line,'\t');j++) { cout << "line\t" << line << endl;
			double v;
			int chk = sscanf(line.c_str(),"%lf",&v);
			if(chk <= 0) printf("sscanf in inPESdata() failed, i : %d, j : %d\n",i,j);

			if(j == 0) vec_f[i] = v; // vec_f.push_back(v);			
			else if(j > 0) vals.push_back(v);
		}
	
		mat_f[i] = vals; // mat_f.push_back(vals);
	}
	printf("inPESdata() ok\n");
}


void outPESdata(ofstream& ofs, vector< vector<double> > mat_f, vector<double> vec_f) {
	// vec_f[i]  \t  mat_f[i][0]  \t  mat_f[i][1]  \t  ...  \n
	
	int n = (int)mat_f.size(); // vector of (vector of distance)
	
	if( n != (int)vec_f.size() ) { cout << "number of rows not matched !\n"; return; }
	else printf("natom chk --> ok (n : %d)\n",n);

	for(int i = 0;i < n;i++) {
		ofs << vec_f[i]; cout << vec_f[i];
		for(int j = 0;j < (int)mat_f[i].size();j++) { ofs << "\t" << mat_f[i][j]; cout << "\t" << mat_f[i][j]; }
		ofs << endl; cout << endl;
	}

	printf("outPESdata() ok\n");
}
*/

int main(int argc, char* argv[]) {

	// vals 

	int nref, nbase, natom, npoint;
	string refs, fit, type, line, all;
	ofstream ofs;
	ifstream ifsref, ifsall;

	
	// load input file xxx.in // ok

	if(argc <= 1) { cout << "No input file\n"; return -1; }

	string filename = argv[1];
	ifstream input( filename.c_str() );

	while( getline(input, line) ) {
		const char* pt = strstr(line.c_str(),":");
		char buf[256];
		if( strstr(line.c_str(),"nref") ) { sscanf(pt + 2,"%d",&nref); } // number of ref. point
		else if( strstr(line.c_str(),"nbase") ) { sscanf(pt + 2,"%d",&nbase); } // number of basis sets
		else if( strstr(line.c_str(),"natom") ) { sscanf(pt + 2,"%d",&natom); } // number of atoms
		else if( strstr(line.c_str(),"npoint") ) { sscanf(pt + 2,"%d",&npoint); } // number of points
		else if( strstr(line.c_str(),"refs") ) { sscanf(pt + 2,"%s", buf); refs = buf; ifsref.open( refs.c_str() ); }
		else if( strstr(line.c_str(),"fit") ) { sscanf(pt + 2,"%s", buf); fit = buf; ofs.open( fit.c_str() ); }
		else if( strstr(line.c_str(),"type") ) { sscanf(pt + 2,"%s", buf); type = buf; }
		else if( strstr(line.c_str(),"all") ) { sscanf(pt + 2,"%s", buf); all = buf; ifsall.open( all.c_str() ); }
	}


	// load data sets // ok

	vector< vector<double> > distref( nref, vector<double>( natom * (natom - 1) ) ), distlist( npoint, vector<double>( natom * (natom - 1) ) );
	vector<double> eneref(nref), enelist(npoint);
	inPESdata(ifsref,distref,eneref); // cout << "load ok\n";
	inPESdata(ifsall,distlist,enelist); // cout << "load ok\n"; outPESdata(ofs,distlist,enelist); // output here ok


	// set matrix for fitting

	vector< vector<double> > mat_f( nref, vector<double>(nbase) );
	for(int i = 0;i < nref;i++) { // for each ref. point
		vector<BaseFunc> vec_bf = setBaseFuncvec(distref[i],type,nbase); // set basis sets for this ref. point
		for(int j = 0;j < nbase;j++) mat_f[i][j] = vec_bf[j].func();
	}


	// fitting

	LSFitting lsf( (int)mat_f.size(), (int)mat_f[0].size() );
	lsf.SetMatf(mat_f); // set matrix
	lsf.SetValf(eneref); // set left vector value


	// calc. value from fitting function

	vector<FittingFunc> vec_ff( npoint, lsf.fit() ); // initializing with least square fitting coefficient
	for(int i = 0;i < npoint;i++) vec_ff[i].SetBasissets( setBaseFuncvec(distlist[i],type,nbase) );


	// file out

	vector<double> enefits(npoint);
	for(int i = 0;i < npoint;i++) enefits[i] = vec_ff[i].func();
	outPESdata(ofs,distlist,enefits);
	
	return 0;
}


