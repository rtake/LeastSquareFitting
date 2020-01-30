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

# include <Eigen/Dense>

using namespace std;
using namespace Eigen;


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
		FittingFunc(VectorXf a) : a(a) {}
		FittingFunc(const FittingFunc& ff) : a(ff.a), vec_b(ff.vec_b) {}
		void SetCoeff(VectorXf a) { this->a = a; }
		void SetBasissets(vector<BaseFunc> vec_bf) { this->vec_b = vec_b; }

		double func() {
			double v = 0; 
			for(int i = 0;i < (int)vec_b.size();i++) v += a(i) * vec_b[i].func();
			return v;
		} // get value

	private:
		VectorXf a; // Coefficients
		vector<BaseFunc> vec_b; // basis sets
};


class LSFitting {
	public:
		LSFitting() {}
		LSFitting(const LSFitting& lsf) : g(lsf.g), F(lsf.F) {}
		void SetValf(vector<double> vals) { for(int i = 0;i < (int)vals.size();i++) F(i) = vals[i]; }

		void SetMatf(vector< vector<double> > mat) { 
			for(int i = 0;i < (int)mat.size();i++) { for(int j = 0;j < (int)mat[i].size();j++) g(i,j) = mat[i][j]; } 
		}

		VectorXf fit() { return g.bdcSvd(ComputeThinU | ComputeThinV).solve(F); } // fitting here

	private:
		MatrixXf g;
		VectorXf F;
};


vector<BaseFunc> setBaseFuncvec(vector<double> distlist, string type) { // set vector of BaseFunc
	vector<BaseFunc> vec_b;
	return vec_b;
}


void inPESdata(const ifstream& ifs, vector< vector<double> > mat_f, vector<double> vec_f) {}

void outPESdata(const ofstream& ofs, vector< vector<double> > mat_f, vector<double> vec_f) {}


int main(int argc, char* argv[]) {
	int nref, nbase, natom, npoint;
	string refs, fit, type, line, all;
	ofstream ofs;
	ifstream ifsref, ifsall;

	if(argc <= 1) { cout << "No input file\n"; return -1; }

	string filename = argv[1];
	ifstream input( filename.c_str() );

	while( getline(input, line) ) { cout << cnt++ << " " << line << endl;
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


	// load data sets

	vector< vector<double> > distref( nref, vector<double>( natom * (natom - 1) ) ), distlist( npoint, vector<double>( natom * (natom - 1) ) );
	vector<double> eneref(nref), enelist(npoint);
	inPESdata(ifsref,distref,eneref); inPESdata(ifsall,distlist,enelist);


	// set matrix for fitting

	vector< vector<double> > mat_f( nref, vector<double>(nbase) );
	for(int i = 0;i < nref;i++) { // for each ref. point
		vector<BaseFunc> vec_b = setBaseFuncvec(distref[i],type); // set basis sets for this ref. point
		for(int j = 0;j < nbase;j++) mat_f[i][j] = vec_b[j].func();
	}


	// fitting

	LSFitting lsf;
	lsf.SetMatf(mat_f); // set matrix
	lsf.SetValf(eneref); // set lef vector value


	// calc. value from fitting function

	vector<FittingFunc> vec_ff( npoint, lsf.fit() ); // initializing with least square fitting coefficient
	for(int i = 0;i < npoint;i++) vec_ff[i].SetBasissets( setBaseFuncvec(distlist[i],type) );


	// file out

	vector<double> enefits(npoint);
	for(int i = 0;i < npoint;i++) enefits[i] = vec_ff[i].func();
	outPESdata(ofs,distlist,enefits);
	
	return 0;
}


