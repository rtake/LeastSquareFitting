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


class Distfunc {
	public:
		Distfunc() {}
		Distfunc(double r, int j) : r(r), j(j) {}
		Distfunc(const Distfunc& df) : r(df.r), j(df.j) {}
		void Set(double r, int j) { this->r = r; this->j = j; }

		double func(string type) {
			if( strstr(type.c_str(),"inverse") ) return pow(r,-j);
			else if( strstr(type.c_str(),"exp") ) return pow( (1 - exp(-0.5 * r) ),j);
			else return -1;
		}

	private:
		double r; // distance
		int j; // degree
};


/*
class InverseDist {
	public:
		InverseDist() {} // constructor
		InverseDist(int deg, double dist) : deg(deg), dist(dist) {} // constructor
		InverseDist(const InverseDist &i) : deg(i.deg), dist(i.dist) {} // copy constructor
		void Set(int deg, double dist) { this->deg = deg; this->dist = dist; } // Setter
		double func() { return pow(dist, -deg); }
	private:
		int deg; // degree
		double dist; // distance
};


class Exptype {
	public:	
		Exptype() {}
		Exptype(int deg, double dist) : deg(deg), dist(dist) {} // constructor
		Exptype(const Exptype &e) : deg(e.deg), dist(e.dist) {} // copy constructor
		void Set(int deg, double dist) { this->deg = deg; this->dist = dist; } // Setter
		double func() { return pow( (1 - exp(-0.5 * dist) ), deg); }
	private:
		int deg;
		double dist;
};
*/


class BaseFunc {
	public:
		BaseFunc() {} // constructor
		// BaseFunc(int deg, double dist) : deg(deg), dist(dist) {} // constructor
		BaseFunc(const BaseFunc &bf) : vec_df(bf.vec_df) {} // copy constructor
		void Set(vector<Distfunc> vec_df) { this->vec_df = vec_df; } // Setter

		double func(string type) {
			double v = 1;
			for(int i = 0;i < (int)vec_df.size();i++) v *= vec_df[i].func(type);
			return v;
		}

	private:
		vector<Distfunc> vec_df;
};


class Fittingfunc {
	public:
		Fittingfunc() {}
		Fittingfunc(int nref, int nbase) : g(nref,nbase), F(nref) {}
		Fittingfunc(const Fittingfunc& ff) : g(ff.g), F(ff.F), a(ff.a) {}
		void SetValf(vector<double> vals) { for(int i = 0;i < (int)vals.size();i++) { F(i) = vals[i]; } }

		void SetMatf(vector< vector<double> > mat) {
			for(int i = 0;i < (int)mat.size();i++) {
				for(int j = 0;j < (int)mat[i].size();j++) g(i,j) = mat[i][j];
			}
		}

		void fit() { a = g.bdcSvd(ComputeThinU | ComputeThinV).solve(F); }

		double GetVal(vector<double> vec) {
			double v = 0;
			for(int i = 0;i < (int)vec.size();i++) v += a(i) * vec[i];
			return v;
		}

	private:
		MatrixXf g; // 
		VectorXf F; // real value
		VectorXf a; // coefficients	
};


int main(int argc, char* argv[]) {
	int nref, nbase;
	string refs, fit, type, filename, line;
	ofstream ofs;
	ifstream ifs;

	if(argc <= 1) { cout << "No input file\n"; return -1; }

	filename = argv[1];
	ifstream input( filename.c_str() );

	int cnt = 0;
	while( getline(input, line) ) { cout << cnt++ << " " << line << endl;
		const char* pt = strstr(line.c_str(),":");
		char buf[256];
		if( strstr(line.c_str(),"nref") ) { sscanf(pt + 2,"%d",&nref); /* cout << nref << endl; */ }
		else if( strstr(line.c_str(),"nbase") ) { sscanf(pt + 2,"%d",&nbase); /* cout << nbase << endl; */ }
		else if( strstr(line.c_str(),"refs") ) { sscanf(pt + 2,"%s", buf); refs = buf; /* cout << refs << endl;*/  ifs.open( refs.c_str() ); }
		else if( strstr(line.c_str(),"fit") ) { sscanf(pt + 2,"%s", buf); fit = buf; /* cout << fit << endl; */  ofs.open( fit.c_str() ); }
		else if( strstr(line.c_str(),"type") ) { sscanf(pt + 2,"%s", buf); type = buf; /* cout << type << endl; */ }
	}

	Fittingfunc ff(nref,nbase);

	double distmax = 3.0, distmin = 0.3, step = 0.005;
	int npoint = (distmax - distmin) / step;

	vector<double> distref(nref),eneref(nref);
	for(int i = 0;i < nref;i++) {
		getline(ifs,line);
		sscanf(line.c_str(),"%lf\t%lf\n", &distref[i], &eneref[i]);
	}

	ff.SetValf(eneref);

	vector< vector<Distfunc> > mat_df( nref, vector<Distfunc>(nbase) );
	vector< vector<double> > mat_f( nref, vector<double>(nbase) );
	for(int i = 0;i < nref;i++) {
		for(int j = 0;j < nbase;j++) {
			mat_df[i][j].Set(distref[i],j + 1);
			mat_f[i][j] = mat_df[i][j].func(type);
		}	
	}

	ff.SetMatf(mat_f);
	ff.fit(); // fitting

	vector<double> vals(npoint), distlist(npoint);
	for(int i = 0;i < npoint;i++) {
		distlist[i] = distmin + step * i;

		vector<double> argvec(nbase,0);
		for(int j = 0;j < nbase;j++) {
			Distfunc df(distlist[i],j + 1);
			argvec[j] = df.func(type);
		}

		vals[i] = ff.GetVal(argvec);
	}

	for(int i = 0;i < npoint;i++) ofs << distlist[i] << "\t" << vals[i] << endl;


/*
	if(strstr(type.c_str(),"inverse") != NULL) { cout << "InverseDist type" << endl;
		vector< vector<InverseDist> > imat( nref, vector<InverseDist>(nbase) );
		for(int i = 0;i < nref;i++) {
			for(int j = 0;j < nbase;j++) {
				imat[i][j].Set(j + 1,distlist[i]);
				g(i,j) = imat[i][j].func();
			}
		}

		a = g.bdcSvd(ComputeThinU | ComputeThinV).solve(F); // fitting here

		for(int i = 0;i < npoint;i++) {
			double energy = 0, distance = distmin + step * i;
			for(int j = 0;j < nbase;j++) { energy += a(j) * pow(distance, - (j + 1) ); }
			ofs << distance << "\t" << energy << endl;
		}

	} else if(strstr(type.c_str(),"exp") != NULL) { cout << "Exponential type" << endl;
		vector< vector<Exptype> > emat( nref, vector<Exptype>(nbase + 1) );
		for(int i = 0;i < nref;i++) {
			for(int j = 0;j < nbase + 1;j++) {
				emat[i][j].Set(j,distlist[i]);
				g(i,j) = emat[i][j].func();
			}
		}

		a = g.bdcSvd(ComputeThinU | ComputeThinV).solve(F); // fitting here

		for(int i = 0;i < npoint;i++) {
			double energy = 0, distance = distmin + step * i;
			for(int j = 0;j < nbase + 1;j++) { energy += a(j) * pow( (1 - exp(-0.5 * distance) ), j); }
			ofs << distance << "\t" << energy << endl;
		}

	} else cout << "No type specified\n";		

	cout << "Here is the matrix g:\n" << g << endl;
	cout << "Here is the right hand side F:\n" << F << endl;
	cout << "The least-squares solution is:\n" << a << endl;
*/
/*	for(int i = 0;i < nref;i++) {
		double F = 0;
		for(int j = 0;j < nbase + 1;j++) { F += a(j) * g(i,j); }
		printf("F %lf\n",F);
	}*/

	return 0;
}


