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


/*
class BaseFunc {
	public:
		BaseFunc() {} // constructor
		BaseFunc(int deg, double dist) : deg(deg), dist(dist) {} // constructor
		BaseFunc(const BaseFunc &i) : deg(i.deg), dist(i.dist) {} // copy constructor
		void Set(int deg, double dist) { this->deg = deg; this->dist = dist; } // Setter
		double func() {}
	private:
		int deg; // degree
		double dist; // distance
};


class InverseDist : BaseFunc {
	public:
		double func() { return pow(dist, -deg); }
};


class Exptype : BaseFunc {
	public:
		double func() { return pow( (1 - exp(-0.5 * dist) ), deg); }
};
*/

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

	MatrixXf g(nref,nbase);
	VectorXf F(nref), a;
	double distmax = 3.0, distmin = 0.3, step = 0.005;
	int npoint = (distmax - distmin) / step;

	vector<double> distlist(nref);
	for(int i = 0;i < nref;i++) {
		string line;
		double distance,energy;
		
		getline(ifs,line);
		sscanf(line.c_str(),"%lf\t%lf\n", &distance, &energy);

		F(i) = energy; distlist[i] = distance;
	}

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

/*	for(int i = 0;i < nref;i++) {
		double F = 0;
		for(int j = 0;j < nbase + 1;j++) { F += a(j) * g(i,j); }
		printf("F %lf\n",F);
	}*/

	return 0;
}


