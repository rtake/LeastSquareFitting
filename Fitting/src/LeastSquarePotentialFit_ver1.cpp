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

# include <eigen-3.3.7/Eigen/Dense>
# include "Rtlib/Rtlib.h"

using namespace std;
using namespace Eigen;


vector< vector<int> > MakeCmb(int size, int sum) {
	vector< vector<int> > mat;

	if(size == 1) {
		for(int i = 0;i <= sum;i++) {
			vector<int> vec(1,i);
			mat.push_back(vec);
		}
	} else {
		for(int i = 0;i <= sum;i++) { // 
			vector< vector<int> > buf = MakeCmb(size - 1,i);
			for(int j = 0;j < buf.size();j++) { // for each vector
				vector<int> vec(1,sum - i);
				for(int k = 0;k < buf[j].size();k++) { vec.push_back(buf[j][k]); }
				mat.push_back(vec);
			}
		}
	}

	return mat;
}


MatrixXf SetMatrixXf(vector< vector<int> > mat_i, vector< vector<double> > datalist, string type) {
	int nbase = mat_i.size(); // number of basis functions
	int npoint = datalist.size(); // number of ref. point

	MatrixXf g(npoint,nbase);

	for(int i = 0;i < npoint;i++) { // for each ref. point, ...

		for(int j = 0;j < nbase;j++) { // for each base function, ...
			double val = 1;
			if( strstr(type.c_str(),"inverse") ) for(int k = 0;k < mat_i[j].size();k++) { val *= pow(datalist[i][k], -(mat_i[j][k]) ); }
			else if( strstr(type.c_str(),"exp") ) for(int k = 0;k < mat_i[j].size();k++) { val *= pow(1 - exp(-0.5 * datalist[i][k]) ,mat_i[j][k]); }
			g(i,j) = val;
		} // set function value

	} // set value

	return g; 
} // here calc. g_all / g_ref	// error


typedef struct ModelFunc_{
	VectorXf a; // coefficients
	vector< vector<int> > mat_i; // set of degrees of pow
	string btype; // type of basis sets
}ModelFunc;


ModelFunc LeastSquarePotentialFitting(ifstream& input) {

	ModelFunc mf;

	// vals 

	int i, j, nref, rst, natom, npoint, npair, nbase;
	string refs, fit, line, all, solver, log;
	ofstream ofs;
	ifstream ifsref, ifsall;

	FILE *fpfit, *fplog;

	
	// load input file xxx.in // ok

	while( getline(input, line) ) {
		const char* pt = strstr(line.c_str(),":");
		char buf[256];
		if( strstr(line.c_str(),"nref") ) { sscanf(pt + 2,"%d",&nref); } // number of ref. point
		else if( strstr(line.c_str(),"rst") ) { sscanf(pt + 2,"%d",&rst); } // number of basis sets
		else if( strstr(line.c_str(),"natom") ) { sscanf(pt + 2,"%d",&natom); } // number of atoms
		else if( strstr(line.c_str(),"npoint") ) { sscanf(pt + 2,"%d",&npoint); } // number of points
		else if( strstr(line.c_str(),"refs") ) { sscanf(pt + 2,"%s", buf); refs = buf; ifsref.open( refs.c_str() ); }
		else if( strstr(line.c_str(),"fit") ) { sscanf(pt + 2,"%s", buf); fit = buf; fpfit = fopen( fit.c_str(), "w" ); /* ofs.open( fit.c_str() ); */ }
		else if( strstr(line.c_str(),"type") ) { sscanf(pt + 2,"%s", buf); mf.btype = buf; }
		else if( strstr(line.c_str(),"all") ) { sscanf(pt + 2,"%s", buf); all = buf; ifsall.open( all.c_str() ); }
		else if( strstr(line.c_str(),"solver") ) { sscanf(pt + 2,"%s", buf); solver = buf; }
		else if( strstr(line.c_str(),"log") ) { sscanf(pt + 2,"%s", buf); log = buf; fplog = fopen( log.c_str(), "w" ); }
	}

	
	// vals

	npair = natom * (natom - 1) / 2;

	vector<double> eneref(nref);
	vector<double> enelist(npoint);
	vector< vector<double> > distref( nref, vector<double>(npair) );
	vector< vector<double> > distlist( npoint, vector<double>(npair) );

	VectorXf F(nref); // 
	MatrixXf g_ref; // matrix for ref.
	MatrixXf g_all; // matrix for all points

	// load data sets // ok

	inPESdata(ifsref,distref,eneref);
	fprintf( fplog, "1. load ref. points OK\n" );
	for(int i = 0;i < nref;i++) {
		fprintf( fplog, "# Ref. %d Energy : %17.12lf\n", i, eneref[i] );
		for(int j = 0;j < npair;j++) { fprintf( fplog, "dist[%d] : %17.12lf\n", j, distref[i][j] ); }
		fprintf( fplog, "\n" );
	}


	inPESdata(ifsall,distlist,enelist); 
	fprintf( fplog, "2. load all points OK\n" );

	for(int i = 0;i < nref;i++) F(i) = eneref[i]; // set VectorXf


	// set matrix for fitting

	fprintf( fplog, "3. SetMatrixXf for fitting START\n" );
	mf.mat_i = MakeCmb(npair,rst);
	nbase = mf.mat_i.size();

	g_ref = SetMatrixXf(mf.mat_i, distref, mf.btype); 
	fprintf( fplog, "SetMatrixXf for fitting OK\n" );


	// fitting

	if( strstr(solver.c_str(),"svd") ) { mf.a = g_ref.bdcSvd(ComputeThinU | ComputeThinV).solve(F); }
	if( strstr(solver.c_str(),"normal") ) { mf.a = (g_ref.transpose() * g_ref).ldlt().solve(g_ref.transpose() * F); }
	else { mf.a = g_ref.colPivHouseholderQr().solve(F); }
	fprintf( fplog, "4. Fitting OK\n" );


	// calc. value from fitting function

	fprintf( fplog, "5. SetMatrixXf for calc. value START\n" );
	g_all = SetMatrixXf(mf.mat_i, distlist, mf.btype);	// g_all(i, j) is calculated here	// error

	for(int i = 0;i < npoint;i++) { // for each data point
		enelist[i] = 0;
		fprintf( stdout, "# Point. %d\n", i );
		for(int j = 0;j < nbase;j++) {
			enelist[i] += g_all(i,j) * mf.a(j);	// energy calculated here	// error in g_all(i, j)
			fprintf( stdout, "coeff[%d] : %17.12lf, BaseFunc : %17.12lf\n", j, mf.a(j), g_all(i,j) );
		}
	}


	// file out

	fprintf( fpfit, "Coeff.\n" );
	for(i = 0;i < nbase;i++) { fprintf( fpfit, "%17.12lf\n", mf.a(i) ); }

	fprintf( fpfit, "Fit\n");
	for(i = 0;i < npoint;i++) {
		fprintf( fpfit, "# Point. %d : Energy %17.12lf\n", i, enelist[i] );
		for(j = 0;j < npair;j++) {
			fprintf( fpfit, "dist[%d] : %17.12lf\n", j, distlist[i][j] );
		}
		fprintf( fpfit, "\n" );
	}


	// 1. energy, configulation
	// 2. coeff, basis sets

	fprintf( fplog, "-- Fitting data --\n" );
	fprintf( fplog, "Energy\tData\t\n" );
	for(int i = 0;i < npoint;i++) {
		fprintf( fplog, "# Point %d Energy : %17.12lf\n", i, enelist[i] );	// energy chk	// error
		for(int j = 0;j < npair;j++) { fprintf( fplog, "dist[%d] : %17.12lf\n", j, distlist[i][j] ); }	// data chk	// ok
		fprintf( fplog, "\n" );
	}

	fprintf( fplog, "-- Fitting function --\n");
	fprintf( fplog, "Coeff.\tBasis set\t\n");
	for(int i = 0;i < npoint;i++) { // for each point
		fprintf( fplog, "# Point. %d\n",i);
		for(int j = 0;j < nbase;j++) {
			fprintf( fplog, "%17.12lf(%17.12lf)\n", mf.a(j), g_all( i, j ) );
			for(int k = 0;k < npair;k++) { fprintf( fplog, "dist[%d] : %17.12lf(%d)\n", k, distlist[i][k], mf.mat_i[j][k] ); }
			fprintf( fplog, "\n" );
		}
	}

	return mf;
}


int main(int argc, char* argv[]) {
	ifstream input(argv[1]);
	LeastSquarePotentialFitting(input);

	return 0;
}



