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

using namespace std;

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

	return 0;
}


