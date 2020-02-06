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

# include "Rtlib.h"

using namespace std;

/*
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

int GetfromMINlog(ifstream& ifs, vector< vector<Atom> >& mols, vector<double>& elist) {
	int nmol = 0;
	string line;

	while( getline(ifs,line) ) { // cout << line << endl;
		
		if( strstr( line.c_str(), "# ITR.") ) { // cout << line << endl;

			vector<Atom> mol;
			while( getline(ifs,line) ) { // get molecule
				Atom a;
				if(a.SetfromString(line) == 0) {
					mol.push_back(a);
					// cout << line << endl; 
				} else break;
			}

			while( getline(ifs,line) ) { // get energy
				double e;
				const char* pt = strstr( line.c_str(), "ENERGY");
				if(pt) { // cout << line << endl; 
					sscanf(pt + 11,"%17lf",&e); elist.push_back(e); 
					break; 
				}
			}

			mols.push_back(mol);
			nmol++;
		}

	}

	return nmol;
}


int GetfromMINlog(ifstream& ifs, vector< vector<Atom> >& mols, vector<double>& elist) {
	int nmol = 0;
	string line;

	while( getline(ifs,line) ) {




	}

	return nmol;
}


int main(int argc, char* argv[]) {

	// vals

	ifstream ifs; // log file
	ofstream ofs; // data file
	vector< vector<double> > dlist;
	vector<double> elist;
	vector< vector<Atom> > mols;
	int nmol, natom;

	// args analysis

	if(argc < 5) { cout << "No input file (-i input -o output)\n"; return -1; }
	else if(argc > 5) { cout << "Too much args (-i input -o output)\n"; return -1; }

	for(int i = 0;i < argc;i++) {
		if(argv[i][0] == '-') {
			if(argv[i][1] == 'i') ifs.open(argv[i + 1]);
			else if(argv[i][1] == 'o') ofs.open(argv[i + 1]);
		}
	}


	// loading input xxx.log

	nmol = GetfromMINlog(ifs,mols,elist); // number of moleculuar
	natom = (int)mols[0].size(); // number of atom

	for(int i = 0;i < nmol;i++) { // for each molecule

		vector<double> vec_d; // vector of double
		for(int j = 0;j < natom;j++) {
			for(int k = j + 1;k < natom;k++) {
				vec_d.push_back( Dist(mols[i][j],mols[i][k]) );
			}
		}

		dlist.push_back(vec_d);
	}


	// output data

	outPESdata(ofs,dlist,elist);

	return 0;
}


