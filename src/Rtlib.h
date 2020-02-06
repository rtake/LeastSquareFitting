# ifndef INCLUDE_GUARD_RTLIB_H
# define INCLUDE_GUARD_RTLIB_H

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

using namespace std;

/*
	class Atom contain element and coordinate
	vector<Atom> match to Moleculuar
*/



class Atom {
	public:
		Atom() : crd(3,0) {}; // constructor
		Atom(const Atom &a) : elm(a.elm), crd(a.crd) {} // copy constructor
		vector<double> GetCrd() { return crd; }

		int SetfromString(string line) {
			char buf[3]; vector<double> con(3,0);
			int chk = sscanf(line.c_str(),"%s%17lf%17lf%17lf",buf,&crd[0],&crd[1],&crd[2]);

			if(chk == 4) { string name(buf); this->elm = name; return 0; }
			else return -1;
		}

	private:
		string elm; // element
		vector<double> crd; // coordinate
};


double Dist(Atom a0, Atom a1) {
	double sum = 0;
	vector<double> c0 = a0.GetCrd(), c1 = a1.GetCrd();

	for(int i = 0;i < 3;i++) sum += pow( (c0[i] - c1[i]), 2);
	return sqrt(sum);
}


void inPESdata(ifstream& ifs, vector< vector<double> >& mat_f, vector<double>& vec_f) {
	// vec_f[i]  \t  mat_f[i][0]  \t  mat_f[i][1]  \t  ...  \n

	string s;
	for(int i = 0;getline(ifs, s);i++) { // cout << s << endl;
		stringstream ssline(s);
		string line;
		vector<double> vals; // vector of distance of each atom pair

		for(int j = 0;getline(ssline,line,'\t');j++) { // cout << "line\t" << line << endl;
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
		for(int j = 0;j < (int)mat_f[i].size();j++) {
			ofs << "\t" << mat_f[i][j];
			// cout << "\t" << mat_f[i][j];
		}
		ofs << endl; cout << endl;
	}

	printf("outPESdata() ok\n");
}

# endif
