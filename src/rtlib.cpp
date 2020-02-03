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

# include "Rtlib.h"

using namespace std;

/*
	class Atom contain element and coordinate
	vector<Atom> match to Moleculuar
*/


/*
class Atom {
	public:
		Atom() : crd(3,0) {}; // constructor
		Atom(const Atom &a) : elm(a.elm), crd(a.crd) {} // copy constructor
		vector<double> GetCrd() { return crd; }

		void SetfromString(string line) {
			char buf[3]; vector<double> con(3,0);
			sscanf(line.c_str(),"%s%17lf%17lf%17lf",buf,&crd[0],&crd[1],&crd[2]);
			string name(buf); this->elm = name;
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
*/

// int main(int argc, char* argv[]) { return 0; }
