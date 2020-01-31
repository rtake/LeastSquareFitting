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

class Atom {
	public:
		Atom() {};
		Atom(string name, vector<double> config) { this->name = name; this->config = config; } // constructor
		Atom(const Atom &a) { name = a.name; config = a.config; } // copy constructor
		string getname() { return name; } // get atom name
		vector<double> getconfig() { return config; } // get configulation

		void SetAtomfromString(string line) {
			char buf[3]; vector<double> con(3,0);
			sscanf(line.c_str(),"%s%17lf%17lf%17lf",buf,&con[0],&con[1],&con[2]);
			string name(buf);
			this->name = name; this->config = con; // set
		}

	private:
		string name; // atom name
		vector<double> config; // atom configulation
};


class Config { // Config is a set of configulation of atoms and the energy
	public:
		Config() {}; // default constructor
		Config(vector<Atom> config, double energy) { this->config = config; this->energy = energy; } // constructor
		Config(const Config &c) { config = c.config; energy = c.energy; } // copy constructor
		double getenergy() { return energy; }
		vector<Atom> getconfig() { return config; }

		void AddAtom(Atom a) { config.push_back(a); }

	private:
		vector<Atom> config;
		double energy;
};


vector<Config> GetConfigfromlog(ifstream& log, string type) {
	vector<Config> vec_con; 
	string line; // int natom = -1;
	
	if(strstr(type.c_str(),"ITR") != NULL) { // get all structures of iterations from "# MIN" calc.
		while( getline(log,line) ) {
			if(strstr(line.c_str(),"# ITR.") != NULL) { // start loading
				Config c;
				while( getline(log,line) ) {
					Atom a; a.SetAtomfromString(line);
					if( strstr(line.c_str(),"Item") != NULL ) break;
					c.AddAtom(a);
				}
				vec_con.push_back(c);
			}
		}
	} else if(strstr(type.c_str(),"Opt") != NULL) { // get optimized structure from "# MIN" calc.

	} else { cout << "Error in GetConfigfromlog\n"; }

	return vec_con;
}

/*
int main(int argc, char* argv[]) {
	stringstream sslog; sslog << "h2o" << ".log";
	string type("ITR");
	vector<Config> vec_con;

	ifstream ifslog; ifslog.open( sslog.str().c_str() ); if(ifslog == NULL) { cout << "Not Found log file\n"; return -1; }
	vec_con = GetConfigfromlog(ifslog,type);
	ifslog.close();

	for(int i = 0;i < (int)vec_con.size();i++) {
		printf("ITR. %d\n",i);
		vector<Atom> vec_a = vec_con[i].getconfig();
		for(int j = 0;j < (int)vec_a.size();j++) {
			Atom a = vec_a[j]; 
			printf("%s%17lf%17lf%17lf\n",a.getname().c_str(),a.getconfig()[0],a.getconfig()[1],a.getconfig()[2]);
		}
	}

	return 0;
}
*/
