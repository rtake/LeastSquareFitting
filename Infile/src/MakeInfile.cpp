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

int GetfromMINlog(ifstream& ifs, vector< vector<Atom> >& mols, vector<double>& elist) {
	int nmol = 0;
	double e;
	char *pt;
	string line;

	while( getline(ifs,line) ) {
		
		if( strstr( line.c_str(), "# ITR.") ) { 

			vector<Atom> mol;
			while( getline(ifs,line) ) { 
				Atom a;
				if(a.SetfromString(line) == 0) {
					mol.push_back(a);
				} else break;
			}

			while( getline(ifs,line) ) { 
				pt = strstr( line.c_str(), "ENERGY");
				if(pt) {
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


int GetfromLUPlog(ifstream& ifs, vector< vector<Atom> >& mols, vector<double>& elist) {
        int nmol = 0;
        string line;

        while( getline(ifs,line) ) {

                if( strstr( line.c_str(), "# NODE") ) { // printf("%s\n",line.c_str());

                        vector<Atom> mol;
                        while( getline(ifs,line) ) {
                                Atom a;
                                if(a.SetfromString(line) == 0) {
                                        mol.push_back(a);
                                } else break;
                        }

                        while( getline(ifs,line) ) {
                                double e;
                                const char* pt = strstr( line.c_str(), "ENERGY");
                                if(pt) {
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


int GetfromIRClog(ifstream& ifs, vector< vector<Atom> >& mols, vector<double>& elist) {
	int nmol = 0;
	string line;

	vector<Atom> ts;
	vector< vector<Atom> > fwd,bck;

	double e_ts;
	vector<double> elist_fwd,elist_bck; 

	while( getline(ifs,line) ) {

		if( strstr(line.c_str(),"INITIAL STRUCTURE") ) { // TS

			while( getline(ifs,line) ) {
				Atom a;
				const char* pt = strstr( line.c_str(), "ENERGY");

				if(a.SetfromString(line) == 0) { ts.push_back(a); }
				if(pt) {
					sscanf(pt + 11,"%17lf",&e_ts);
					break;				
				}				

			}

			nmol++;

		} 
		
		if( strstr(line.c_str(),"IRC FOLLOWING (FORWARD)") ) { cout << "FWD start\n";

			while( getline(ifs,line) ) {
				if( strstr(line.c_str(),"IRC FOLLOWING (BACKWARD)") ) break;

				if( strstr(line.c_str(),"# STEP") ) { 
					vector<Atom> mol; // buffer	

					while( getline(ifs,line) ) {
						Atom a;
						double e;
						const char* pt = strstr( line.c_str(), "ENERGY");

						if(a.SetfromString(line) == 0) { mol.push_back(a); }
						if(pt) { // cout << line << endl;
							sscanf(pt + 11,"%17lf",&e);
							elist_fwd.push_back(e);
							break;
						}
					}

					fwd.push_back(mol);
					nmol++;
				}

			}

		}
		
		if( strstr(line.c_str(),"IRC FOLLOWING (BACKWARD)") ) { cout << "BCK start\n";

			while( getline(ifs,line) ) { // cout << line << endl;
				if( strstr(line.c_str(),"Energy profile along IRC") ) break;

				if( strstr(line.c_str(),"# STEP") ) {
					vector<Atom> mol;

					while( getline(ifs,line) ) {
						Atom a;
						double e;
						const char* pt = strstr( line.c_str(), "ENERGY");

						if(a.SetfromString(line) == 0) { mol.push_back(a); }
						if(pt) { // cout << line << endl;
							sscanf(pt + 11,"%17lf",&e);
							elist_bck.push_back(e);
							break;
						}
					}

					bck.push_back(mol);
					nmol++;
				}

			}

		}

	}

	for(int i = bck.size() - 1;i >= 0;i--) {
		mols.push_back(bck[i]);
		elist.push_back(elist_bck[i]); // cout << elist_bck[i] << endl;
	}

	mols.push_back(ts);
	elist.push_back(e_ts);

	for(int i = 0;i < fwd.size();i++) {
		mols.push_back(fwd[i]);
		elist.push_back(elist_fwd[i]); // cout << elist_fwd[i] << endl;
	}

	cout << "GetfromIRClog() ok\n";
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
	string logtype;

	// args analysis

	if(argc < 7) { cout << "No input file (-i input -o output -t logtype(min or irc))\n"; return -1; }
	else if(argc > 7) { cout << "Too much args (-i input -o output -t logtype(min or irc))\n"; return -1; }

	for(int i = 0;i < argc;i++) {
		if(argv[i][0] == '-') {
			if(argv[i][1] == 'i') ifs.open(argv[i + 1]);
			else if(argv[i][1] == 'o') ofs.open(argv[i + 1]);
			else if(argv[i][1] == 't') {
				if( strstr(argv[i + 1],"min") ) logtype = "min";
				else if( strstr(argv[i + 1],"irc") ) logtype = "irc";
				else if( strstr(argv[i + 1],"lup") ) logtype = "lup"; // xxx_LUPOUTt.log
			}
		}
	}


	// loading input xxx.log

	if(logtype == "min") { nmol = GetfromMINlog(ifs,mols,elist); }
	else if(logtype == "irc") { nmol = GetfromIRClog(ifs,mols,elist); }
	else if(logtype == "lup") { nmol = GetfromLUPlog(ifs,mols,elist); }
	else cout << "can't find log file\n";

	natom = (int)mols[0].size();
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


