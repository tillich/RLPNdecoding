//============================================================================
// Name        : DecodingCost.cpp
// Author      : Kévin Carrier
// Version     :
// Copyright   : open source
// Description : Optimization of various method for decoding linear code.
//               C++, Ansi-style
//============================================================================

#include "Tools.h"
#include "Prange.h"
#include "SternMayOzerov.h"
#include "BothMayOld.h"
#include "BothMay.h"
#include "MayOzerov.h"

using namespace std;
using namespace tools;


string algo = "Prange";
Decoder *isd;

// code parameters
double R = 0;
double t = 0;

// optimization parameters
bool optim_flag = true;
double radius = 0.032;
double step = 0.016;
double limit = 0.001;
double precision = 0.0000001;

// algorithm paramters
string path_params = "";


int parse_cmdline(int argc, char *argv[]){
	int i;

	for (i = 1 ; i < argc ; ++i){
		if (argv[i][0] == '-') {
			switch (argv[i][1]){
			case 'R':
				if (!argv[i][2]){
					++i;
					R = atof(argv[i]);
				}else{
					R = atof(argv[i]+2);
				}
				break;
			case 't':
				if (!argv[i][2]){
					++i;
					t = atof(argv[i]);
				}else{
					t = atof(argv[i]+2);
				}
				break;
			case 'A':
				if (!argv[i][2]){
					++i;
					algo = argv[i];
				}else{
					algo = argv[i]+2;
				}
				break;
			case 'O':
				if (!argv[i][2]){
					++i;
					radius = atof(argv[i]);
				}else{
					radius = atof(argv[i]+2);
				}
				++i;
				step = atof(argv[i]);
				++i;
				limit = atof(argv[i]);
				++i;
				precision = atof(argv[i]);
				break;
			case 'X':
				optim_flag = false;
				break;
			case 'P':
				if (!argv[i][2]){
					++i;
					path_params = argv[i];
				}else{
					path_params = argv[i]+2;
				}
				break;
			case 'h' :
				cout << "Usage : " << argv[0] << " <arguments>" << endl;
				cout << "\t\"-h\" : help" << endl;
				cout << "\t\"-R\" : Rate of the code to decode (double REQUIRED)" << endl;
				cout << "\t\"-t\" : Weight of the error vector to decode (double OPTIONAL : DEFAULT = Hinv(1 - R))" << endl;
				cout << "\t\"-A\" : The decoding algorithm to consider (string OPTIONAL : DEFAULT = \"Prange\")" << endl;
				cout << "\t\t Possible values : \"Prange\", \"SternMayOzerov\", \"BothMayOzerov\", \"BothMay\", \"BothMayOld\"" << endl;
				cout << "\t\"-O\" : The optimization parameters [radius, step, limit, precision] (double double double double OPTIONAL : DEFAULT = 0.032 0.016 0.001 0.0000001)" << endl;
				cout << "\t\"-X\" : To compute the complexity without optimization (void OPTIONAL)" << endl;
				cout << "\t\"-P\" : The path of the .txt file containing the initial parameters for the considered decoding algorithm (string OPTIONAL)" << endl;
				return 2;
				break;
			default :
				cerr << "Arguments parse error : argument \"-" << argv[i][1] << "\" unknown !" << endl;
				cerr << " -h for help !" << endl << endl;
				return 1;
			}
		}
	}
	if (!R){
		cerr << "Arguments parse error : " << endl;
		cerr << "\t argument -R is required !" << endl;
		cerr << " -h for help !" << endl << endl;
		return 1;
	}

	if (!t){
		t = Hinv(1.0 - R);
	}

	if (algo == "Prange") {
		isd = new Prange(R, t);
	} else if (algo == "SternMayOzerov") {
		isd = new SternMayOzerov(R, t);
	} else if (algo == "BothMay"){
		if (path_params == ""){
			isd = new BothMay(R, t);
		} else {
			double r[5] = {0,0,0,0,0};
			double p[5] = {0,0,0,0,0};
			double ww[5][5] = {{0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}};
			ifstream fd(path_params);
			string param_line;
			string param_name;
			double param_value;
			if(fd) {
				while(getline(fd, param_line)){
					if (param_line.find(';') < param_line.length()){
						param_name = param_line.substr(0,param_line.find_first_of("="));
						param_value = stod(param_line.substr(param_line.find_first_of("=") + 1, param_line.find_first_of(";")));
						if (param_name == "r[1]") { r[1] = param_value; }
						if (param_name == "r[2]") { r[2] = param_value; }
						if (param_name == "r[3]") { r[3] = param_value; }
						if (param_name == "r[4]") { r[4] = param_value; }
						if (param_name == "p[0]") { p[0] = param_value; }
						if (param_name == "p[1]") { p[1] = param_value; }
						if (param_name == "p[2]") { p[2] = param_value; }
						if (param_name == "p[3]") { p[3] = param_value; }
						if (param_name == "p[4]") { p[4] = param_value; }
						if (param_name == "ww[1][1]") { ww[1][1] = param_value; }
						if (param_name == "ww[1][2]") { ww[1][2] = param_value; }
						if (param_name == "ww[1][3]") { ww[1][3] = param_value; }
						if (param_name == "ww[1][4]") { ww[1][4] = param_value; }
						if (param_name == "ww[2][2]") { ww[2][2] = param_value; }
						if (param_name == "ww[2][3]") { ww[2][3] = param_value; }
						if (param_name == "ww[2][4]") { ww[2][4] = param_value; }
						if (param_name == "ww[3][3]") { ww[3][3] = param_value; }
						if (param_name == "ww[3][4]") { ww[3][4] = param_value; }
						if (param_name == "ww[4][4]") { ww[4][4] = param_value; }
					}
				}
			}
			else {
				cerr << "ERROR: Cannot open the file \'" << path_params << "\'" << endl;
				exit(EXIT_FAILURE);
			}
			fd.close();
			isd = new BothMay(R, t, r, p, ww);
		}
		((BothMay *)isd)->set_optim_params(radius, step, limit, precision);
	} else if (algo == "BothMayOld"){
		if (path_params == ""){
			isd = new BothMayOld(R, t);
		} else {
			double r[5] = {0,0,0,0,0};
			double p[5] = {0,0,0,0,0};
			double ww[5][5] = {{0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}, {0,0,0,0,0}};
			ifstream fd(path_params);
			string param_line;
			string param_name;
			double param_value;
			if(fd) {
				while(getline(fd, param_line)){
					if (param_line.find(';') < param_line.length()){
						param_name = param_line.substr(0,param_line.find_first_of("="));
						param_value = stod(param_line.substr(param_line.find_first_of("=") + 1, param_line.find_first_of(";")));
						if (param_name == "r[1]") { r[1] = param_value; }
						if (param_name == "r[2]") { r[2] = param_value; }
						if (param_name == "r[3]") { r[3] = param_value; }
						if (param_name == "r[4]") { r[4] = param_value; }
						if (param_name == "p[0]") { p[0] = param_value; }
						if (param_name == "p[1]") { p[1] = param_value; }
						if (param_name == "p[2]") { p[2] = param_value; }
						if (param_name == "p[3]") { p[3] = param_value; }
						if (param_name == "p[4]") { p[4] = param_value; }
						if (param_name == "ww[1][1]") { ww[1][1] = param_value; }
						if (param_name == "ww[1][2]") { ww[1][2] = param_value; }
						if (param_name == "ww[1][3]") { ww[1][3] = param_value; }
						if (param_name == "ww[1][4]") { ww[1][4] = param_value; }
						if (param_name == "ww[2][2]") { ww[2][2] = param_value; }
						if (param_name == "ww[2][3]") { ww[2][3] = param_value; }
						if (param_name == "ww[2][4]") { ww[2][4] = param_value; }
						if (param_name == "ww[3][3]") { ww[3][3] = param_value; }
						if (param_name == "ww[3][4]") { ww[3][4] = param_value; }
						if (param_name == "ww[4][4]") { ww[4][4] = param_value; }
					}
				}
			}
			else {
				cerr << "ERROR: Cannot open the file \'" << path_params << "\'" << endl;
				exit(EXIT_FAILURE);
			}
			fd.close();
			isd = new BothMayOld(R, t, r, p, ww);
		}
		((BothMayOld *)isd)->set_optim_params(radius, step, limit, precision);
	} else if (algo == "MayOzerov"){
		isd = new MayOzerov(R, t);
		((MayOzerov *)isd)->set_optim_param(precision);
	} else {
		cerr << "Arguments parse error : unknown decoding algorithm " << algo << " !" << endl;
		cerr << " -h for help !" << endl << endl;
		return 1;
	}

	return 0;
}



int main(int argc, char * argv[]){
	// double precision for printing
	cout.setf(ios::fixed);
	cout.precision(12);

	// parsing arguments
	switch(parse_cmdline(argc, argv)){
	case 0 :
		break;
	case 1 :
		exit(EXIT_FAILURE);
	case 2 :
		exit(EXIT_SUCCESS);
	}

	if (optim_flag){
		isd->optimization();
	}

	isd->set_verbose(true);
	isd->complexity();

	return 0;
}
