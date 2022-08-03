
struct informationIteration{
    size_t weight_on_U;
    size_t weight_on_N;
    size_t number_distinct_parityChecks;
	double theoric_error;
    double average_proba_error;


	size_t position_solution;
	long int walsh_value_solution;

	long int first_value_walsh;
	long int second_value_walsh;
    string value_print_python(){
        string str = "[" +to_string(weight_on_U) + ","+to_string(weight_on_N) + "," +to_string(number_distinct_parityChecks) + "," + to_string(theoric_error) + ","+ to_string(average_proba_error) + "," + to_string(position_solution) + "," + to_string(walsh_value_solution) + "," + to_string(first_value_walsh) + ","+to_string(second_value_walsh) + "] \n";
        return str;
    }
    void print_python(ofstream file){
        string str = value_print_python();
        file << str;
    }

};
struct informationDecode{
    long int treshhold;
    string value_print_python(){
        string str = "[";
       
        for(int i = 0; i < iterations.size(); i++)
        {
            str += iterations[i].value_print_python();
            if(i != iterations.size() - 1){
                str += ",";
            }
        }
        str+="]\n \n";
        return str;
    }
    std::vector<informationIteration> iterations;
};
struct informationAverage{
    long int treshhold;
    string value_print_python(){
        string str = "[";
       
        for(int i = 0; i <codes.size(); i++)
        {
            str += codes[i].value_print_python();
            if(i != codes.size() - 1){
                str += ",";
            }
        }
        str+="]";
        return str;
    }
    std::vector<informationDecode> codes;
};