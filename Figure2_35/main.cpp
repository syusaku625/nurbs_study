#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<sstream>
#include<cmath>

using namespace std;

double return_N0(int index, double c, vector<double> knot_vector)
{
    if(knot_vector[index] > c && c>=knot_vector[index-1.0]){
        return 1.0;
    }
    else{
        return 0.0;
    }
}

double p_order_N(vector<double> knot_vector, int index, int order, double c)
{
    if(order==0){
        double N = return_N0(index, c, knot_vector);
        return N;
    }
    double N1, N2;
    if(knot_vector[index-1+order]-knot_vector[index-1]==0){
        N1 = 0.0;
    }
    else{
        N1 = (c-knot_vector[index-1])/(knot_vector[index-1+order]-knot_vector[index-1]);
    }
    if(knot_vector[index-1+order+1]-knot_vector[index-1+1]==0){
        N2 = 0.0;
    }
    else{
        N2 = (knot_vector[index-1+order+1]-c)/(knot_vector[index-1+order+1]-knot_vector[index-1+1]);
    }
    double N = N1 * p_order_N(knot_vector, index, order-1, c) + N2 * p_order_N(knot_vector, index+1, order-1, c);
    return N;
}

int main()
{
    int order = 2;
    ifstream ifs("control_point.dat");
    string str;
    vector<vector<double>> control_point;
    while(getline(ifs, str)){
        istringstream ss(str);
        vector<double> tmp_point;
        for(int i=0; i<2; i++){
            getline(ss, str, ' ');
            tmp_point.push_back(stod(str));
        }
        control_point.push_back(tmp_point);
    }
    ifs.close();

    ifs.open("weight.dat");
    vector<double> weight;
    while(getline(ifs,str)){
        if(stoi(str)==1){
            weight.push_back(1.0);
        }
        else{
            weight.push_back(1.0/sqrt(0.5));
        }
    }

    ofstream ofs("output.dat");
    //p_order_N(knot_vector, index, order, coodinate)
    vector<double> knot_vector = {0,0,0,1,1,2,2,3,3,4,4,4};
    for(double i=0.0; i<=4.0; i+=0.01){
        double r_sum_X=0.0;
        double r_sum_Y=0.0;
        for(int j=0; j<control_point.size(); j++){
            double lower_sum=0.0;
            for(int k=0; k<control_point.size(); k++){
                lower_sum += p_order_N(knot_vector, k+1, order, i) * weight[k];
            }
            double tmp_R = p_order_N(knot_vector, j+1, order, i) * weight[j];
            tmp_R /= lower_sum;

            r_sum_X += tmp_R*control_point[j][0];
            r_sum_Y += tmp_R*control_point[j][1];
        }
        ofs << r_sum_X << " " << r_sum_Y << endl;
    }
}