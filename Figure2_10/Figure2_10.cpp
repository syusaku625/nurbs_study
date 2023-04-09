#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<sstream>

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
    //cout << N << endl;
    return N;
}

int main()
{
    ifstream ifs("Fig_2_12_control_point.dat");
    vector<double> knot_vector = {0,0,0,1,2,3,4,4,5,5,5};
    vector<vector<double>> control_point;
    string str;
    while(getline(ifs,str)){
        istringstream ss(str);
        vector<double> point;
        for(int i=0; i<2; i++){
            getline(ss,str,' ');
            point.push_back(stod(str));
        }
        control_point.push_back(point);
    }
    ifs.close();

    //p_order_N(knot_vector, index, order, coodinate)
    int order = 2;
    ofstream ofs("fig_2_12_order2.dat");
    for(double i=0.0; i<=5.0; i+=0.01){
        cout << i << endl;
        double sum_X = 0.0;
        double sum_Y = 0.0;
        for(int j=0; j<control_point.size(); j++){
            sum_X += p_order_N(knot_vector, j+1, order, i)*control_point[j][0];
            sum_Y += p_order_N(knot_vector, j+1, order, i)*control_point[j][1];
        }
        ofs << sum_X << " " << sum_Y << endl;
    }
}