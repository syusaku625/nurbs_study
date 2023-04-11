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
    vector<vector<vector<double>>> control_point(9, vector<vector<double>>(2, vector<double>(2)));
    for(int i=0; i<2; i++){
        for(int j=0; j<9; j++){
            getline(ifs,str); 
            istringstream ss(str);
            for(int k=0; k<2; k++){
                getline(ss, str, ' ');
                control_point[j][i][k]=stod(str);
            }
        }
    }
    ifs.close();

    ifs.open("weight.dat");
    vector<vector<double>> weight(9, vector<double>(2));
    for(int i=0; i<2; i++){
        for(int j=0; j<9; j++){
            getline(ifs,str);
            if(str=="circle"){
                weight[j][i] = 1.0;
            }
            else{
                weight[j][i] = 1.0/sqrt(2.0);
            }
            cout << weight[j][i] << endl;
        }
    }

    ofstream ofs("output2.dat");
    
    //p_order_N(knot_vector, index, order, coodinate)
    int order_xi = 2;
    int order_eta = 1;
    vector<double> knot_vector_xi = {0,0,0,1,1,2,2,3,3,4,4,4};
    vector<double> knot_vector_eta = {0,0,1,1};

    for(double i=0.0; i<=4.0; i+=0.01){
        cout << i << endl;
        for(double j=0.0; j<=1.0; j+=0.01){
            double r_sum_X=0.0;
            double r_sum_Y=0.0;
            for(int k=0; k<control_point.size(); k++){
                for(int l=0; l<control_point[k].size(); l++){
                    double lower_sum=0.0;
                    for(int m=0; m<control_point.size(); m++){
                        for(int n=0; n<control_point[m].size(); n++){
                            lower_sum += p_order_N(knot_vector_xi, m+1, order_xi, i) * p_order_N(knot_vector_eta, n+1, order_eta, j) * weight[m][n];
                        }
                    }
                    double tmp_R = p_order_N(knot_vector_xi, k+1, order_xi, i)*p_order_N(knot_vector_eta, l+1, order_eta, j) * weight[k][l];
                    tmp_R /= lower_sum;

                    r_sum_X += tmp_R * control_point[k][l][0];
                    r_sum_Y += tmp_R * control_point[k][l][1];
                }
            }
            if(fabs(j)<0.01) ofs << r_sum_X << " " << r_sum_Y << endl;
            if(fabs(j-1.0)<0.01) ofs << r_sum_X << " " << r_sum_Y << endl;
        }
    }
    ofs.close();
}