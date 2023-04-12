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
    ifstream ifs("control_point.dat");
    string str;
    vector<vector<vector<vector<double>>>> control_point(
        3, vector<vector<vector<double>>>(
            9, vector<vector<double>>(
                2, vector<double>(
                    3
                )
            )
        )
    );

    for(int i=0; i<2; i++){
        for(int j=0; j<9; j++){
            for(int k=0; k<3; k++){
                getline(ifs,str); 
                istringstream ss(str);
                for(int l=0; l<3; l++){
                    getline(ss, str, ' ');
                    control_point[k][j][i][l]=stod(str);
                }
            }
        }
    }

    ifs.close();
    ifs.open("weight.dat");

    vector<vector<vector<double>>> weight(
        3, vector<vector<double>>(
            9, vector<double>(
                2
            )
        )
    );

    for(int i=0; i<2; i++){
        for(int j=0; j<9; j++){
            for(int k=0; k<3; k++){
                getline(ifs,str);
                if(str=="circle"){
                    weight[k][j][i] = 1.0;
                }
                else{
                    weight[k][j][i] = 1.0/sqrt(2);
                }
            }
        }
    }

    ofstream ofs("output.dat");

    int order_xi=2;
    int order_eta=1;
    int order_jita=2;

    vector<double> knot_vector_xi = {0,0,0,1,1,2,2,3,3,4,4,4};
    vector<double> knot_vector_eta = {0,0,1,1};
    vector<double> knot_vector_jita = {0,0,0,1,1,2,2,2};

    for(double i=0.0; i<=4.0; i+=0.1){
        cout << i << endl;
        for(double j=0.0; j<=1.0; j+=0.1){
            for(double k=0.0; k<=2.0; k+=0.1){
                double r_sum_X=0.0;
                double r_sum_Y=0.0;
                double r_sum_Z=0.0;
                for(int l=0; l<control_point.size(); l++){
                    for(int m=0; m<control_point[l].size(); m++){
                        for(int n=0; n<control_point[l][m].size(); n++){
                            double lower_sum=0.0;
                            for(int o=0; o<control_point.size(); o++){ //z
                                for(int p=0; p<control_point[o].size(); p++){ //xy
                                    for(int q=0; q<control_point[o][p].size(); q++){  //r
                                        lower_sum += p_order_N(knot_vector_jita, o+1, order_jita, k) * \
                                        p_order_N(knot_vector_xi, p+1, order_xi, i) * \
                                        p_order_N(knot_vector_eta, q+1, order_eta, j) * weight[o][p][q];
                                    }
                                }
                            }
                            double tmp_R = p_order_N(knot_vector_jita, l+1, order_jita, k) * \
                            p_order_N(knot_vector_xi, m+1, order_xi, i) * \
                            p_order_N(knot_vector_eta, n+1, order_eta, j) * weight[l][m][n];
                            tmp_R /= lower_sum;

                            r_sum_X += tmp_R * control_point[l][m][n][0];
                            r_sum_Y += tmp_R * control_point[l][m][n][1];
                            r_sum_Z += tmp_R * control_point[l][m][n][2];
                        }
                    }
                }
                if(fabs(j)<0.1) ofs << r_sum_X << " " << r_sum_Y << " " << r_sum_Z << endl;
                if(fabs(j-1.0)<0.1) ofs << r_sum_X << " " << r_sum_Y << " "  << r_sum_Z << endl;
            }
        }
    }
    ofs.close();
}