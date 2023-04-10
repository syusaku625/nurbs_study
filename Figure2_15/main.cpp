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
    vector<double> knot_vector_x = {0,0,0,0.5,1,1,1};
    vector<double> knot_vector_y = {0,0,0,1,1,1};
    int order_x = 2;
    int order_y = 2;
    double to_x = 1.0;
    double to_y = 1.0;
    vector<vector<vector<double>>> control_point(4, vector<vector<double>>(3, vector<double>(2,0)));
    string str;
    for(int i=0; i<4; i++){
        getline(ifs,str);
        istringstream ss(str);
        for(int j=0; j<3; j++){
            for(int k=0; k<2; k++){
                getline(ss,str,' ');
                control_point[i][j][k] = stod(str);
            }
        }
    }
    ifs.close();

    for(int i=0; i<4; i++){
        for(int j=0; j<3; j++){
            cout << control_point[i][j][0] << " " << control_point[i][j][1] << endl;
        }
    }


    ofstream ofs1("order1.dat");
    ofstream ofs2("order2.dat");
    ofstream ofs3("order3.dat");
    ofstream ofs4("order4.dat");
    ofstream ofs5("order5.dat");

    for(double i=0.0; i<=to_x; i+=0.0001){
        cout << i << endl;
        for(double j=0.0; j<=to_y; j+=0.0001){
            double sum_X_1 = 0.0;
            double sum_Y_1 = 0.0;
            for(int k=0; k<4; k++){
                for(int l=0; l<3; l++){
                    sum_X_1 += p_order_N(knot_vector_x, k+1, order_x, i) * p_order_N(knot_vector_y, l+1, order_y, j) * control_point[k][l][0];
                    sum_Y_1 += p_order_N(knot_vector_x, k+1, order_x, i) * p_order_N(knot_vector_y, l+1, order_y, j) * control_point[k][l][1];
                }
            }
            if(fabs(i)<0.0001) ofs1 << sum_X_1 << " " << sum_Y_1 << endl;
            if(fabs(i-1.0)<0.0001)  ofs2 << sum_X_1 << " " << sum_Y_1 << endl;
            if(fabs(j)<0.0001) ofs3 << sum_X_1 << " " << sum_Y_1 << endl;
            if(fabs(j-1.0)<0.0001) ofs4 << sum_X_1 << " " << sum_Y_1 << endl;
            if(fabs(i-0.5)<0.00001) ofs5 << sum_X_1 << " " << sum_Y_1 << endl;
        }
    }

    //ofs.open("control_point_plot.dat");
    //for(int i=0; i<4; i++){
    //    for(int j=0; j<3; j++){
    //        ofs << control_point[i][j][0] << " " << control_point[i][j][1] << endl;
    //    }
    //}
    //ofs.close();
}
