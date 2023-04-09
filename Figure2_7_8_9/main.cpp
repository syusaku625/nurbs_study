#include<iostream>
#include<vector>
#include<math.h>
#include<string>
#include<algorithm>
#include<fstream>

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

    cout << N << endl;

    return N;
}

int main()
{
    vector<double> knot_vector = {0,0,0,1,2,3,4,4,5,5,5};
    int order = 4;
    double to = 5.0;
    ofstream ofs("N_1_2.dat");
    for(double i=0.0; i<=to; i+=0.01){
        ofs << i << " " << p_order_N(knot_vector, 1, order, i) << endl; //p_order_N(knot_vector, index, order, coodinate)
    }
    ofs.close();

    ofs.open("N_2_2.dat");
    for(double i=0.0; i<=to; i+=0.01){
        ofs << i << " " << p_order_N(knot_vector, 2, order, i) << endl; //p_order_N(knot_vector, index, order, coodinate)
    }
    ofs.close();

    ofs.open("N_3_2.dat");
    for(double i=0.0; i<=to; i+=0.01){
        ofs << i << " " << p_order_N(knot_vector, 3, order, i) << endl; //p_order_N(knot_vector, index, order, coodinate)
    }
    ofs.close();

    ofs.open("N_4_2.dat");
    for(double i=0.0; i<=to; i+=0.01){
        ofs << i << " " << p_order_N(knot_vector, 4, order, i) << endl; //p_order_N(knot_vector, index, order, coodinate)
    }
    ofs.close();

    ofs.open("N_5_2.dat");
    for(double i=0.0; i<=to; i+=0.01){
        ofs << i << " " << p_order_N(knot_vector, 5, order, i) << endl; //p_order_N(knot_vector, index, order, coodinate)
    }
    ofs.close();

    ofs.open("N_6_2.dat");
    for(double i=0.0; i<=to; i+=0.01){
        ofs << i << " " << p_order_N(knot_vector, 6, order, i) << endl; //p_order_N(knot_vector, index, order, coodinate)
    }
    ofs.close();

    ofs.open("N_7_2.dat");
    for(double i=0.0; i<=to; i+=0.01){
        ofs << i << " " << p_order_N(knot_vector, 7, order, i) << endl; //p_order_N(knot_vector, index, order, coodinate)
    }
    ofs.close();

    ofs.open("N_8_2.dat");
    for(double i=0.0; i<=to; i+=0.01){
        ofs << i << " " << p_order_N(knot_vector, 8, order, i) << endl; //p_order_N(knot_vector, index, order, coodinate)
    }
    ofs.close();
}