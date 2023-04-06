#include<iostream>
#include<vector>
#include<math.h>
#include<string>
#include<algorithm>
#include<fstream>

using namespace std;

double return_N0(int index, double c)
{
    if(fabs(c-index)<0.00001){
        return 1.0;
    }
    else if(c-(index-1.0)<0.00001){
        return 0.0;
    }
    else if(c<index){
        return 1.0;
    }
    else{
        return 0.0;
    }
}

double p_order_N(vector<double> knot_vector, int index, int order, double c)
{
    if(order==0){
        double N = return_N0(index, c);
        return N;
    }

    double N = (c-knot_vector[index-1])/(knot_vector[index-1+order]-knot_vector[index-1])*p_order_N(knot_vector, index, order-1, c)\
    +(knot_vector[index-1+order+1]-c)/(knot_vector[index-1+order+1]-knot_vector[index-1+1])*p_order_N(knot_vector, index+1, order-1, c);
    return N;
}

int main()
{
    vector<double> knot_vector = {0,1,2,3,4,5};
    ofstream ofs("test.dat");
    for(double i=0.0; i<=5.0; i+=0.01){
        ofs << i << " " << p_order_N(knot_vector, 2, 3, i) << endl; //p_order_N(knot_vector, index, order, coodinate)
    }
    ofs.close();
}