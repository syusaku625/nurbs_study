#include<iostream>
#include<vector>
#include<fstream>
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

double p_order_dNdr(vector<double> knot_vector, int index, int order, double c)
{
    if(order==0){
        return 0.0;
    }
    double N1, N2;
    if(knot_vector[index-1+order]-knot_vector[index-1]==0){
        N1 = 0.0;
    }
    else{
        N1 = order/(knot_vector[index-1+order]-knot_vector[index-1]);
    }
    if(knot_vector[index-1+order+1]-knot_vector[index-1+1]==0){
        N2 = 0.0;
    }
    else{
        N2 = order/(knot_vector[index-1+order+1]-knot_vector[index-1+1]);
    }
    double N = N1 * p_order_N(knot_vector, index, order-1, c) - N2 * p_order_N(knot_vector, index+1, order-1, c);
    return N;
}

void Jacobi_method(std::vector<double> &p, vector<double> &R, vector<vector<double>> &G, double convergence)
{
    vector<double> tmp_p;
    int NN=p.size();
    while(1){
        tmp_p.resize(NN);
        for (int i = 0; i < NN; i++)
        {
            tmp_p[i]=R[i];
            for (int j = 0; j < NN; j++){
                if(i!=j){
                    tmp_p[i] -= G[i][j] * tmp_p[j];
                }
            }
            tmp_p[i] /= G[i][i];
        }
        double err = 0.0;
        for (int i = 0; i < NN; i++){
            err += fabs(tmp_p[i] - p[i]);
            p[i] = tmp_p[i];
        }
        if (err < convergence){
            break;
        }
    }
}

void calc_dRdr(vector<double> &dRdr, int ic, vector<vector<int>> element, vector<double> knot_vector, vector<double> xi_coordinate, double gauss, int order)
{
    for(int l=0; l<element[ic].size(); l++){
        double gauss_coordinate = 0.5*((xi_coordinate[ic+1]-xi_coordinate[ic])*gauss+(xi_coordinate[ic]+xi_coordinate[ic+1]));
        double sum_W=0.0;
        double sum_dWdr=0.0;
        for (int p=0; p<element[ic].size(); p++){
            sum_dWdr += p_order_dNdr(knot_vector,element[ic][p]+1,order,gauss_coordinate);
            sum_W += p_order_N(knot_vector,element[ic][p]+1,order,gauss_coordinate);
        }
        double tmp_dRdr = p_order_dNdr(knot_vector,element[ic][l]+1,order,gauss_coordinate)*sum_W-p_order_N(knot_vector,element[ic][l]+1,order,gauss_coordinate)*sum_dWdr;
        tmp_dRdr /= (sum_W*sum_W);
        dRdr[l] = tmp_dRdr;
    }
}

int main()
{
    int order = 2;
    vector<double> knot_vector = {0.0,0.0,0.0,0.5,1.0,1.0,1.0};

    vector<vector<double>> control_point(4, vector<double>(2));
    control_point[0][0] = 0.0;
    control_point[1][0] = 0.2;
    control_point[2][0] = 0.8;
    control_point[3][0] = 1.0;
    
    vector<vector<int>> element(2, vector<int>(3));
    element[0][0] = 0;
    element[0][1] = 1;
    element[0][2] = 2;
    element[1][0] = 1;
    element[1][1] = 2;
    element[1][2] = 3;

    vector<double> gauss = {-0.86113, -0.33998, 0.33998, 0.86113};
    vector<double> weight = {0.34785, 0.65214, 0.65214, 0.34785};

    vector<double> xi_coordinate = {0.0, 0.5, 1.0};

    vector<vector<double>> K(4, vector<double>(4));

    for(int i=0; i<element.size(); i++){//elementループ
        for(int row=0; row<element[i].size(); row++){//control_pointのループ
            for(int col=0; col<element[i].size(); col++){//control_pointのループ
                for(int j=0; j<gauss.size(); j++){//ガウスポイントループ
                    vector<double> dNdx(element[i].size());
                    vector<double> dRdr(element[i].size());

                    calc_dRdr(dRdr, i, element, knot_vector, xi_coordinate, gauss[j], order);
                    
                    double detJ1=0.0;
                    for(int l=0; l<element[i].size(); l++){
                        detJ1 += dRdr[l] * control_point[element[i][l]][0];
                    }
                    double detJ2 = 0.25*(xi_coordinate[i+1]-xi_coordinate[i]);

                    for(int l=0; l<element[i].size(); l++){
                        dNdx[l] = (1.0/detJ1)*dRdr[l];
                    }
                    K[element[i][row]][element[i][col]] += dNdx[row]*dNdx[col]*detJ1*detJ2*weight[j];
                }
            }
        }
    }    
    
    vector<double> b = {0.0,0.0,0.0,0.5};
    for(int i=0; i<K.size(); i++){
        K[3][i]=0.0;
        K[0][i]=0.0;
        K[0][0]=1.0;
        K[3][3]=1.0;
    }

    double convegence=0.00001;

    vector<double> u(4);
    Jacobi_method(u, b, K, convegence);
    for(int i=0; i<u.size(); i++){
        cout << u[i] << endl;
    }

    ofstream ofs("output_p.dat");
    for(double i=0.0; i<=1.0; i+=0.01){
        double p_sum=0.0;
        for(int l=0; l<control_point.size(); l++){
            double lower_sum=0.0;
            for(int o=0; o<control_point.size(); o++){
                lower_sum += p_order_N(knot_vector, o+1, order, i);
            }
            double tmp_R = p_order_N(knot_vector, l+1, order, i);
            tmp_R /= lower_sum;
            
            p_sum += tmp_R * u[l];
        }
        ofs << p_sum <<  endl;
    }
    ofs.close();
}