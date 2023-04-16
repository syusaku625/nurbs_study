#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
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
        N1 = double(order)/(knot_vector[index-1+order]-knot_vector[index-1]);
    }
    if(knot_vector[index-1+order+1]-knot_vector[index-1+1]==0){
        N2 = 0.0;
    }
    else{
        N2 = double(order)/(knot_vector[index-1+order+1]-knot_vector[index-1+1]);
    }
    double N = N1 * p_order_N(knot_vector, index, order-1, c) - N2 * p_order_N(knot_vector, index+1, order-1, c);
    return N;
}

void calc_dxdr(int ic, vector<vector<double>> node, vector<vector<int>> element, vector<vector<double>> &dxdr, vector<vector<double>> dRdr)
{
    for(int k=0;k<2;k++){
        for(int l=0;l<2;l++){
            dxdr[k][l] = 0e0;
            for(int p=0;p<element[ic].size();p++){
                dxdr[k][l] += dRdr[p][k] * node[element[ic][p]][l];
            }
        }
    }
}

void calc_inverse_matrix_2x2(std::vector<std::vector<double>> dxdr, std::vector<std::vector<double>> &drdx)
{
  double det = dxdr[0][0]*dxdr[1][1]-dxdr[0][1]*dxdr[1][0];
  drdx[0][0] = 1.0/det*dxdr[1][1];
  drdx[1][1] = 1.0/det*dxdr[0][0];
  drdx[0][1] = -1.0/det*dxdr[1][0];
  drdx[1][0] = -1.0/det*dxdr[0][1];
}

void calc_dNdx(int ic, vector<vector<int>> IGA_element, vector<vector<double>> &dNdx, vector<vector<double>> dNdr, vector<vector<double>> drdx)
{
    for(int k=0; k<IGA_element[ic].size(); k++){
      for(int l=0; l<2; l++){
        dNdx[k][l] = 0.0;
        for(int p=0; p<2; p++){
          dNdx[k][l] += dNdr[k][p]*drdx[p][l];
        }
      }
    }
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

void calc_dRdr(vector<vector<double>> &dRdr, int ic, vector<vector<int>> element, vector<vector<int>> element_xi, \
vector<vector<int>> element_eta, vector<double> knot_vector_xi, vector<double> knot_vector_eta, \
vector<double> xi_coordinate, vector<double> eta_coordinate, double gauss_1, double gauss_2, int order_xi, int order_eta)
{   
    for(int i=0; i<3; i++){//eta dir loop 3
        for(int j=0; j<2; j++){//xi dir loop 4
            //xi direction
            double sum_W=0.0;
            double sum_dWdr=0.0;
            for(int k=0; k<3; k++){ //eta dir loop 3
                for(int l=0; l<2; l++){ //xi dir loop 3
                    cout << l+1 << " " << element_eta[ic][k]+1 << endl;
                    sum_dWdr +=  p_order_dNdr(knot_vector_xi,l+1,order_xi,gauss_1)\
                    *p_order_N(knot_vector_eta,element_eta[ic][k]+1,order_eta,gauss_2);
                    sum_W += p_order_N(knot_vector_xi,l+1,order_xi,gauss_1)\
                     * p_order_N(knot_vector_eta,element_eta[ic][k]+1,order_eta,gauss_2);
                }
            }
            double tmp_dRdr = p_order_dNdr(knot_vector_xi,j+1,order_xi,gauss_1)\
            *p_order_N(knot_vector_eta,element_eta[ic][i]+1,order_eta,gauss_2)*sum_W\
            -p_order_N(knot_vector_xi,j+1,order_xi,gauss_1)*\
            p_order_N(knot_vector_eta,element_eta[ic][i]+1,order_eta,gauss_2)*sum_dWdr;

            tmp_dRdr/=(sum_W*sum_W);
            dRdr[i*2+j][0] = tmp_dRdr;

            //eta direction
            sum_W=0.0;
            sum_dWdr=0.0;
            for(int k=0; k<3; k++){ //eta dir loop 3
                for(int l=0; l<2; l++){ //xi dir loop 3
                    sum_dWdr +=  p_order_N(knot_vector_xi,l+1,order_xi,gauss_1)\
                    *p_order_dNdr(knot_vector_eta,element_eta[ic][k]+1,order_eta,gauss_2);
                    sum_W += p_order_N(knot_vector_xi,l+1,order_xi,gauss_1)\
                     * p_order_N(knot_vector_eta,element_eta[ic][k]+1,order_eta,gauss_2);
                }
            }
            tmp_dRdr = p_order_N(knot_vector_xi,j+1,order_xi,gauss_1)\
            *p_order_dNdr(knot_vector_eta,element_eta[ic][i]+1,order_eta,gauss_2)*sum_W\
            -p_order_N(knot_vector_xi,j+1,order_xi,gauss_1)*\
            p_order_N(knot_vector_eta,element_eta[ic][i]+1,order_eta,gauss_2)*sum_dWdr;
            tmp_dRdr/=(sum_W*sum_W);
            dRdr[i*2+j][1] = tmp_dRdr;
        }
    }
}

void input_control(string filename, vector<vector<double>> &control_point)
{
    ifstream ifs(filename);
    string str;
    while(getline(ifs,str)){
        istringstream ss(str);
        vector<double> tmp_control;
        for(int i=0; i<2; i++){
            getline(ss, str, ' ');
            tmp_control.push_back(stod(str));
        }
        control_point.push_back(tmp_control);
    }
}

void input_element(string filename, vector<vector<int>> &element)
{
    ifstream ifs(filename);
    string str;
    
    while(getline(ifs,str)){
        vector<int> tmp_element;
        istringstream ss(str);
        while(getline(ss, str, ' ')){
            tmp_element.push_back(stoi(str));
        }
        element.push_back(tmp_element);
    }
}

void input_knot_span(string knot_span_xi_file, string knot_span_eta_file, vector<vector<int>> &element_xi, vector<vector<int>> &element_eta)
{
    string str;
    ifstream ifs(knot_span_xi_file);
    while(getline(ifs,str)){
        vector<int> tmp_knot_span;
        istringstream ss(str);
        while(getline(ss, str, ' ')){
            tmp_knot_span.push_back(stoi(str));
        }
        element_xi.push_back(tmp_knot_span);
    }
    ifs.close();

    ifs.open(knot_span_eta_file);
    while(getline(ifs,str)){
        vector<int> tmp_knot_span;
        istringstream ss(str);
        while(getline(ss, str, ' ')){
            tmp_knot_span.push_back(stoi(str));
        }
        element_eta.push_back(tmp_knot_span);
    }
    ifs.close();
}

void input_parametric_coordinate(string xi_coordinate_file, string eta_coordinate_file, vector<double> &xi_coordinate, vector<double> &eta_coordinate)
{
    string str;
    ifstream ifs(xi_coordinate_file);
    while(getline(ifs, str)){
        xi_coordinate.push_back(stod(str));
    }
    ifs.close();
    ifs.open(eta_coordinate_file);
    while (getline(ifs,str)){
        eta_coordinate.push_back(stod(str));
    }
    ifs.close();
}

void input_knot_connectivity(string knot_connectivity_file, vector<int> &knot_connectivity_xi, vector<int> &knot_connectivity_eta)
{
    ifstream ifs(knot_connectivity_file);
    string str;
    while(getline(ifs,str)){
        istringstream ss(str);
        getline(ss, str, ' ');
        knot_connectivity_xi.push_back(stoi(str));
        getline(ss, str, ' ');
        knot_connectivity_eta.push_back(stoi(str));
    }
    ifs.close();
}

int main()
{
    int order_xi = 1;
    int order_eta = 2;
    vector<double> knot_vector_xi={0.0, 0.0, 1.0, 1.0};
    vector<double> knot_vector_eta={0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};
    vector<vector<double>> control_point;
    input_control("control.dat",control_point);
    vector<vector<int>> element;
    input_element("element.dat",element);
    vector<vector<int>> element_xi, element_eta;
    input_knot_span("knot_span_xi.dat", "knot_span_eta.dat", element_xi, element_eta);
    vector<double> xi_coordinate, eta_coordinate;
    input_parametric_coordinate("xi_coordinate.dat", "eta_coordinate.dat", xi_coordinate, eta_coordinate);
    vector<int> knot_connectivity_xi, knot_connectivity_eta;
    input_knot_connectivity("knot_connectivity.dat", knot_connectivity_xi, knot_connectivity_eta);

    vector<double> gauss = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
    vector<double> gauss_weight = {5.0/9.0, 5.0/9.0, 5.0/9.0};
    
    vector<vector<double>> K(control_point.size(), vector<double>(control_point.size()));

    for(int i=0; i<element.size(); i++){//elementループ
        int eta_loop_num = knot_connectivity_eta[i];
        int xi_loop_num = knot_connectivity_xi[i];
        for(int j=0; j<element[i].size(); j++){//control_pointのループ(row)
            for(int k=0; k<element[i].size(); k++){//control_pointのループ(column)
                for(int l=0; l<gauss.size(); l++){//gauss loop 1
                    for(int m=0; m<gauss.size(); m++){//gauss loop 2
                        double gauss_eta = 0.5*((eta_coordinate[eta_loop_num+1]-eta_coordinate[eta_loop_num])*\
                            gauss[l]+(eta_coordinate[eta_loop_num]+eta_coordinate[eta_loop_num+1]));
                        double gauss_xi = 0.5*((xi_coordinate[xi_loop_num+1]-xi_coordinate[xi_loop_num])*\
                            gauss[m]+(xi_coordinate[xi_loop_num]+xi_coordinate[xi_loop_num+1]));
                        vector<vector<double>> dxdr(2, vector<double>(2));
                        vector<vector<double>> drdx(2, vector<double>(2));
                        vector<vector<double>> dNdx(element[i].size(), vector<double>(2));
                        vector<vector<double>> dRdr(element[i].size(), vector<double>(2));
                        calc_dRdr(dRdr, i, element, element_xi, element_eta, knot_vector_xi, knot_vector_eta, xi_coordinate, \
                        eta_coordinate, gauss_xi, gauss_eta, order_xi, order_eta);
                        calc_dxdr(i, control_point, element, dxdr, dRdr);
                        calc_inverse_matrix_2x2(dxdr, drdx);
                        double detJ1 = dxdr[0][0] * dxdr[1][1]  - dxdr[1][0] * dxdr[0][1];
                        calc_dNdx(i, element, dNdx, dRdr, drdx);

                        double detJ2 = 1.0/4.0*(eta_coordinate[knot_connectivity_eta[i]+1]-eta_coordinate[knot_connectivity_eta[i]])*\
                        (xi_coordinate[knot_connectivity_xi[i]+1]-xi_coordinate[knot_connectivity_xi[i]]);

                        K[element[i][j]][element[i][k]] += (dNdx[j][0]*dNdx[k][0]+dNdx[j][1]*dNdx[k][1])*detJ1*detJ2*gauss_weight[l]*gauss_weight[m];
                    }
                }
            }
        }
    }

    for(int i=0; i<8; i++){
        for(int j=0; j<8; j++){
            cout << K[i][j] << " ";
        }
        cout << endl;
    }

    vector<double> b(8);
    b[0] = 1.0;
    b[1] = 1.0;
    b[6] = 0.0;
    b[7] = 0.0;

    for(int i=0; i<K.size(); i++){
        K[0][i] = 0.0;
        K[1][i] = 0.0;
        K[6][i] = 0.0;
        K[7][i] = 0.0;
    }
    K[0][0] = 1.0;
    K[1][1] = 1.0;
    K[6][6] = 1.0;
    K[7][7] = 1.0;

    double convegence=0.00001;

    vector<double> u(8);
    Jacobi_method(u, b, K, convegence);
    for(int i=0; i<u.size(); i++){
        cout << u[i] << endl;
    }

    ofstream ofs("p.dat");
    for(double i=0; i<=1.0; i+=0.01){
        //cout << i << endl;
        for(double j=0; j<=1.0; j+=0.01){
            double r_sum_u=0.0;
            double r_sum_X=0.0;
            double r_sum_Y=0.0;
            for(int k=0; k<4; k++){
                for(int l=0; l<2; l++){
                    double lower_sum=0.0;
                    for(int m=0; m<4; m++){
                        for(int n=0; n<2; n++){
                            lower_sum += p_order_N(knot_vector_xi, n+1, order_xi, i) * p_order_N(knot_vector_eta, m+1, order_eta, j);
                        }
                    }
                    double tmp_R = p_order_N(knot_vector_xi, l+1, order_xi, i)*p_order_N(knot_vector_eta, k+1, order_eta, j);
                    tmp_R /= lower_sum;

                    r_sum_u += tmp_R * u[k*2+l];
                    r_sum_X += tmp_R * control_point[k*2+l][0];
                    r_sum_Y += tmp_R * control_point[k*2+l][1];
                }
            }
            ofs << r_sum_X << " " << r_sum_Y << " " << r_sum_u << endl;
        }
    }
}