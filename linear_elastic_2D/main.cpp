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
                cout << node[element[ic][p]][0] << endl;
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
        //cout << err << endl;
        if (err < convergence){
            break;
        }
    }
}
void calc_dRdr(vector<vector<double>> &dRdr, int ic, vector<vector<int>> element, vector<vector<int>> element_xi, \
vector<vector<int>> element_eta, vector<double> knot_vector_xi, vector<double> knot_vector_eta, \
vector<double> xi_coordinate, vector<double> eta_coordinate, double gauss_1, double gauss_2, int order_xi, int order_eta, int xi_loop_num, int eta_loop_num, vector<double> control_weight)
{   
    for(int dir=0; dir<2; dir++){
        for(int i=0; i<element_eta[eta_loop_num].size(); i++){//eta dir loop 3
            for(int j=0; j<element_xi[xi_loop_num].size(); j++){//xi dir loop 4
                //xi direction
                double sum_W=0.0;
                double sum_dWdr=0.0;
                for(int k=0; k<element_eta[eta_loop_num].size(); k++){ //eta dir loop 3
                    for(int l=0; l<element_xi[xi_loop_num].size(); l++){ //xi dir loop 3
                        if(dir==0){
                            sum_dWdr += p_order_dNdr(knot_vector_xi,element_xi[xi_loop_num][l]+1,order_xi,gauss_1)\
                                *p_order_N(knot_vector_eta,element_eta[eta_loop_num][k]+1,order_eta,gauss_2)*control_weight[element[ic][k*2+l]];
                            sum_W += p_order_N(knot_vector_xi,element_xi[xi_loop_num][l]+1,order_xi,gauss_1)\
                                * p_order_N(knot_vector_eta,element_eta[eta_loop_num][k]+1,order_eta,gauss_2)*control_weight[element[ic][k*2+l]];
                        }
                        if(dir==1){
                            sum_dWdr +=  p_order_N(knot_vector_xi,element_xi[xi_loop_num][l]+1,order_xi,gauss_1)\
                                *p_order_dNdr(knot_vector_eta,element_eta[eta_loop_num][k]+1,order_eta,gauss_2)*control_weight[element[ic][k*2+l]];
                            sum_W += p_order_N(knot_vector_xi,element_xi[xi_loop_num][l]+1,order_xi,gauss_1)\
                                * p_order_N(knot_vector_eta,element_eta[eta_loop_num][k]+1,order_eta,gauss_2)*control_weight[element[ic][k*2+l]];
                        }
                    }
                }
                double tmp_dRdr = 0.0;
                if(dir==0){
                    tmp_dRdr = p_order_dNdr(knot_vector_xi,element_xi[xi_loop_num][j]+1,order_xi,gauss_1)\
                        *p_order_N(knot_vector_eta,element_eta[eta_loop_num][i]+1,order_eta,gauss_2)*sum_W\
                        -p_order_N(knot_vector_xi,element_xi[xi_loop_num][j]+1,order_xi,gauss_1)*\
                        p_order_N(knot_vector_eta,element_eta[eta_loop_num][i]+1,order_eta,gauss_2)*sum_dWdr;
                }
                if(dir==1){
                    tmp_dRdr = p_order_N(knot_vector_xi,element_xi[xi_loop_num][j]+1,order_xi,gauss_1)\
                        *p_order_dNdr(knot_vector_eta,element_eta[eta_loop_num][i]+1,order_eta,gauss_2)*sum_W\
                        -p_order_N(knot_vector_xi,element_xi[xi_loop_num][j]+1,order_xi,gauss_1)*\
                        p_order_N(knot_vector_eta,element_eta[eta_loop_num][i]+1,order_eta,gauss_2)*sum_dWdr;
                }
                tmp_dRdr/=(sum_W*sum_W);
                dRdr[i*element_xi[xi_loop_num].size()+j][dir] = tmp_dRdr*control_weight[element[ic][i*2+j]];
            }
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
void input_knot_vector(string knot_vector_xi_file,string knot_vector_eta_file,vector<double> &knot_vector_xi, vector<double> &knot_vector_eta)
{
    ifstream ifs(knot_vector_xi_file);
    string str;
    while(getline(ifs,str)){
        knot_vector_xi.push_back(stod(str));
    }
    ifs.close();
    ifs.open(knot_vector_eta_file);
    while(getline(ifs,str)){
        knot_vector_eta.push_back(stod(str));
    }
    ifs.close();
}
void input_control_weight(string control_weight_file, vector<double> &control_weight)
{
    ifstream ifs(control_weight_file);
    string str;
    while(getline(ifs,str)){
        control_weight.push_back(stod(str));
    }
    ifs.close();
}

void return_D_marix(vector<vector<double>> &D)
{
    double Young_modulus = 2e9;
    double Poisson_ratio = 0.3;

    double E = Young_modulus/(1e0-(Poisson_ratio*Poisson_ratio));

    D[0][0] = E; D[0][1] = E*Poisson_ratio; D[0][2] = 0e0;
    D[1][0] = E*Poisson_ratio; D[1][1] = E; D[1][2] = 0e0;
    D[2][0] = 0e0; D[2][1] = 0e0; D[2][2] = E*(1e0-Poisson_ratio)/2e0;
}

vector<vector<double>> calc_B_matrix(vector<int> element, vector<vector<double>> dNdx)
{
    vector<vector<double>> B(3, vector<double>(element.size()*2, 0.0));
    for(int i=0; i<element.size(); i++){
        B[0][2*i] = dNdx[i][0];
        B[2][2*i] = dNdx[i][1];
        B[1][2*i+1] = dNdx[i][1];
        B[2][2*i+1] = dNdx[i][0];
    }
    return B;
}

vector<vector<double>> Transposed_mat(vector<vector<double>> C)
{
    vector<vector<double>> tmp_C(C[0].size(), vector<double>(C.size()));
    for(int i=0; i<C.size(); i++){
        for(int j=0; j<C[0].size(); j++){
            tmp_C[j][i] = C[i][j];
        }
    }
    return tmp_C;
}

vector<vector<double>> mat_product(vector<vector<double>> C1, vector<vector<double>> C2)
{
    vector<vector<double>> tmp_C(C1.size(), vector<double>(C2[0].size()));
    for(int i=0; i<C1.size(); i++){
        for(int j=0; j<C2[0].size(); j++){
            for(int k=0; k<C1[0].size(); k++){
                tmp_C[i][j] += C1[i][k] * C2[k][j];
            }
        }
    }
    return tmp_C;
}

void update_control_coordinate(vector<vector<double>> &x, vector<double> u)
{
    for(int i=0; i<x.size(); i++){
        x[i][0] = x[i][0] + u[2*i];
        x[i][1] = x[i][1] + u[2*i+1];
    }
}

//knot vectorの個数=制御点の個数+order+1
int main()
{
    int order_xi = 1;
    int order_eta = 2;
    string input_dir = "input_file";
    vector<double> knot_vector_xi, knot_vector_eta;
    string knot_vector_xi_file = input_dir + "/" + "knot_vector_xi.dat";
    string knot_vector_eta_file = input_dir + "/" + "knot_vector_eta.dat";
    string control_file = input_dir + "/" + "control.dat";
    string element_file = input_dir + "/" + "element.dat";
    string knot_span_xi_file = input_dir + "/" + "knot_span_xi.dat";
    string knot_span_eta_file = input_dir + "/" + "knot_span_eta.dat";
    string knot_connectivity_file = input_dir + "/" + "knot_connectivity.dat";
    string xi_coordinate_file = input_dir + "/" + "xi_coordinate.dat";
    string eta_coordinate_file = input_dir + "/" + "eta_coordinate.dat";
    string control_weight_file = input_dir + "/" + "control_weight.dat";
    input_knot_vector(knot_vector_xi_file,knot_vector_eta_file,knot_vector_xi, knot_vector_eta);
    vector<vector<double>> control_point;
    input_control(control_file,control_point);
    vector<vector<int>> element;
    input_element(element_file,element);
    vector<vector<int>> element_xi, element_eta;
    input_knot_span(knot_span_xi_file, knot_span_eta_file, element_xi, element_eta);
    vector<double> xi_coordinate, eta_coordinate;
    input_parametric_coordinate(xi_coordinate_file, eta_coordinate_file, xi_coordinate, eta_coordinate);
    vector<int> knot_connectivity_xi, knot_connectivity_eta;
    input_knot_connectivity(knot_connectivity_file, knot_connectivity_xi, knot_connectivity_eta);
    vector<double> control_weight;
    input_control_weight(control_weight_file, control_weight);
    vector<double> gauss = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
    vector<double> gauss_weight = {5.0/9.0, 5.0/9.0, 5.0/9.0};
    vector<vector<double>> K(control_point.size()*2, vector<double>(control_point.size()*2));
    vector<vector<double>> D(3, vector<double>(3, 0.0));
    return_D_marix(D);

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
                        eta_coordinate, gauss_xi, gauss_eta, order_xi, order_eta,  xi_loop_num, eta_loop_num, control_weight);
                        calc_dxdr(i, control_point, element, dxdr, dRdr);
                        calc_inverse_matrix_2x2(dxdr, drdx);
                        double detJ1 = fabs(dxdr[0][0] * dxdr[1][1]  - dxdr[1][0] * dxdr[0][1]);
                        calc_dNdx(i, element, dNdx, dRdr, drdx);
                        double detJ2 = 0.5*(eta_coordinate[eta_loop_num+1]-eta_coordinate[eta_loop_num])*\
                            0.5*(xi_coordinate[xi_loop_num+1]-xi_coordinate[xi_loop_num]);
                        vector<vector<double>> B = calc_B_matrix(element[i], dNdx);
                        vector<vector<double>> B_t = Transposed_mat(B);
                        vector<vector<double>> K_e = mat_product(mat_product(B_t, D), B);
                        for(int n=0; n<2; n++){
                            for(int o=0; o<2; o++){
                                K[element[i][j]*2+o][element[i][k]*2+n] += K_e[j*2+o][k*2+n]*detJ1*detJ2*gauss_weight[l]*gauss_weight[m];
                            }
                        }
                    }
                }
            }
        }
    }

    for(int i=0; i<K.size(); i++){
        for(int j=0; j<K[i].size(); j++){
            cout << K[i][j] << " ";
        }
        cout << endl;
    }

    ofstream check("check.csv");
    for(int i=0; i<K.size(); i++){
        for(int j=0; j<K[i].size(); j++){
            check << K[i][j] << ",";
        }
        check << endl;
    }
    check.close();

    vector<double> f(16, 0.0);
    f[13] = -30000.0;
    f[15] = -30000.0;

    vector<double> u(16, 0.0);

    for(int i=0; i<K.size(); i++){
        K[0][i] = 0e0;
        K[1][i] = 0e0;
        K[2][i] = 0e0;
        K[3][i] = 0e0;
    }
    K[0][0] = 1e0;
    K[1][1] = 1e0;
    K[2][2] = 1e0;
    K[3][3] = 1e0;

    double convegence=0.00001;
    Jacobi_method(u, f, K, convegence);
    for(int i=0; i<u.size(); i++){
        cout << u[i] << endl;
    }   

    
    update_control_coordinate(control_point, u);

    ofstream output("output.dat");
    for(double i=0; i<=1.0; i+=0.01){
        cout << i << endl;
        for(double j=0; j<=1.0; j+=0.01){
            double r_sum_ux=0.0;
            double r_sum_uy=0.0;
            double r_sum_X=0.0;
            double r_sum_Y=0.0;
            for(int k=0; k<4; k++){
                for(int l=0; l<2; l++){
                    double lower_sum=0.0;
                    for(int m=0; m<4; m++){
                        for(int n=0; n<2; n++){
                            lower_sum += p_order_N(knot_vector_xi, n+1, order_xi, i) * p_order_N(knot_vector_eta, m+1, order_eta, j)*control_weight[m*2+n];
                        }
                    }
                    double tmp_R = p_order_N(knot_vector_xi, l+1, order_xi, i)*p_order_N(knot_vector_eta, k+1, order_eta, j)*control_weight[k*2+l];
                    tmp_R /= lower_sum;

                    r_sum_ux += tmp_R * u[(k*2+l)*2];
                    r_sum_uy += tmp_R * u[(k*2+l)*2+1];
                    r_sum_X += tmp_R * (control_point[k*2+l][0]+u[(k*2+l)*2]);
                    r_sum_Y += tmp_R * (control_point[k*2+l][1]+u[(k*2+l)*2+1]);
                }
            }
            output << r_sum_X << " " << r_sum_Y << " " << r_sum_ux << " " << r_sum_uy << endl;
        }
    }

    double r_sum_ux=0.0;
    double r_sum_uy=0.0;
    double r_sum_X=0.0;
    double r_sum_Y=0.0;
    for(int k=0; k<4; k++){
        for(int l=0; l<2; l++){
            double lower_sum=0.0;
            for(int m=0; m<4; m++){
                for(int n=0; n<2; n++){
                    lower_sum += p_order_N(knot_vector_xi, n+1, order_xi, 0.0) * p_order_N(knot_vector_eta, m+1, order_eta, 0.0)*control_weight[m*2+n];
                }
            }
            double tmp_R = p_order_N(knot_vector_xi, l+1, order_xi, 0.0)*p_order_N(knot_vector_eta, k+1, order_eta, 0.0)*control_weight[k*2+l];
            tmp_R /= lower_sum;
            r_sum_ux += tmp_R * u[(k*2+l)*2];
            r_sum_uy += tmp_R * u[(k*2+l)*2+1];
            r_sum_X += tmp_R * (control_point[k*2+l][0]+u[(k*2+l)*2]);
            r_sum_Y += tmp_R * (control_point[k*2+l][1]+u[(k*2+l)*2+1]);
        }
    }
    cout << r_sum_X << " " << r_sum_Y << " " << r_sum_ux << " " << r_sum_uy << endl;
}