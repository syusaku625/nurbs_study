#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include<sstream>

using namespace std;

void input_knot_vector_2D(vector<double> &xi, vector<double> &eta)
{
    string str;
    ifstream ifs("knot_vector_xi.dat");
    while(getline(ifs,str)){
        xi.push_back(stod(str));
    }
    ifs.close();
    ifs.open("knot_vector_eta.dat");
    while(getline(ifs,str)){
        eta.push_back(stod(str));
    }
    ifs.close();
}

void input_control_points_2D(vector<vector<vector<double>>> &control_points, int numOf_fxi, int numOf_feta)
{
    string str;
    ifstream ifs("control.dat");
    for(int i=0; i<numOf_feta; i++){
        for(int j=0; j<numOf_fxi; j++){
            control_points[i][j].resize(2);
            getline(ifs,str);
            istringstream ss(str);
            for(int k=0; k<2; k++){
                getline(ss,str,' ');
                control_points[i][j][k] = stod(str);
            }
        }
    }
    ifs.close();
}

void input_weight(vector<vector<double>> &weight, int numOf_fxi, int numOf_feta)
{
    string str;
    ifstream ifs("weight.dat");
    for(int i=0; i<numOf_feta; i++){
        for(int j=0; j<numOf_fxi; j++){
            getline(ifs,str);
            weight[i][j]=stod(str);
        }
    }
    ifs.close();
}

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

double calc_R_2D(vector<double> knot_vector_xi, vector<double> knot_vector_eta, vector<vector<double>> weight, \
int numOf_f_xi, int numOf_f_eta, int order_xi, int order_eta, int index_xi, int index_eta, double param_xi, double param_eta)
{
    double lower_sum=0.0;
    for(int m=0; m<numOf_f_eta; m++){
        for(int n=0; n<numOf_f_xi; n++){
            lower_sum += p_order_N(knot_vector_xi, n+1, order_xi, param_xi) * p_order_N(knot_vector_eta, m+1, order_eta, param_eta) * weight[m][n];
        }
    }
    double tmp_R = p_order_N(knot_vector_xi, index_xi+1, order_xi, param_xi)*p_order_N(knot_vector_eta, index_eta+1, order_eta, param_eta) * weight[index_eta][index_xi];
    tmp_R /= lower_sum;

    return tmp_R;
}

double calc_dRdr_2D_xi(vector<double> knot_vector_xi, vector<double> knot_vector_eta, vector<vector<double>> weight, \
int numOf_f_xi, int numOf_f_eta, int order_xi, int order_eta, int index_xi, int index_eta, double param_xi, double param_eta)
{
    double lower_sum = 0.0;
    double upper_sum = 0.0;
    for(int m=0; m<numOf_f_eta; m++){
        for(int n=0; n<numOf_f_xi; n++){
            lower_sum += p_order_N(knot_vector_xi, n+1, order_xi, param_xi)*p_order_N(knot_vector_eta,m+1,order_eta,param_eta)*weight[m][n];
            upper_sum += p_order_dNdr(knot_vector_xi, n+1, order_xi, param_xi)*p_order_N(knot_vector_eta, m+1, order_eta, param_eta)*weight[m][n];
        }
    }
    double tmp_dRdr = weight[index_eta][index_xi]*p_order_dNdr(knot_vector_xi, index_xi+1, order_xi, param_xi)*p_order_N(knot_vector_eta, index_eta+1, order_eta, param_eta)\
    -upper_sum*calc_R_2D(knot_vector_xi,knot_vector_eta,weight,numOf_f_xi,numOf_f_eta,order_xi,order_eta,index_xi,index_eta,param_xi,param_eta);
    if(index_xi==3 && index_eta==0){
        cout << p_order_N(knot_vector_xi, index_xi+1, order_xi, param_xi) << endl;
    }
    tmp_dRdr /= lower_sum;

    return tmp_dRdr;
}

double calc_dRdr_2D_eta(vector<double> knot_vector_xi, vector<double> knot_vector_eta, vector<vector<double>> weight, \
int numOf_f_xi, int numOf_f_eta, int order_xi, int order_eta, int index_xi, int index_eta, double param_xi, double param_eta)
{
    double lower_sum = 0.0;
    double upper_sum = 0.0;
    for(int m=0; m<numOf_f_eta; m++){
        for(int n=0; n<numOf_f_xi; n++){
            lower_sum += p_order_N(knot_vector_xi, n+1, order_xi, param_xi)*p_order_N(knot_vector_eta,m+1,order_eta,param_eta)*weight[m][n];
            upper_sum += p_order_N(knot_vector_xi, n+1, order_xi, param_xi)*p_order_dNdr(knot_vector_eta, m+1, order_eta, param_eta)*weight[m][n];
        }
    }
    double tmp_dRdr = weight[index_eta][index_xi]*p_order_N(knot_vector_xi, index_xi+1, order_xi, param_xi)*p_order_dNdr(knot_vector_eta, index_eta+1, order_eta, param_eta)\
    -upper_sum*calc_R_2D(knot_vector_xi,knot_vector_eta,weight,numOf_f_xi,numOf_f_eta,order_xi,order_eta,index_xi,index_eta,param_xi,param_eta);

    tmp_dRdr /= lower_sum;

    return tmp_dRdr;
}

vector<vector<double>> calc_dRdr_2D(int ic, vector<vector<int>> IGA_element, vector<double> knot_vector_xi, vector<double> knot_vector_eta, vector<vector<double>> weight, \
int numOf_f_xi, int numOf_f_eta, int order_xi, int order_eta, double gauss_xi, double gauss_eta)
{
    vector<vector<double>> dRdr_2D(IGA_element[ic].size(), vector<double>(2));
    for(int i=0; i<IGA_element[ic].size(); i++){
        int index_eta, index_xi;
        int control_number = IGA_element[ic][i];
        if(control_number<=4) index_eta=0;
        else if(control_number>=5 && control_number<=8) index_eta=1;
        else if(control_number>=9 && control_number<=12) index_eta=2;
        if(control_number==1 || control_number==5 || control_number==9) index_xi=0;
        else if(control_number  ==2 || control_number==6 || control_number==10) index_xi=1;
        else if(control_number==3 || control_number==7 || control_number==11) index_xi=2;
        else if(control_number==4 || control_number==8 || control_number==12) index_xi=3;
        if(IGA_element[ic][i]-1==3){
            cout << index_xi << " " << index_eta << endl;
        }
        double dRdr_xi=calc_dRdr_2D_xi(knot_vector_xi, knot_vector_eta, weight, \
            numOf_f_xi, numOf_f_eta, order_xi, order_eta, index_xi, index_eta, gauss_xi, gauss_eta);
        double dRdr_eta=calc_dRdr_2D_eta(knot_vector_xi, knot_vector_eta, weight, \
            numOf_f_xi, numOf_f_eta, order_xi, order_eta, index_xi, index_eta, gauss_xi, gauss_eta);
        
        
        dRdr_2D[i][0] = dRdr_xi;
        dRdr_2D[i][1] = dRdr_eta;
        //if(IGA_element[ic][i]-1==3){
        //    cout << dRdr_xi << " " << dRdr_eta << endl;
        //}
    }
    
    return dRdr_2D;
}

void calc_dxdr(int ic, vector<vector<vector<double>>> node, vector<vector<int>> IGA_element, vector<vector<double>> &dxdr, vector<vector<double>> dRdr)
{
    for(int k=0;k<2;k++){
        for(int l=0;l<2;l++){
            dxdr[k][l] = 0e0;
            for(int m=0; m<IGA_element[ic].size(); m++){
                int control_number = IGA_element[ic][m];
                int eta_dir, xi_dir;
                int numOf_f_xi = 4;
                if(control_number%4==0){
                    eta_dir = control_number/4-1;
                }
                else{
                    eta_dir = control_number/4;
                }
                xi_dir = control_number - eta_dir * numOf_f_xi -1;
                dxdr[k][l] += dRdr[m][k] * node[eta_dir][xi_dir][l];
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
    dNdx.resize(IGA_element[ic].size());
    for(int i=0; i<dNdx.size(); i++){
        dNdx[i].resize(2);
    }
    for(int k=0; k<IGA_element[ic].size(); k++){
      for(int l=0; l<2; l++){
        dNdx[k][l] = 0.0;
        for(int p=0; p<2; p++){
          dNdx[k][l] += dNdr[k][p]*drdx[p][l];
        }
      }
    }
}

int main()
{
    int order_xi = 2;
    int order_eta = 2;

    vector<double> knot_vector_xi={0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};
    vector<double> knot_vector_eta={0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    vector<vector<double>> control_point = {
        {-1, 0},
        {-1, 0.41421356},
        {-0.41421356, 1},
        {0, 1},
        {-2.5, 0},
        {-2.5, 0.75},
        {-0.75, 2.5},
        {0, 2.5},
        {-4, 0},
        {-4, 4},
        {-4, 4},
        {0, 4}
    };
    vector<vector<double>> nurbs_weight = {

    };
    vector<vector<int>> element = {
        {0,1,2,4,5,6,8,9,10},
        {1,2,3,5,6,7,9,10,11}
    };
    vector<double> gauss = {-0.86113, -0.33998, 0.33998, 0.86113};
    vector<double> gauss_weight = {0.34785, 0.65214, 0.65214, 0.34785};

    vector<double> xi_coordinate = {0.0, 0.5, 1.0};
    vector<double> eta_coordinate = {0.0, 1.0};

    vector<vector<double>> K(4, vector<double>(4));

    for(int i=0; i<element.size(); i++){//elementループ
        for(int j=0; j<element[i].size(); j++){//control_pointのループ(row)
            for(int k=0; k<element[i].size(); j++){//control_pointのループ(column)
                for(int l=0; l<gauss.size(); l){
                    
                }
            }
        }
    }

    //vector<vector<vector<double>>> control_points;
    //vector<vector<double>> weight;
//
    //input_knot_vector_2D(knot_vector_xi, knot_vector_eta);
//
    //int order_xi=2, order_eta=2;
    //int numOf_f_xi = knot_vector_xi.size()-1-order_xi;
    //int numOf_f_eta = knot_vector_eta.size()-1-order_eta;
    //control_points.resize(numOf_f_eta);
    //weight.resize(numOf_f_eta);
    //for(int i=0; i<numOf_f_eta; i++){
    //    control_points[i].resize(numOf_f_xi);
    //    weight[i].resize(numOf_f_xi);
    //}
    //input_control_points_2D(control_points, numOf_f_xi, numOf_f_eta);
//
    //input_weight(weight, numOf_f_xi, numOf_f_eta);
//
    //cout << "numOfeta: " << numOf_f_eta << endl;
    //cout << "numOfxi: " << numOf_f_xi << endl;
    //ofstream ofs("output.dat");
    //for(double i=0.0; i<=1.0; i+=0.5){
    //   for(double j=0.0; j<=1.0; j+=0.5){
    //       double r_sum_X=0.0;
    //       double r_sum_Y=0.0;
    //       for(int index_eta=0; index_eta<numOf_f_eta; index_eta++){
    //           for(int index_xi=0; index_xi<numOf_f_xi; index_xi++){
    //               double R = calc_R_2D(knot_vector_xi,knot_vector_eta,weight,numOf_f_xi,numOf_f_eta,order_xi,order_eta,index_xi,index_eta,j,i);
    //               r_sum_X += R * control_points[index_eta][index_xi][0];
    //               r_sum_Y += R * control_points[index_eta][index_xi][1];
    //           }
    //       }
    //       ofs << r_sum_X << " " << r_sum_Y << endl;
    //   }
    //}

    //vector<vector<double>> K(control_points.size()*control_points[0].size(), vector<double>(control_points.size()*control_points[0].size(), 0.0));
    //vector<double> p(control_points.size());
//
    //vector<vector<int>> IGA_element = {{1,2,3,5,6,7,9,10,11},{2,3,4,6,7,8,10,11,12}};
//
    //vector<double> gauss(3,0.0);
    //vector<double> gauss_weight(3, 0.0);
    //gauss[0] = -sqrt(3.0/5.0); gauss[1] = 0.0; gauss[2] = sqrt(3.0/5.0); 
    //gauss_weight[0] = 0.55556; gauss_weight[1] = 0.88889; gauss_weight[2] = 0.55556;
    //vector<vector<double>> dxdr(2, vector<double>(2)), drdx(2, vector<double>(2));
    //vector<vector<double>> dNdx;;

    //for(int i=0; i<2; i++){ //for element
    //    for(int j=0; j<9; j++){
    //        for(int k=0; k<9; k++){ 
    //            for(int l=0; l<3; l++){//gauss
    //                for(int m=0; m<3; m++){//gauss      
    //                    vector<vector<double>> dRdr =  calc_dRdr_2D(i, IGA_element, knot_vector_xi, knot_vector_eta, weight, \
    //                        numOf_f_xi, numOf_f_eta, order_xi, order_eta, (gauss[m]*0.5+0.5)*0.5, (gauss[l]*1.0+1.0)*0.5);
    //                    
    //                    calc_dxdr(i, control_points, IGA_element, dxdr, dRdr);
    //                    
    //                    calc_inverse_matrix_2x2(dxdr, drdx);
    //                    
    //                    calc_dNdx(i, IGA_element, dNdx, dRdr, drdx);
    //                    
    //                    double detJ1 = fabs(dxdr[0][0] * dxdr[1][1]  - dxdr[1][0] * dxdr[0][1]);
    //                    double detJ2 = 1.0 / 8.0;
    //                    K[IGA_element[i][j]-1][IGA_element[i][k]-1] += (dNdx[j][0]*dNdx[k][0] + dNdx[j][1] * dNdx[k][1]) * detJ1 * detJ2 * gauss_weight[l];
    //                }
    //            }
    //        }
    //    }
    //}

    //for(int i=0; i<K.size(); i++){
    //    for(int j=0; j<K[i].size(); j++){
    //        cout << K[3][j] << " ";
    //    }
    //    cout << endl;
    //}
//
    //exit(1);
    //for(int i=0; i<K[i].size(); i++){
    //    K[11][i] = 0.0;
    //}
    //K[11][11] = 1.0;
//
    //vector<double> b={0,0,0,0,0,0,0,0,0,0,0,1.0};
    //vector<double> u(12,0.0);
    ////連立方程式ソルバー
    //int count=1;
    //while(1){
    //    vector<double> u_tmp(u.size(),0.0);
    //    for(int j=0; j<u_tmp.size(); j++){
    //        u_tmp[j]=b[j];
    //    }
//
    //    for(int j=0; j<u.size(); j++){
    //        for(int k=0; k<u.size(); k++){
    //            if(j!=k){
    //                u_tmp[j]-=K[j][k]*u[k];
    //            }
    //        }
    //        u_tmp[j]/=K[j][j];
    //    }
    //    double error=0;
    //    for(int j=0; j<u.size(); j++){
    //        error+=fabs(u[j]-u_tmp[j]);
    //    }
    //    cout << "iter : " << count << "  " << "error : " << error << endl;
    //    for(int j=0; j<u.size(); j++){
    //        u[j]=u_tmp[j];
    //    }
    //    u[0]=1.0;
    //    if(error<1e-5){
    //        break;
    //    }
    //    count++;
    //    break;
    //}
}