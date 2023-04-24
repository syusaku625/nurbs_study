#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include<sstream>

using namespace std;

double return_N0(int index, double c, vector<double> knot_vector)
{
    if(knot_vector[index] >= c && c>=knot_vector[index-1.0] || fabs(c-knot_vector[index-1.0])<0.000001){
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
    if(fabs(knot_vector[index-1+order]-knot_vector[index-1])<0.00001){
        N1 = 0.0;
    }
    else{
        N1 = (c-knot_vector[index-1])/(knot_vector[index-1+order]-knot_vector[index-1]);
    }
    if(fabs(knot_vector[index-1+order+1]-knot_vector[index-1+1])<0.00001){
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

void calc_dxdr_3D(int ic, vector<vector<double>> node, vector<vector<int>> element, vector<vector<double>> &dxdr, vector<vector<double>> dRdr)
{
    for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
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
      for(int l=0; l<3; l++){
        dNdx[k][l] = 0.0;
        for(int p=0; p<3; p++){
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
        cout << err << endl;
        if (err < convergence){
            break;
        }
    }
}

void calc_dRdr_3D(vector<vector<double>> &dRdr, int ic, vector<vector<int>> element, vector<vector<int>> element_xi, \
vector<vector<int>> element_eta, vector<vector<int>> element_zeta, vector<double> knot_vector_xi, vector<double> knot_vector_eta, \
vector<double> knot_vector_zeta, vector<double> xi_coordinate, vector<double> eta_coordinate, vector<double> zeta_coordinate,\
double gauss_1, double gauss_2, double gauss_3, int order_xi, int order_eta, int order_zeta,\
int xi_loop_num, int eta_loop_num, int zeta_loop_num, vector<double> control_weight)
{   
    for(int dir=0; dir<3; dir++){
        for(int i=0; i<element_eta[eta_loop_num].size(); i++){//eta dir loop 3
            for(int j=0; j<element_xi[xi_loop_num].size(); j++){//xi dir loop 4
                for(int k=0; k<element_zeta[zeta_loop_num].size(); k++){//xi dir loop 4
                    //xi direction
                    double sum_W=0.0;
                    double sum_dWdr=0.0;
                    for(int l=0; l<element_eta[eta_loop_num].size(); l++){ //eta dir loop 3
                        for(int m=0; m<element_xi[xi_loop_num].size(); m++){ //xi dir loop 3
                            for(int n=0; n<element_zeta[zeta_loop_num].size(); n++){
                                if(dir==0){
                                    sum_dWdr += p_order_N(knot_vector_zeta,element_zeta[zeta_loop_num][n]+1,order_zeta,gauss_3)\
                                        *p_order_dNdr(knot_vector_xi,element_xi[xi_loop_num][m]+1,order_xi,gauss_1)\
                                        *p_order_N(knot_vector_eta,element_eta[eta_loop_num][l]+1,order_eta,gauss_2)*control_weight[element[ic][l*6+m*3+n]];

                                    sum_W += p_order_N(knot_vector_zeta,element_zeta[zeta_loop_num][n]+1,order_zeta,gauss_3)\
                                        *p_order_N(knot_vector_xi,element_xi[xi_loop_num][m]+1,order_xi,gauss_1)\
                                        * p_order_N(knot_vector_eta,element_eta[eta_loop_num][l]+1,order_eta,gauss_2)*control_weight[element[ic][l*6+m*3+n]];
                                    
                                }
                                if(dir==1){
                                    sum_dWdr +=  p_order_N(knot_vector_zeta,element_zeta[zeta_loop_num][n]+1,order_zeta,gauss_3)\
                                        *p_order_N(knot_vector_xi,element_xi[xi_loop_num][m]+1,order_xi,gauss_1)\
                                        *p_order_dNdr(knot_vector_eta,element_eta[eta_loop_num][l]+1,order_eta,gauss_2)*control_weight[element[ic][l*6+m*3+n]];

                                    sum_W += p_order_N(knot_vector_zeta,element_zeta[zeta_loop_num][n]+1,order_zeta,gauss_3)\
                                        *p_order_N(knot_vector_xi,element_xi[xi_loop_num][m]+1,order_xi,gauss_1)\
                                        * p_order_N(knot_vector_eta,element_eta[eta_loop_num][l]+1,order_eta,gauss_2)*control_weight[element[ic][l*6+m*3+n]];
                                }
                                if(dir==2){
                                    sum_dWdr +=  p_order_dNdr(knot_vector_zeta,element_zeta[zeta_loop_num][n]+1,order_zeta,gauss_3)\
                                        *p_order_N(knot_vector_xi,element_xi[xi_loop_num][m]+1,order_xi,gauss_1)\
                                        *p_order_N(knot_vector_eta,element_eta[eta_loop_num][l]+1,order_eta,gauss_2)*control_weight[element[ic][l*6+m*3+n]];

                                    sum_W += p_order_N(knot_vector_zeta,element_zeta[zeta_loop_num][n]+1,order_zeta,gauss_3)\
                                        *p_order_N(knot_vector_xi,element_xi[xi_loop_num][m]+1,order_xi,gauss_1)\
                                        * p_order_N(knot_vector_eta,element_eta[eta_loop_num][l]+1,order_eta,gauss_2)*control_weight[element[ic][l*6+m*3+n]];
                                }
                            }
                        }
                    }
                    double tmp_dRdr = 0.0;
                    if(dir==0){
                        tmp_dRdr = p_order_dNdr(knot_vector_xi,element_xi[xi_loop_num][j]+1,order_xi,gauss_1)\
                            *p_order_N(knot_vector_eta,element_eta[eta_loop_num][i]+1,order_eta,gauss_2)\
                            *p_order_N(knot_vector_zeta,element_zeta[zeta_loop_num][k]+1,order_zeta,gauss_3)*sum_W\
                            -p_order_N(knot_vector_xi,element_xi[xi_loop_num][j]+1,order_xi,gauss_1)\
                            *p_order_N(knot_vector_eta,element_eta[eta_loop_num][i]+1,order_eta,gauss_2)\
                            *p_order_N(knot_vector_zeta,element_zeta[zeta_loop_num][k]+1,order_zeta,gauss_3)*sum_dWdr;
                    }
                    if(dir==1){
                        tmp_dRdr = p_order_N(knot_vector_xi,element_xi[xi_loop_num][j]+1,order_xi,gauss_1)\
                            *p_order_dNdr(knot_vector_eta,element_eta[eta_loop_num][i]+1,order_eta,gauss_2)\
                            *p_order_N(knot_vector_zeta,element_zeta[zeta_loop_num][k]+1,order_zeta,gauss_3)*sum_W\
                            -p_order_N(knot_vector_xi,element_xi[xi_loop_num][j]+1,order_xi,gauss_1)\
                            *p_order_N(knot_vector_eta,element_eta[eta_loop_num][i]+1,order_eta,gauss_2)\
                            *p_order_N(knot_vector_zeta,element_zeta[zeta_loop_num][k]+1,order_zeta,gauss_3)*sum_dWdr;
                    }
                    if(dir==2){
                        tmp_dRdr = p_order_N(knot_vector_xi,element_xi[xi_loop_num][j]+1,order_xi,gauss_1)\
                            *p_order_N(knot_vector_eta,element_eta[eta_loop_num][i]+1,order_eta,gauss_2)\
                            *p_order_dNdr(knot_vector_zeta,element_zeta[zeta_loop_num][k]+1,order_zeta,gauss_3)*sum_W\
                            -p_order_N(knot_vector_xi,element_xi[xi_loop_num][j]+1,order_xi,gauss_1)\
                            *p_order_N(knot_vector_eta,element_eta[eta_loop_num][i]+1,order_eta,gauss_2)\
                            *p_order_N(knot_vector_zeta,element_zeta[zeta_loop_num][k]+1,order_zeta,gauss_3)*sum_dWdr;
                    }

                    tmp_dRdr/=(sum_W*sum_W);
                    dRdr[i*6+j*3+k][dir] = tmp_dRdr*control_weight[element[ic][i*6+j*3+k]];
                }
            }
        }
    }
}

void input_control(string filename, vector<vector<double>> &control_point)
{
    ifstream ifs(filename);
    if(!ifs){
        cout << "Error: control point file not opened." << endl;
        exit(1);
    }
    string str;
    while(getline(ifs,str)){
        istringstream ss(str);
        vector<double> tmp_control;
        for(int i=0; i<3; i++){
            getline(ss, str, ' ');
            tmp_control.push_back(stod(str));
        }
        control_point.push_back(tmp_control);
    }
}

void input_element(string filename, vector<vector<int>> &element)
{
    ifstream ifs(filename);
    if(!ifs){
        cout << "Error: element file not opened." << endl;
        exit(1);
    }
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

void input_knot_span_3D(string knot_span_xi_file, string knot_span_eta_file, string knot_span_zeta_file,\
vector<vector<int>> &element_xi, vector<vector<int>> &element_eta, vector<vector<int>> &element_zeta)
{
    string str;
    
    ifstream ifs(knot_span_xi_file);
    if(!ifs){
        cout << "Error: knot span xi file not opened." << endl;
        exit(1);
    }
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
    if(!ifs){
        cout << "Error: knot span eta file not opened." << endl;
        exit(1);
    }
    while(getline(ifs,str)){
        vector<int> tmp_knot_span;
        istringstream ss(str);
        while(getline(ss, str, ' ')){
            tmp_knot_span.push_back(stoi(str));
        }
        element_eta.push_back(tmp_knot_span);
    }
    ifs.close();

    ifs.open(knot_span_zeta_file);
    if(!ifs){
        cout << "Error: knot span zeta file not opened." << endl;
        exit(1);
    }
    while(getline(ifs,str)){
        vector<int> tmp_knot_span;
        istringstream ss(str);
        while(getline(ss, str, ' ')){
            tmp_knot_span.push_back(stoi(str));
        }
        element_zeta.push_back(tmp_knot_span);
    }
    ifs.close();
}

void input_parametric_coordinate_3D(string xi_coordinate_file, string eta_coordinate_file, string zeta_coordinate_file,\
vector<double> &xi_coordinate, vector<double> &eta_coordinate, vector<double> &zeta_coordinate)
{
    string str;
    ifstream ifs(xi_coordinate_file);
    if(!ifs){
        cout << "Error: xi coordinate file not opened." << endl;
        exit(1);
    }
    while(getline(ifs, str)){
        xi_coordinate.push_back(stod(str));
    }
    ifs.close();
    ifs.open(eta_coordinate_file);
    if(!ifs){
        cout << "Error: eta coordinate file not opened." << endl;
        exit(1);
    }
    while (getline(ifs,str)){
        eta_coordinate.push_back(stod(str));
    }
    ifs.close();
    ifs.open(zeta_coordinate_file);
    if(!ifs){
        cout << "Error: eta coordinate file not opened." << endl;
        exit(1);
    }
    while (getline(ifs,str)){
        zeta_coordinate.push_back(stod(str));
    }
    ifs.close();
}

void input_knot_connectivity_3D(string knot_connectivity_file, vector<int> &knot_connectivity_xi, vector<int> &knot_connectivity_eta, vector<int> &knot_connectivity_zeta)
{
    ifstream ifs(knot_connectivity_file);
    if(!ifs){
        cout << "Error: knot connectivity file not opened." << endl;
        exit(1);
    }
    string str;
    while(getline(ifs,str)){
        istringstream ss(str);
        getline(ss, str, ' ');
        knot_connectivity_xi.push_back(stoi(str));
        getline(ss, str, ' ');
        knot_connectivity_eta.push_back(stoi(str));
        getline(ss, str, ' ');
        knot_connectivity_zeta.push_back(stoi(str));
    }
    ifs.close();
}

void input_knot_vector_3D(string knot_vector_xi_file,string knot_vector_eta_file,string knot_vector_zeta_file, \
vector<double> &knot_vector_xi, vector<double> &knot_vector_eta, vector<double> &knot_vector_zeta)
{
    ifstream ifs(knot_vector_xi_file);
    if(!ifs){
        cout << "Error: knot vector xi file not opened." << endl;
        exit(1);
    }
    string str;
    while(getline(ifs,str)){
        knot_vector_xi.push_back(stod(str));
    }
    ifs.close();
    
    ifs.open(knot_vector_eta_file);
    if(!ifs){
        cout << "Error: knot vector eta file not opened." << endl;
        exit(1);
    }
    while(getline(ifs,str)){
        knot_vector_eta.push_back(stod(str));
    }
    ifs.close();

    ifs.open(knot_vector_zeta_file);
    if(!ifs){
        cout << "Error: knot vector zeta file not opened." << endl;
        exit(1);
    }
    while(getline(ifs,str)){
        knot_vector_zeta.push_back(stod(str));
    }
    ifs.close();
}

void input_control_weight(string control_weight_file, vector<double> &control_weight)
{
    ifstream ifs(control_weight_file);
    if(!ifs){
        cout << "Error: control weight file not opened." << endl;
        exit(1);
    }
    string str;
    while(getline(ifs,str)){
        control_weight.push_back(stod(str));
    }
    ifs.close();
}

void export_vtu(const std::string &file, vector<vector<int>> element, vector<vector<double>> x, vector<vector<double>> displacement)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", x.size(), element.size());
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * x.size() * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    for (int j = 0; j < element[i].size(); j++) fprintf(fp, "%d ", element[i][j]);
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < element.size(); i++)
  {
    num += element[i].size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++)
    fprintf(fp, "%d\n", 12);
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"displacement[m/s]\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * x.size() * 3;

  //fprintf(fp, "<DataArray type=\"Float64\" Name=\"pressure[Pa]\" NumberOfComponents=\"1\" format=\"appended\" offset=\"%d\"/>\n",offset);
  //offset += sizeof(int) + sizeof(double) * x.size();
  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[x.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < x.size(); ic++){
    for(int j=0;j<3;j++){
      data_d[num] = x[ic][j];
      num++;
    }
  }
  size=sizeof(double)*x.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  num=0;
  for (int ic = 0; ic < x.size(); ic++){
      data_d[num]   = displacement[ic][0];
      data_d[num+1] = displacement[ic][1];
      data_d[num+2] = displacement[ic][2];
      num=num+3;
  }
  size=sizeof(double)*x.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);
//
  //num=0;
  //for (int ic = 0; ic < numOfNode; ic++){
  //    data_d[num]   = p(ic);
  //    num++;
  //}
  //size=sizeof(double)*x.size();
  //ofs.write((char *)&size, sizeof(size));
  //ofs.write((char *)data_d, size);
//
  delete data_d;
  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

void export_vtu_preprpcess(vector<vector<double>> x, vector<vector<double>> displacement, int interpolate_element_xi, int interpolate_element_eta, int interpolate_element_zeta)
{
    int line_node_zeta = interpolate_element_zeta+1;
    int line_node_eta = interpolate_element_eta+1;
    int line_node_xi = interpolate_element_xi+1;
    
    vector<vector<int>> element;
    //for(int i = 0; i<interpolate_element_zeta; i++){
    //    for(int j=0; j<interpolate_element_eta; j++){
    //        for(int k=0; k<interpolate_element_xi; k++){
    //            vector<int> tmp_element;
    //            tmp_element.push_back((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i);
    //            
    //            tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+(line_node_zeta*line_node_xi)+(line_node_zeta));
    //            tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+(line_node_zeta*line_node_xi));
    //            tmp_element.push_back((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i+1);
    //            tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+(line_node_zeta)+1);
    //            tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+(line_node_zeta*line_node_xi)+(line_node_zeta)+1);
    //            tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+(line_node_zeta*line_node_xi)+1);
    //            element.push_back(tmp_element);
    //        }
    //    }
    //}
    for(int i = 0; i<interpolate_element_zeta; i++){
        for(int j=0; j<interpolate_element_eta; j++){
            for(int k=0; k<interpolate_element_xi; k++){
                vector<int> tmp_element;
                tmp_element.push_back((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i);
                tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+(line_node_zeta));
                tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+(line_node_zeta)+(line_node_zeta*line_node_xi));
                tmp_element.push_back((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i+(line_node_zeta*line_node_xi));
                tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+1);
                tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+(line_node_zeta)+1);
                tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+(line_node_zeta)+(line_node_zeta*line_node_xi)+1);
                tmp_element.push_back(((line_node_zeta*line_node_xi)*j+line_node_zeta*k+i)+(line_node_zeta*line_node_xi)+1);
                element.push_back(tmp_element);
            }
        }
    }
    export_vtu("test.vtu",element,x,displacement);
}

void return_D_marix_3D(vector<vector<double>> &D)
{
    double Young_modulus = 2e9;
    double Poisson_ratio = 0.3;

    double Ex = (1.0 - Poisson_ratio) * Young_modulus / ((1.0 + Poisson_ratio) * (1.0 - 2.0 * Poisson_ratio));
    double Ey = Poisson_ratio * Young_modulus / ((1.0 + Poisson_ratio) * (1.0 - 2.0 * Poisson_ratio));
    double G = Young_modulus / (2.0 * (1.0 + Poisson_ratio));
    
    D[0][0] = Ex; D[0][1] = Ey; D[0][2] = Ey;
    D[1][0] = Ey; D[1][1] = Ex; D[1][2] = Ey;
    D[2][0] = Ey; D[2][1] = Ey; D[2][2] = Ex;
    D[3][3] = G; D[4][4] = G; D[5][5] = G;
}

double calcDeterminant_3x3(vector<vector<double>> a)
{
  double det  = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2]
              - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
  return det;
}

void calcInverseMatrix_3x3(vector<vector<double>> &inv_a,vector<vector<double>> a)
{
  double det;

  det = calcDeterminant_3x3(a);

  inv_a[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  inv_a[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
  inv_a[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
  inv_a[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  inv_a[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
  inv_a[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
  inv_a[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  inv_a[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
  inv_a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) inv_a[i][j] = inv_a[i][j] / det;
  }
}

vector<vector<double>> calc_B_matrix_3D(vector<int> element, vector<vector<double>> dNdx)
{
    vector<vector<double>> B(6, vector<double>(element.size()*3, 0.0));
    for (int i = 0; i < element.size(); i++){
        B[0][i*3] = dNdx[i][0];
        B[1][i*3+1] = dNdx[i][1];
        B[2][i*3+2] = dNdx[i][2];
        B[3][i*3+1] = dNdx[i][2];
        B[3][i*3+2] = dNdx[i][1];
        B[4][i*3] = dNdx[i][2];
        B[4][i*3+2] = dNdx[i][0];
        B[5][i*3] = dNdx[i][1];
        B[5][i*3+1] = dNdx[i][0];
    }
    return B;
}

vector<vector<double>> Transposed_mat_3D(vector<vector<double>> C)
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

int main()
{
    int order_xi = 1;
    int order_eta = 1;
    int order_zeta = 2;
    string input_dir = "input_file2";
    vector<double> knot_vector_xi, knot_vector_eta, knot_vector_zeta;
    string knot_vector_xi_file = input_dir + "/" + "knot_vector_xi.dat";
    string knot_vector_eta_file = input_dir + "/" + "knot_vector_eta.dat";
    string knot_vector_zeta_file = input_dir + "/" + "knot_vector_zeta.dat";

    string control_file = input_dir + "/" + "control.dat";
    string element_file = input_dir + "/" + "element.dat";

    string knot_span_xi_file = input_dir + "/" + "knot_span_xi.dat";
    string knot_span_eta_file = input_dir + "/" + "knot_span_eta.dat";
    string knot_span_zeta_file = input_dir + "/" + "knot_span_zeta.dat";

    string knot_connectivity_file = input_dir + "/" + "knot_connectivity.dat";

    string xi_coordinate_file = input_dir + "/" + "xi_coordinate.dat";
    string eta_coordinate_file = input_dir + "/" + "eta_coordinate.dat";
    string zeta_coordinate_file = input_dir + "/" + "zeta_coordinate.dat";

    string control_weight_file = input_dir + "/" + "control_weight.dat";

    input_knot_vector_3D(knot_vector_xi_file,knot_vector_eta_file,knot_vector_zeta_file,knot_vector_xi,knot_vector_eta,knot_vector_zeta);
    vector<vector<double>> control_point;
    input_control(control_file,control_point);
    vector<vector<int>> element;
    input_element(element_file,element);
    vector<vector<int>> element_xi, element_eta, element_zeta;
    input_knot_span_3D(knot_span_xi_file,knot_span_eta_file,knot_span_zeta_file,element_xi,element_eta,element_zeta);
    vector<double> xi_coordinate, eta_coordinate, zeta_coordinate;
    input_parametric_coordinate_3D(xi_coordinate_file,eta_coordinate_file,zeta_coordinate_file,xi_coordinate,eta_coordinate,zeta_coordinate);
    vector<int> knot_connectivity_xi, knot_connectivity_eta, knot_connectivity_zeta;
    input_knot_connectivity_3D(knot_connectivity_file, knot_connectivity_xi, knot_connectivity_eta, knot_connectivity_zeta);
    vector<double> control_weight;
    input_control_weight(control_weight_file, control_weight);
    vector<double> gauss = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
    vector<double> gauss_weight = {5.0/9.0, 5.0/9.0, 5.0/9.0};
    vector<vector<double>> K(control_point.size()*3, vector<double>(control_point.size()*3));
    vector<vector<double>> D(6, vector<double>(6, 0.0));
    return_D_marix_3D(D);

    cout << control_weight.size() << endl;

    for(int i=0; i<element.size(); i++){//elementループ
        cout << i << endl;
        int eta_loop_num = knot_connectivity_eta[i];
        int xi_loop_num = knot_connectivity_xi[i];
        int zeta_loop_num = knot_connectivity_zeta[i];
        for(int j=0; j<element[i].size(); j++){//control_pointのループ(row)
            for(int k=0; k<element[i].size(); k++){//control_pointのループ(column)
                for(int l=0; l<gauss.size(); l++){//gauss loop 1
                    for(int m=0; m<gauss.size(); m++){//gauss loop 2
                        for(int n=0; n<gauss.size(); n++){//gauss loop 2
                            double gauss_eta = 0.5*((eta_coordinate[eta_loop_num+1]-eta_coordinate[eta_loop_num])*\
                                gauss[l]+(eta_coordinate[eta_loop_num]+eta_coordinate[eta_loop_num+1]));
                            double gauss_xi = 0.5*((xi_coordinate[xi_loop_num+1]-xi_coordinate[xi_loop_num])*\
                                gauss[m]+(xi_coordinate[xi_loop_num]+xi_coordinate[xi_loop_num+1]));
                            double gauss_zeta = 0.5*((zeta_coordinate[zeta_loop_num+1]-zeta_coordinate[zeta_loop_num])*\
                                gauss[n]+(zeta_coordinate[zeta_loop_num]+zeta_coordinate[zeta_loop_num+1]));
                            vector<vector<double>> dxdr(3, vector<double>(3));
                            vector<vector<double>> drdx(3, vector<double>(3));
                            vector<vector<double>> dNdx(element[i].size(), vector<double>(3));
                            vector<vector<double>> dRdr(element[i].size(), vector<double>(3));
                            calc_dRdr_3D(dRdr, i, element, element_xi, element_eta, element_zeta, knot_vector_xi, knot_vector_eta, knot_vector_zeta, xi_coordinate, \
                                eta_coordinate, zeta_coordinate, gauss_xi, gauss_eta, gauss_zeta, order_xi, order_eta, order_zeta, xi_loop_num, eta_loop_num, \
                                zeta_loop_num, control_weight);
                            
                            calc_dxdr_3D(i, control_point, element, dxdr, dRdr);
                            calcInverseMatrix_3x3(drdx,dxdr);
                            double detJ1 = calcDeterminant_3x3(dxdr);
                            calc_dNdx(i, element, dNdx, dRdr, drdx);
                            double detJ2 = 0.5*(eta_coordinate[eta_loop_num+1]-eta_coordinate[eta_loop_num])\
                            *0.5*(xi_coordinate[xi_loop_num+1]-xi_coordinate[xi_loop_num])\
                            *0.5*(zeta_coordinate[zeta_loop_num+1]-zeta_coordinate[zeta_loop_num]);

                            vector<vector<double>> B = calc_B_matrix_3D(element[i], dNdx);
                            vector<vector<double>> B_t = Transposed_mat_3D(B);
                            vector<vector<double>> K_e = mat_product(mat_product(B_t, D), B);
                            for(int o=0; o<3; o++){
                                for(int p=0; p<3; p++){
                                    K[element[i][j]*3+p][element[i][k]*3+o] += \
                                    K_e[j*3+p][k*3+o]*detJ1*detJ2*gauss_weight[l]*gauss_weight[m]*gauss_weight[n];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    vector<double> f(control_point.size()*3, 0.0);    

    vector<double> u(control_point.size()*3, 0.0);

    //vector<int> boundary_wall = {0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68};
    vector<int> boundary_wall = {0,4,8,12};
    vector<int> boundary_f = {3,7,11,15};

    for(int i=0; i<boundary_f.size(); i++){
        f[boundary_f[i]*3+1] = -3000000;
    }


    for(int i=0; i<boundary_wall.size(); i++){
        for(int j=0; j<K[boundary_wall[i]].size(); j++){
            K[boundary_wall[i]*3][j] = 0.0;
            K[boundary_wall[i]*3+1][j] = 0.0;
            K[boundary_wall[i]*3+2][j] = 0.0;
        }
        K[boundary_wall[i]*3][boundary_wall[i]*3] = 1.0;
        K[boundary_wall[i]*3+1][boundary_wall[i]*3+1] = 1.0;
        K[boundary_wall[i]*3+2][boundary_wall[i]*3+2] = 1.0;
    }

    double convegence=0.00001;
    Jacobi_method(u, f, K, convegence);
    for(int i=0; i<u.size(); i++){
        cout << u[i] << endl;
    }

    //可視化プロセス
    vector<vector<double>> x;
    vector<vector<double>> displacement;

    int interpolate_element_xi = 16;
    int interpolate_element_eta = 16;
    int interpolate_element_zeta = 64;

    double delta_eta =1.0/(double)interpolate_element_eta;
    double delta_xi =1.0/(double)interpolate_element_xi;
    double delta_zeta =1.0/(double)interpolate_element_zeta;

    //cout << "check3" << endl;

    for(double i=0; i<=1.0; i+=delta_eta){//eta
        //cout << i << endl;
        for(double j=0.0; j<=1.0; j+=delta_xi){//xi
            for(double k=0.0; k<=1.0; k+=delta_zeta){//zeta
                double r_sum_ux=0.0;
                double r_sum_uy=0.0;
                double r_sum_uz=0.0;
                double r_sum_X=0.0;
                double r_sum_Y=0.0;
                double r_sum_Z=0.0;
                for(int l=0; l<2; l++){//eta
                    for(int m=0; m<2; m++){//xi
                        for(int n=0; n<4; n++){//zeta
                            double lower_sum=0.0;
                            for(int o=0; o<2; o++){//eta
                                for(int p=0; p<2; p++){//xi
                                    for(int q=0; q<4; q++){//zeta
                                        lower_sum += p_order_N(knot_vector_zeta, q+1, order_zeta, k)*p_order_N(knot_vector_xi, p+1, order_xi, j)\
                                        *p_order_N(knot_vector_eta, o+1, order_eta, i)*control_weight[o*8+p*4+q];
                                    }
                                }
                            }
                            double tmp_R = p_order_N(knot_vector_zeta, n+1, order_zeta, k)*p_order_N(knot_vector_xi, m+1, order_xi, j)\
                                *p_order_N(knot_vector_eta, l+1, order_eta, i)*control_weight[l*8+m*4+n];
                                tmp_R /= lower_sum;
                            r_sum_ux += tmp_R * u[(l*8+m*4+n)*3];
                            r_sum_uy += tmp_R * u[(l*8+m*4+n)*3+1];
                            r_sum_uy += tmp_R * u[(l*8+m*4+n)*3+2];
                            r_sum_X += tmp_R * (control_point[l*8+m*4+n][0]+u[(l*8+m*4+n)*3]);
                            r_sum_Y += tmp_R * (control_point[l*8+m*4+n][1]+u[(l*8+m*4+n)*3+1]);
                            r_sum_Z += tmp_R * (control_point[l*8+m*4+n][2]+u[(l*8+m*4+n)*3+2]);
                            r_sum_X += tmp_R * (control_point[l*8+m*4+n][0]);
                            r_sum_Y += tmp_R * (control_point[l*8+m*4+n][1]);
                            r_sum_Z += tmp_R * (control_point[l*8+m*4+n][2]);
                        }
                    }
                }
                vector<double> tmp_x;
                vector<double> tmp_displacement;

                tmp_x.push_back(r_sum_X);
                tmp_x.push_back(r_sum_Y);
                tmp_x.push_back(r_sum_Z);   
                
                tmp_displacement.push_back(r_sum_ux);
                tmp_displacement.push_back(r_sum_uy);
                tmp_displacement.push_back(r_sum_uz);
//
                x.push_back(tmp_x);
                displacement.push_back(tmp_displacement);
                //]displacement.resize(x.size());
                //for(int i=0; i<x.size(); i++){
                //    displacement[i].resize(3);
                //}
            }
        }
    }
    
    export_vtu_preprpcess(x, displacement, interpolate_element_xi, interpolate_element_eta, interpolate_element_zeta);

}
