#include<iostream>
#include<vector>
#include<fstream>
#include<cmath>
#include<algorithm>
#include<map>

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

    //cout << N << endl;

    return N;
}

int return_insert_position(vector<double> old_knot_vector, double insert_knot)
{
    for(int i=0; i<old_knot_vector.size(); i++){
        if(fabs(old_knot_vector[i]-insert_knot)<0.00001){
            return i;
        }
    }
}

void insert_knot_to_current_knot_vector(vector<double> &old_knot_vector, double insert_knot, int insert_position)
{
    vector<double> tmp_knot_vector;
    int count = 1;
    for(int i=0; i<old_knot_vector.size(); i++){
        if(count==insert_position){
            tmp_knot_vector.push_back(insert_knot);
            tmp_knot_vector.push_back(old_knot_vector[i]);
            count++;
        }
        else{
            tmp_knot_vector.push_back(old_knot_vector[i]);
            count++;
        }
    }
    old_knot_vector = tmp_knot_vector;
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

vector<vector<double>> one_d_ectraction_operator(int numOfk, int order, vector<double> &knot_vector, vector<double> insert_knot)
{
    vector<vector<vector<double>>> operator_C;
    for(int i=1; i<=insert_knot.size(); i++){
        int insert_position = return_insert_position(knot_vector, insert_knot[i-1]);
        if(i!=1){
            if(insert_knot[i-1]==insert_knot[i-2]){
                insert_position++;
            }
        }
        
        vector<double> alpha;
        for(int j=0; j<numOfk+i; j++){
            //cout << numOfk+1 << " " << j << " " << insert_position+1 << " ";

            if(j+1<=insert_position-order && j+1>=1){
                alpha.push_back(1.0);
            }
            else if(j+1<=insert_position && j+1>=insert_position-order+1){
                alpha.push_back((insert_knot[i-1]-knot_vector[j])/(knot_vector[j+order]-knot_vector[j]));
            }
            else if(j+1>=insert_position+1){
                alpha.push_back(0.0);
            }
            else{
            }
        }
        vector<vector<double>> C(numOfk+i-1, vector<double>(numOfk+i));
        for(int j=0; j<numOfk+i-1; j++){
            C[j][j] = alpha[j];
            C[j][j+1] = 1.0-alpha[j+1];
        }

        operator_C.push_back(C);
        insert_knot_to_current_knot_vector(knot_vector, insert_knot[i-1], insert_position+1);
    }

    reverse(operator_C.begin(), operator_C.end());

    vector<vector<double>> final_operator = mat_product(Transposed_mat(operator_C[0]), Transposed_mat(operator_C[1]));
    
    for(int i=2; i<operator_C.size(); i++){
        final_operator = mat_product(final_operator, Transposed_mat(operator_C[i]));
    }

    final_operator = Transposed_mat(final_operator);

    return final_operator;
}

int main()
{
    //n=k+order+1
    ofstream ofs("result1.dat");
    int order = 3;
    vector<double> knot_vector = {0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0};
    int numOfk=knot_vector.size()-order-1;
    for(double i=0.0; i<=4.0; i+=0.01){
        for(int j=0; j<numOfk; j++){
            ofs << j << " " << i << " " << p_order_N(knot_vector, j+1, order, i) << endl;
        }
    }
    ofs.close();

    vector<double> insert_knot = {1.0, 1.0, 2.0, 2.0, 3.0, 3.0};
    

    vector<vector<double>> final_operator = one_d_ectraction_operator(numOfk, order, knot_vector, insert_knot);
    
    
    final_operator = Transposed_mat(final_operator);
    
    ofs.open("result3.dat");
    for(double i=0.0; i<=4.0; i+=0.001){
        for(int j=0; j<final_operator.size(); j++){
            double sum =0.0;
            for(int k=0; k<final_operator[0].size(); k++){
                sum += final_operator[j][k] * p_order_N(knot_vector, j+1, order, i);
            }
            ofs << j << " " << i << " " << sum << endl;
        }
    }
    ofs.close();
}