#include<iostream>
#include<cmath>
#include<vector>
#include <omp.h>
#include<algorithm>
#include<set>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class mod2system_solver{

    private:

        int n;
        int m;
        Matrix<short,Dynamic,Dynamic> A;
        vector< vector< int>> pos_arr;


        void swap_rows_positional(int i){
            int minPos = pos_arr[i][0];
            int minRow = i;
            for(int k=i; k<n; k++){
                if (pos_arr[k][0] < minPos){
                    minPos = pos_arr[k][0];
                    minRow = k;
                }
            }
            pos_arr[i].swap(pos_arr[minRow]);       
        }

        void swap_rows(int i){
            int maxEl = abs(A(i,i));
            int maxRow = i;
            for (int k=i+1; k<n; k++) {
                if (abs(A(k,i)) > maxEl) {
                    maxEl = abs(A(k,i));
                    maxRow = k;
                }
            }
            for (int k=i; k<m+1;k++) {
                int tmp = A(maxRow,k);
                A(maxRow,k) = A(i,k);
                A(i,k) = tmp;
            }
        }

        vector<int> find_difference(vector<int> v1, vector<int> v2){
            vector<int> difference;
            for(int i =0; i < v1.size(); i++){
                bool match_v1 = false;
                for(int j = 0; j < v2.size(); j++){
                    if(v1[i] == v2[j]){
                        match_v1 = true;
                    }
                }
                if(match_v1 == false){
                    difference.push_back(v1[i]);
                }
            }

            for(int j =0; j < v2.size(); j++){
                bool match_v2 = false;
                for(int i = 0; i < v1.size(); i++){
                    if(v1[i] == v2[j]){
                        match_v2 = true;
                    }
                }
                if(match_v2 == false){
                    difference.push_back(v2[j]);
                }
            }
            sort(difference.begin(), difference.end());
            return difference;
        }

        void eliminate_below_positional(int i){
            #pragma omp parallel for
            for(int k=i+1; k < n; k++){
                //if first position in rows matches, conduct elimination. Equal to 1 in lower row in pivot col
                if(pos_arr[i][0] == pos_arr[k][0]){
                    vector<int> tmp;
                    //finding the difference is eqivalent to operating an XOR on the explicit rows
                    tmp = find_difference(pos_arr[i],pos_arr[k]);
                    pos_arr[k] = tmp;
                }
            }
        }


        void eliminate_below(int i){

            #pragma omp parallel for
            for (int k=i+1; k<n; k++) {
                int c = A(k,i);
                if (c != 0){
                    for (int j=i; j<m+1; j++) {
                        if (i==j) {
                            A(k,j) = 0;
                        } else {
                            A(k,j) = (A(k,j) + c * A(i,j))%2;
                        }
                    }
                }
            }
        }


        void eliminate_above_positional(int i){
            #pragma omp parallel for
            for (int k=0; k < i; k++){
                for (int j = 0; j < pos_arr[k].size(); j++){
                    if (pos_arr[i][0] == pos_arr[k][j]){
                        vector<int> tmp;
                        tmp = find_difference(pos_arr[i],pos_arr[k]);
                        pos_arr[k] = tmp;
                    }
                }
            }
        }

        void eliminate_above(int i, int pivot){
            
            #pragma omp parallel for
            for (int k=0; k < i; k++){
                if (A(k,pivot) == 1){
                    for (int j = pivot; j <m+1; j++){
                        A(k,j) = (A(k,j) + A(i,j))%2; 
                    }    
                }
            }
        }
        

        int find_pivot_col(int row, int prev_pivot){

            int next_pivot = -1;
            for (int j = prev_pivot+1; j < m; j++){
                if (A(row,j) == 1){
                    next_pivot = j; 
                    break;
                }
            }
            return next_pivot;
        }

        void convert_to_positional(){
            for(int i = 0; i < n; i++){
                vector<int> row;
                for(int j=0; j< m+1; j++){
                    if (A(i,j) == 1){
                        row.push_back(j);
                    }
                }
                pos_arr.push_back(row);
            }
        }

        void prune_empty_rows(){
                //remove empty rows
                vector<int> empty_row_index;
                #pragma omp parallel for
                for(int j = 0; j < pos_arr.size(); j++){
                    if(pos_arr[j].size() == 0){
                        empty_row_index.push_back(j);
                    }
                }

                for(int k = 0; k < empty_row_index.size(); k++){
                    pos_arr.erase(pos_arr.begin() + empty_row_index[k] - k);
                }
                n = n - empty_row_index.size();
        }

    public:

        mod2system_solver(){
        }     

        void setSystem(int N, int M, Matrix<short,Dynamic,Dynamic> &a){
            n = N;
            m = M;
            A = a;
            cout << "System Initialised as: " << endl;
            print_matrix();
            convert_to_positional();
        }

        void print_positional(){
            for(int i=0; i<n; i++){
                for(int j=0;j<m+1;j++){
                    if (count(pos_arr[i].begin(), pos_arr[i].end(),j) == 1){
                        cout << "1" << "\t";
                    }
                    else{
                        cout << "0" << "\t";
                    }
                    if(j == m-1){
                        cout << "|";
                    }
                }
                cout << "\n";
            }
            cout << endl;
        }

        void print_matrix(){
            for (int i=0; i<n; i++){
                for(int j=0; j<m+1;j++){
                    cout << A(i,j) << "\t";
                    if (j == m-1){
                        cout << "|";
                    }
                }
                cout << "\n";
            }
            cout << endl;
        }   


        void gauss_jordan_elimination_positional(){

            for (int i=0; i<n; i++) {
                cout << "SWAP" << endl;
                swap_rows_positional(i);
                print_positional();

                cout << "Eliminate Below" << endl;
                eliminate_below_positional(i);
                print_positional();

                cout << "Eliminate Above" << endl;
                eliminate_above_positional(i);
                print_positional();


                cout << "After pruning empty rows" << endl;
                prune_empty_rows();
                print_positional();


            }
        }


        void gauss_jordan_elimination() {

            int prev_pivot = -1;
            int next_pivot;
            int final_row = n-1;
            

            for (int i=0; i<n; i++) {
                
                swap_rows(i);
                cout << "SWAP" << endl;      
                print_matrix();

                next_pivot = find_pivot_col(i, prev_pivot);
                cout << "Pivot coloumn: " << next_pivot << endl;

                //condition is satisfied when a full zero row is encountered
                if (next_pivot != -1){
                    eliminate_below(i);
                    cout << "Eliminate below" << endl;
                    print_matrix();

                    eliminate_above(i,next_pivot);
                    cout << "Eliminate above" << endl;
                    print_matrix();
                }
                //the swap condition should ensure that there are only zero rows below
                else{
                    final_row = i-1;
                    break;                    
                }
                prev_pivot = next_pivot;
            }

            Matrix<short, Dynamic, Dynamic> temp(final_row+1,m+1);
            temp = A(seq(0,final_row),seq(0,m));
            A = temp;
            n = final_row + 1;
            print_matrix();

            cout << "RREF finished" << endl;
        }

        vector<vector <int>> binary_combinations(int k){

            vector<vector<int> > combinations;
            int total_combs = pow(2,k);

            #pragma omp parallel for
            for(int i = 0; i < total_combs; i++){
                vector<int> comb;
                for(int j = k-1; j >= 0; j --){
                    if(i & (1 << j)){
                        comb.push_back(1);
                    }
                    else{
                        comb.push_back(0);
                    }
                }
                combinations.push_back(comb);
            }            
            return combinations;
        }

        Matrix<int,Dynamic,Dynamic> find_solutions(){
            
            //find free variables, number k
            //create solution space of 2^k x m
            //find all bimary combinations of free variables
            //solve for the resulting solutions

            vector<int> basic_vars;
            vector<int> free_vars;
            for(int i = 0; i < n; i++){
                basic_vars.push_back(pos_arr[i][0]);

                //if not a basic variable, then its a free variable
                //second condition stops considering b vector as in A matrix
                if(pos_arr.size() > 1){
                    for(int j = 1; j < pos_arr.size(); j ++){
                        if (pos_arr[i][j] != m && count(free_vars.begin(), free_vars.end(), pos_arr[i][j]) == 0){
                            free_vars.push_back(pos_arr[i][j]);
                        }
                    }
                }
            }
            
            Matrix<short,free_vars.size(),m> solutions;
            vector<vector<int>> bins;
            bins = binary_combinations(free_vars.size());
            //for each combination
            for(int i = 0; i < pow(2,free_vars.size()), i++){
                //set combination parameters
                for(int j = 0; j < free_vars.size(); j ++){
                    solutions[i][free_vars[j]] = bins[i][j];
                }
                
                //solve basic variables
                for(int j = 0; j < basic_vars.size(); j++){
                    int solved_value = 0;
                    for(int k=1; k < pos_arr[basic_vars[j]].size() ; k++){
                        if()
                    }
                }
            }



        }
        
};



int main1(){

    int n = 11;
    int m = 6;


    int lines[11][7] = {
	{1,1,0,0,0,0,0}, //swap first element to 1
	{0,0,1,1,0,0,0},
	{0,0,0,0,1,1,1},
	{1,0,0,0,0,0,0},
	{0,0,1,0,0,0,1},
	{0,0,1,0,0,0,1},
	{1,0,0,1,1,1,0},
	{0,1,1,0,0,0,1},
	{1,1,0,0,0,0,0},
	{0,0,1,1,0,0,0},
	{0,0,0,0,1,1,1}};

    Matrix<short,Dynamic,Dynamic> A(n,m+1);
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m+1; j++){
            A(i,j) = lines[i][j];
        }
    }

    mod2system_solver lin_system;
    lin_system.setSystem(n,m,A);
    lin_system.gauss_jordan_elimination_positional();
    //lin_system.find_solutions();
    return 0;
}



int main() {
    

    int n = 3;
    int m = 6;

    int lines[3][7] = {
        {1,1,0,1,0,0,1},
        {1,0,1,0,1,0,0},
        {0,1,1,0,0,1,0}
    };
    
    Matrix<short,Dynamic,Dynamic> A(n,m+1);
    
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m+1; j++){
            A(i,j) = lines[i][j];
        }
    }

    mod2system_solver lin_system;
    lin_system.setSystem(n,m,A);
    lin_system.gauss_jordan_elimination_positional();
    lin_system.find_solutions();

}