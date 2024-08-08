#include<iostream>
#include<cmath>
#include<vector>
#include <omp.h>
#include <eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class mod2system_solver{

    private:

        int n;
        int m;
        Matrix<short,Dynamic,Dynamic> A;

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


        void cut_empty_rows(){

            int last_nonZero_row = -1;
            for(int i=n-1; i >= 0; i--){
                bool isZeroRow = true;
                for(int j = 0; j < m+1; j++){
                    if (A(i,j) == 1){
                        isZeroRow = false;
                        last_nonZero_row = i;
                        break;
                    }
                }
                if (isZeroRow == false){
                    break;
                }
            }

            Matrix<short, Dynamic, Dynamic> temp(last_nonZero_row+1,m+1);
            temp = A(seq(0,last_nonZero_row),seq(0,m));
            A = temp;
            n = last_nonZero_row;

            print_matrix();
            
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

        void gauss_jordan_elimination() {

            int prev_pivot = -1;
            int next_pivot;
            int final_row = -1;
            

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

        Matrix<int,Dynamic,Dynamic> find_solutions(){
            
        }
        
};



int main(){

    int n = 11;
    int m = 6;


    int lines[11][7] = {
	{1,1,0,0,0,0,0},
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
    lin_system.gauss_jordan_elimination();

}



int main1() {
    

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
    lin_system.gauss_jordan_elimination();
}