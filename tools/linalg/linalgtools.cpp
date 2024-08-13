#include<iostream>
#include<cmath>
#include<vector>
#include <omp.h>
#include<algorithm>
#include<set>
#include <eigen3/Eigen/Dense>
//#include <pybind11/pybind11.h>
//#include <pybind11/eigen.h>
#include<chrono>
#include<fstream>
#include<sstream>


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
           //     cout << omp_get_num_threads() << endl;
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
               // cout << omp_get_num_threads() << endl;
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
                cout << "Num Threads: " << omp_get_num_threads() << "\n" << endl;
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
            for(int j = 0; j < n; j++){
                if(pos_arr[j].size() == 0){
                    empty_row_index.push_back(j);
            //        cout << j << " empty" << endl;
                }
            }

       /*     for(int i = 0; i < pos_arr.size(); i++){
                cout << i << ":";
                for(int j = 0; j < pos_arr[i].size(); j++){
                    cout << pos_arr[i][j] << " ";
                }
                cout << endl;
            } */
            
            for(int k = 0; k < empty_row_index.size(); k++){
           //     cout << empty_row_index[k] << endl;
                pos_arr.erase(pos_arr.begin() + empty_row_index[k] - k);
            }
            n = n - empty_row_index.size();
        }

        vector<vector <int>> binary_combinations(int k){

            vector<vector<int> > combinations;
            int total_combs = pow(2,k);

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

        int get_bv_value(int i, vector<int> &free_vars, vector<int> fvs, vector<vector<int>> &bins){
            int fv_index = -1;
            int binary_value = 0;
            for (int k = 0; k < fvs.size(); k++){
                for (int j = 0; j < free_vars.size(); j++){
                    if(fvs[k] == free_vars[j]){
                        binary_value = (binary_value + bins[i][j])%2;
                        break;
                  }
                }
            }
            if(fvs.back() == m){
                binary_value = (binary_value + 1)%2;
            }
            return binary_value;
        }

    public:

        mod2system_solver(){
        }     

        void setSystem(int N, int M, Matrix<short,Dynamic,Dynamic> a){
            n = N;
            m = M;
         //   cout << N << endl;
            A = a;

        //    cout << A << endl;
//            cout << "System Initialised as: " << endl;
          //  print_matrix();
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
                //cout << "SWAP" << endl;
                swap_rows_positional(i);
               // print_positional();

               // cout << "Eliminate Below" << endl;
                eliminate_below_positional(i);
                //print_positional();

                //cout << "Eliminate Above" << endl;
                eliminate_above_positional(i);
                //print_positional();


//                cout << "After pruning empty rows" << endl;
                prune_empty_rows();
  //              print_positional();


            }
        }


        void gauss_jordan_elimination() {

            int prev_pivot = -1;
            int next_pivot;
            int final_row = n-1;
            

            for (int i=0; i<n; i++) {
                
                swap_rows(i);
             //  cout << "SWAP" << endl;      
               // print_matrix();

                next_pivot = find_pivot_col(i, prev_pivot);
               // cout << "Pivot coloumn: " << next_pivot << endl;

                //condition is satisfied when a full zero row is encountered
                if (next_pivot != -1){
                    eliminate_below(i);
                 //   cout << "Eliminate below" << endl;
                   // print_matrix();

                    eliminate_above(i,next_pivot);
                    //cout << "Eliminate above" << endl;
                //    print_matrix();
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
          //  print_matrix();

           // cout << "RREF finished" << endl;
        }

        

        Matrix<short,Dynamic,Dynamic> find_solutions(){
            
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
            }

            for(int i = 0; i < m; i++){
                if(count(basic_vars.begin(), basic_vars.end(), i) == 0){
                    free_vars.push_back(i);
                }
            }

            
            vector<vector<int>> bins;
            bins = binary_combinations(free_vars.size());
            Matrix<short,Dynamic,Dynamic> solutions;
            solutions.resize(bins.size(),m);
            
            #pragma omp parallel for
            for(int i = 0; i < bins.size(); i++){
                for(int bv = 0; bv < basic_vars.size(); bv++){
                    solutions(i,basic_vars[bv]) = get_bv_value(i,free_vars,pos_arr[bv],bins);
                }
                for(int fv = 0; fv < free_vars.size(); fv++){
                    solutions(i,free_vars[fv]) = bins[i][fv];
                }
            }
            return solutions;
        }
};

/*
namespace py = pybind11;
//constexpr auto byref = py::return_value_policy::reference_internal;

PYBIND11_MODULE(LinAlg, m) {
    m.doc() = "optional module docstring";

    py::class_<mod2system_solver>(m, "mod2system_Solver")
        .def(py::init<>())
        .def("setSystem",&mod2system_solver::setSystem,
        py::arg("N"),py::arg("M"),py::arg("a"))
        .def("print_positional",&mod2system_solver::print_positional)
        .def("gauss_jordan_elimination_positional",&mod2system_solver::gauss_jordan_elimination_positional)
        .def("find_solutions",&mod2system_solver::find_solutions,py::return_value_policy::take_ownership);
    ;
}*/

vector<vector<int>> read_traces_from_file(string file){

    ifstream f;
    vector<std::vector<int>> array;    /* vector of vector<int>  */

    f.open (file);   /* open file with filename as argument */
    if (! f.is_open()) {    /* validate file open for reading */
        std::cerr << "error: file open failed " << "'.\n";
    }
    else{
        string line, val;                  /* string for line & value */

        while (std::getline (f, line)) {        /* read each line */
            vector<int> v;                 /* row vector v */
            stringstream s (line);         /* stringstream line */
            while (getline (s, val, ','))       /* get each value (',' delimited) */
                v.push_back (stoi (val));  /* add to row vector */
            array.push_back (v);                /* add row vector to array */
        }
    }
    return array;
}

int main () {

    string files[] = {"system128.csv", "system256.csv", "system512.csv", "system1024.csv", "system2048.csv", "system4096.csv"};
    for(string file: files){

        vector<vector<int>> traces = read_traces_from_file(file);
        int n = traces.size();
        int m = traces[0].size()-1;
        Matrix<short,Dynamic,Dynamic> A(n,m+1);
        
        for(int i = 0; i < n; i++){
            for(int j = 0; j < m+1; j++){
                A(i,j) = traces[i][j];
            }
        }

        mod2system_solver lin_system;
        lin_system.setSystem(n,m,A);
        auto start = chrono::high_resolution_clock::now();
        lin_system.gauss_jordan_elimination();
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop-start);
        cout << duration.count() << endl;
        //Matrix<short,Dynamic,Dynamic> solutions;
        //lin_system.find_solutions(solutions);
        //cout << solutions << endl;
    }
    return 0;
}

/*
int main1(){
    int n = 11;
    int m = 10;

    int lines[11][11]{
        {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1},
        {0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0},
        {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0},
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0},
        {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1}
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

    Matrix<short,Dynamic,Dynamic> solutions;
    //lin_system.find_solutions(solutions);
    cout << solutions << endl;
    return 0;
}



int main2(){

    int n = 12;
    int m = 6;


    int lines[12][7] = {
	{1,1,0,0,0,0,0}, //swap first element to 1
	{0,0,1,1,0,0,0},
	{0,0,0,0,1,1,1},
	{1,0,0,0,0,0,0},
	{0,0,1,0,0,0,1},
	{0,0,1,0,0,0,1},
    {1,1,0,0,0,0,0},
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

    Matrix<short,Dynamic,Dynamic> solutions;
   // lin_system.find_solutions(solutions);
    cout << solutions << endl;
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

    Matrix<short,Dynamic,Dynamic> solutions;
    solutions = lin_system.find_solutions();
    cout << solutions << endl;
}
*/