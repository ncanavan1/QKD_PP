#include<iostream>
#include<cmath>
#include<vector>

using namespace std;


void print(vector<vector<int> > A){
    int n = A.size();
    int m = A[0].size()-1;
    for (int i=0; i<n; i++){
        for(int j=0; j<m+1;j++){
            cout << A[i][j] << "\t";
            if (j == m-1){
                cout << "|";
            }
        }
        cout << "\n";
    }
    cout << endl;
}

vector< vector<int> > swap_rows(vector< vector<int> > A, int i, int n, int m){
    int maxEl = abs(A[i][i]);
    int maxRow = i;
    for (int k=i+1; k<n; k++) {
        if (abs(A[k][i]) > maxEl) {
            maxEl = abs(A[k][i]);
            maxRow = k;
        }
    }
    for (int k=i; k<m+1;k++) {
        int tmp = A[maxRow][k];
        A[maxRow][k] = A[i][k];
        A[i][k] = tmp;
    }
    return A;
}


vector< vector<int> > eliminate_below(vector< vector<int> > A, int i, int n, int m){
    
    for (int k=i+1; k<n; k++) {
        int c = A[k][i];
        for (int j=i; j<m+1; j++) {
            if (i==j) {
                A[k][j] = 0;
            } else {
                A[k][j] = (A[k][j] + c * A[i][j])%2;
            }
        }
    }
    return A;
}


vector< vector<int> > eliminate_above(vector< vector<int> > A, int i, int pivot, int n, int m){

    for (int k=0; k < i; k++){
        if (A[k][pivot] == 1){
            for (int j = pivot; j <m+1; j++){
                A[k][j] = (A[k][j] + A[i][j])%2; 
            }    
        }
    }
    return A;
}

int find_pivot_col(vector<int> row, int m, int prev_pivot){

    int next_pivot;
    for (int j = prev_pivot+1; j < m; j++){
        if (row[j] == 1){
            next_pivot = j; 
            break;
        }
    }
    return next_pivot;
}

vector<int> gauss(vector< vector<int> > A) {

    int n = A.size();
    int m = A[0].size()-1;
    int prev_pivot = -1;
    int next_pivot;

    for (int i=0; i<n; i++) {
        
        A = swap_rows(A,i,n,m);
        cout << "SWAP" << endl;      
        print(A);

        next_pivot = find_pivot_col(A[i], m, prev_pivot);
        cout << "Pivot coloumn: " << next_pivot << endl;

        A = eliminate_below(A,i,n,m);
        cout << "Eliminate below" << endl;
        print(A);

        A = eliminate_above(A,i,next_pivot,n,m);
        cout << "Eliminate above" << endl;
        print(A);

        prev_pivot = next_pivot;
    }

    

    vector<int> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = (A[i][n]/A[i][i])%2;
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= (A[k][i] * x[i])%2;
        }
    }
    return x;
}


int main() {
    int n = 3;
    int m = 6;

    vector<int> line(m+1,0); //m is matrix size, m+1 to include row solution
    vector< vector<int> > A(n,line);

    int lines[3][7] = {
        {1,1,0,1,0,0,1},
        {1,0,1,0,1,0,0},
        {0,1,1,0,0,1,0}
    };
    
    
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m+1; j++){
            A[i][j] = lines[i][j];
        }
    }


    //Print input
    print(A);
    //Calculate solution
    vector<int> x(n);
    x = gauss(A);
    cout << "Result:\t";
    for (int i=0; i<n; i++) {
        cout << x[i] << " ";
    }
    cout << endl;
}

int main_manual_inpt() {
    int n;
    cin >> n;
    vector<int> line(n+1,0);
    vector< vector<int> > A(n,line);
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            cin >> A[i][j];
            }
    }
    for (int i=0; i<n; i++) {
        cin >> A[i][n];
    }
    //Print input
    print(A);
    //Calculate solution
    vector<int> x(n);
    x = gauss(A);
    cout << "Result:\t";
    for (int i=0; i<n; i++) {
        cout << x[i] << " ";
    }
    cout << endl;
}