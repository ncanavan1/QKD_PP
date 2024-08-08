#include <iostream>
#include <vector>
#include <omp.h>

using namespace std;

void gauss_jordan_elimination(vector<vector<double>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    for (int col = 0; col < cols; col++) {
        // Find the pivot row
        int pivot_row = col;
        for (int row = col + 1; row < rows; row++) {
            if (abs(matrix[row][col]) > abs(matrix[pivot_row][col])) {
                pivot_row = row;
            }
        }

        // Swap rows if necessary
        if (pivot_row != col) {
            swap(matrix[col], matrix[pivot_row]);
        }

        // Eliminate below the pivot
        #pragma omp parallel for
        for (int row = col + 1; row < rows; ++row) {
            double factor = matrix[row][col];
            for (int j = 0; j < cols; ++j) {
                matrix[row][j] -= factor * matrix[col][j];
            }
        }
    }

    // Backward elimination
    for (int col = cols - 1; col >= 0; --col) {
        #pragma omp parallel for
        for (int row = 0; row < col; ++row) {
            double factor = matrix[row][col];
            for (int j = 0; j < cols; ++j) {
                matrix[row][j] -= factor * matrix[col][j];
            }
        }
    }
}

int main() {
    vector<vector<double>> matrix = {
        {2, 1, -1, -3},
        {-3, -1, 2, -3},
        {-2, 1, 2, -3}
    };

    gauss_jordan_elimination(matrix);

    for (const auto& row : matrix) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    return 0;
}
