import numpy as np 
# Function to check if matrix is in REF 

def is_row_echelon_form(matrix): 
	if not matrix.any(): 
		return False

	rows = matrix.shape[0] 
	cols = matrix.shape[1] 
	prev_leading_col = -1

	for row in range(rows): 
		leading_col_found = False
		for col in range(cols): 
			if matrix[row, col] != 0: 
				if col <= prev_leading_col: 
					return False
				prev_leading_col = col 
				leading_col_found = True
				break
		if not leading_col_found and any(matrix[row, col] != 0 for col in range(cols)): 
			return False
	return True

def find_nonzero_row(matrix, pivot_row, col): 
	nrows = matrix.shape[0] 
	for row in range(pivot_row, nrows): 
		if matrix[row, col] != 0: 
			return row 
	return None

# Swapping rows so that we can have our non zero row on the top of the matrix 
def swap_rows(matrix, row1, row2, b, row_order): 
	matrix[[row1, row2]] = matrix[[row2, row1]] 
	b[[row1,row2]] = b[[row2,row1]]
	row_order[[row1,row2]] = row_order[[row2,row1]]

def make_pivot_one(matrix, pivot_row, col,b): 
	pivot_element = matrix[pivot_row, col] 
	matrix[pivot_row] //= pivot_element 
	b[pivot_row] //= pivot_element
	# print(pivot_element) 

def eliminate_below(matrix, pivot_row, col,b): 
	nrows = matrix.shape[0] 
	pivot_element = matrix[pivot_row, col] 
	for row in range(pivot_row + 1, nrows): 
		factor = matrix[row, col] 
		matrix[row] = (matrix[row] - factor * matrix[pivot_row])%2
		b[row] = (b[row] - factor * b[pivot_row])%2

def backward_reduce(matrix, pivot_row, b):
	
    col = np.where(matrix[pivot_row] == 1)
    if not np.size(col,axis=None):
        col = col[0][0]
        for row in range(0,pivot_row):
            if matrix[row,col] == 1:
                matrix[row] = (matrix[row] - matrix[pivot_row])%2
                b[row] = (b[row] - b[pivot_row])%2
		

# Implementing above functions 
def row_echelon_form(matrix,b): 
    nrows = matrix.shape[0] 
    ncols = matrix.shape[1] 
    pivot_row = 0
    row_order = np.arange(nrows)
    # this will run for number of column times. If matrix has 3 columns this loop will run for 3 times 
    for col in range(ncols): 
        nonzero_row = find_nonzero_row(matrix, pivot_row, col) 
        if nonzero_row is not None: 
            swap_rows(matrix, pivot_row, nonzero_row,b,row_order) 
            make_pivot_one(matrix, pivot_row, col,b) 
            eliminate_below(matrix, pivot_row, col,b) 
            if pivot_row != nrows-1:
                backward_reduce(matrix,pivot_row+1, b)
            pivot_row += 1
    return matrix, b, row_order

"""
matrix = np.array([[1,1,0,1,0,0],[1,0,1,0,1,0],[0,1,1,0,0,1]]) 
b = np.array([[1],[0],[0]])
print("Matrix Before Converting:") 
print(matrix)
print(b) 
print() 
matrix,b,row_order = row_echelon_form(matrix,b) 
print("After Converting to Row Echelon Form:") 
print(matrix)
print(b)
if is_row_echelon_form(matrix): 
	print("In REF") 
else: 
	print("Not in REF--------------->")"""
