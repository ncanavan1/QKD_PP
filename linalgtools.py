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

def eliminate_below(matrix, pivot_row, col, b): 
	nrows = matrix.shape[0] 
	pivot_element = matrix[pivot_row, col] 
	for row in range(pivot_row + 1, nrows): 
		factor = matrix[row, col] 
		matrix[row] = (matrix[row] - factor * matrix[pivot_row])%2
		b[row] = (b[row] - factor * b[pivot_row])%2

def eliminate_above(matrix,pivot_row,col,b):
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
			eliminate_above(matrix,pivot_row,col,b)
			pivot_row += 1

	first_empty = nrows
	for row in range(nrows):
		if (matrix[row] == np.zeros(ncols)).all():
			first_empty = row
			break
	matrix = matrix[:first_empty,:]
	b = b[:first_empty]

	return matrix, b, row_order


def find_solution(matrix,b):

	if is_row_echelon_form(matrix):
		basic_variables = []
		free_variables = []

		nrows = matrix.shape[0]
		ncols = matrix.shape[1]

		start_col = 0
		for row in range(nrows):
			for col in range(start_col,ncols):
				if matrix[row,col] == 1:
					if np.count_nonzero(matrix[:,col]) == 1 and np.count_nonzero(matrix[row,:col]) == 0:
						basic_variables.append(col)
						start_col = col
						break
		
		for i in range(ncols):
			if not basic_variables.__contains__(i):
				free_variables.append(i)

		sol_list = []
		for i in range(nrows):
			eqn = []
			eqn.append(b[i][0])
			for col in range(basic_variables[i]+1,ncols):
				if matrix[i,col] == 1:
					eqn.append(col)
			sol_list.append(eqn)

		k=7

	
	else:
		print("Invalid Matrix Form")
	



matrix = np.array([[1,1,0,1,0,0],[1,0,1,0,1,0],[0,1,1,0,0,1]]) 
b = np.array([[1],[0],[0]])

"""
matrix = np.array([
	[1,1,0,0,0,0],
	[0,0,1,1,0,0],
	[0,0,0,0,1,1],
	[1,0,0,0,0,0],
	[0,0,1,0,0,0],
	[0,0,1,0,0,0],
	[1,0,0,1,1,1],
	[0,1,1,0,0,0],
	[1,1,0,0,0,0],
	[0,0,1,1,0,0],
	[0,0,0,0,1,1]]
)
b = np.array([[0],[0],[1],[0],[1],[1],[0],[1],[0],[0],[1]])
"""

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
	print("Not in REF--------------->")



#matrix = np.array([[1,0,1,0,0,0],[0,1,1,0,0,1],[0,0,0,1,1,1]]) 
#b = np.array([[1],[0],[0]])

find_solution(matrix,b)
