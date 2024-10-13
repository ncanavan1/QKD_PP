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


def gen_binary_comb(n):
	bin_strings = []
	def get_bin(n,bs=''):
		if len(bs) == n:
			bin_strings.append(bs)
		else:
			get_bin(n, bs + '0')
			get_bin(n, bs + '1')
	get_bin(n)
	for i in range(len(bin_strings)):
		bs = []
		for d in bin_strings[i]:
			bs.append(int(d))
		bin_strings[i] = bs
	return bin_strings


def read_off_basic_variables(matrix,b):
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

		eqn_list = []
		for i in range(nrows):
			eqn = []
			for col in range(basic_variables[i],ncols):
				if matrix[i,col] == 1:
					eqn.append(col)
			eqn_list.append(eqn)

		x = np.ones(matrix.shape[1])*-1
		for i in range(len(eqn_list)):
			if len(eqn_list[i]) == 1:
				x[int(eqn_list[i][0])] = b[i]
		return x



		




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

		eqn_list = []
		for i in range(nrows):
			eqn = []
			for col in range(basic_variables[i]+1,ncols):
				if matrix[i,col] == 1:
					eqn.append(col)
			eqn_list.append(eqn)
		
		
		free_solutions = gen_binary_comb(len(free_variables))
		free_variables_space = np.asarray([free_variables]+free_solutions)

		basic_solution_space = []
		for i in range(1,free_variables_space.shape[0]):
			basic_solution = []
			for j in range(len(eqn_list)):
				bvar_val = int(b[j])
				for fv in eqn_list[j]:
					fv_index = np.where(free_variables_space[0] == fv)[0][0] ##locate free variable position
					bvar_val = (bvar_val + free_variables_space[i,fv_index])%2
				basic_solution.append(bvar_val)
			basic_solution_space.append(basic_solution)
		
		basic_solution_space = np.asarray([basic_variables] + basic_solution_space)
		
		full_solution_space = np.zeros([len(free_solutions),matrix.shape[1]])
		for i in range(matrix.shape[1]):
			ctr = 0
			for j in free_variables_space[0]:
				full_solution_space[:,j] = free_variables_space[1:,ctr]
				ctr = ctr + 1
			ctr = 0
			for j in basic_solution_space[0]:
				full_solution_space[:,j] = basic_solution_space[1:,ctr]
				ctr = ctr + 1
		return full_solution_space 

	else:
		print("Invalid Matrix Form")
	
##Obtains a solution space for a linear system representation of the list of parity checks
def linear_system_solver(N,Included_matrix,P_A,partial_key,unconf_pos,Y):

    debugging = True
     # Solve for x using linear algebra techniques
    #################

    Included_matrix_new = Included_matrix.copy()
    P_A_new = P_A.copy()
    dels = 0

    ###build a new matrix that includes only the variables (bits) still unconfident
    for i in range(N):
        if not unconf_pos.__contains__(i):
            conf_x_vec = partial_key[i]*Included_matrix[:,i]
            delpos = i - dels
            Included_matrix_new = np.delete(Included_matrix_new,delpos,1)
            P_A_new = (P_A_new - conf_x_vec) % 2
            dels = dels + 1

    ###remove any all zero rows (occours when all bits are confident in the row)
    temp = Included_matrix_new.copy()
    temp_PA = P_A_new.copy()
    dels = 0
    for i in range(Included_matrix_new.shape[0]):
        if (Included_matrix_new[i,:] == np.zeros(len(unconf_pos))).all():
            temp = np.delete(temp,i-dels,0)
            temp_PA = np.delete(temp_PA,i-dels,0)
            dels = dels + 1
    Included_matrix_new = temp
    P_A_new = temp_PA

    ##convert linear system to reduced row echelon form
    print("Converting to RREF...")
    row_ech,b,row_order = row_echelon_form(Included_matrix_new,P_A_new)

  
    ##Get solution space for linear system
    print("Solving solution of size {0}x{1}...".format(row_ech.shape[0],row_ech.shape[1]))
    num_free_var = row_ech.shape[1] - row_ech.shape[0]
    if num_free_var <= 12:
        calc_soln = find_solution(row_ech,b)
    
        if debugging == True:
            ##Alices actual solution (used only to verify our algoithm works, not available in real scenario)
            correct_soln = []
            unconf_pos.sort()
            for i in unconf_pos:
                correct_soln.append(int(Y[i]))

            found = False
            if calc_soln is not None:
                for soln in calc_soln:
                    solnl = soln.tolist()
                    if (solnl == correct_soln):
                        print("Correct solution contained in solution space of length: {0}".format(calc_soln.shape[0]))
                        found = True
                        break
                if found == False:
                    print("No correct solution contained in solution space of length: {0}".format(calc_soln.shape[0]))

    else:
        calc_soln = []
        found = False
        print("Solution space too large with {0} free variables".format(num_free_var))

    return calc_soln, found


#matrix = np.array([[1,1,0,1,0,0],[1,0,1,0,1,0],[0,1,1,0,0,1]]) 
#b = np.array([[1],[0],[0]])

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
"""