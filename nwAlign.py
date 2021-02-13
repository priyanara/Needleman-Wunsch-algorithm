#!/usr/bin/env python3
import sys

#assigning values
match_value=1
mismatch_value=-1
gap_penalty=-1

#access the file to get the length of the sequence
with open(sys.argv[1],'r') as fh1:
	line1=fh1.readlines()
	first_sequence=line1[1].strip()

with open(sys.argv[2],'r') as fh2:
	line2=fh2.readlines()
	second_sequence=line2[1].strip() 
	#remove the leading and trailing spaces
	second_sequence
#calling the sequence
row_matrix=len(first_sequence)
column_matrix=len(second_sequence)

matrix=[]
#empty matrix to which zero value is added for it to be initialized
#checking to see if the character in the string match or not
def assignment_score(a,b):
	if a == b:
		return match_value
	elif a == '-' or b == '-':
		return gap_penalty
	else:
		return mismatch_value

#create a zero matrix which with the same dimensions as the required matrix
def initial_matrix(row_matrix,column_matrix):
	for i in range(row_matrix):
		sub_matrix=[]
		for j in range(column_matrix):
			sub_matrix.append(0)
		matrix.append(sub_matrix)
	for i in range(row_matrix):
		for j in range(column_matrix):
			matrix[i][j]
	return matrix

#to determine the score between two bases in the alignment 

#to assign scores to each cell of the matrix according to the two sequences
def main(first_sequence,second_sequence):

#initialisation
#creating a zero matrix 
	first_matrix=initial_matrix(row_matrix+1,column_matrix+1)

	for i in range(0,row_matrix+1):
		first_matrix[i][0]=gap_penalty*i

	for j in range(0,column_matrix+1):
		first_matrix[0][j]=gap_penalty*j

#filling
	for i in range(1,row_matrix+1):
		for j in range(1,column_matrix+1):
			x = first_matrix[i - 1][j - 1] + assignment_score(first_sequence[i-1], second_sequence[j-1])
			y = first_matrix[i - 1][j] + gap_penalty
			z = first_matrix[i][j - 1] + gap_penalty
			first_matrix[i][j]=max(x,y,z)

	return first_matrix

#empty list to add the matched string
seq1=[]
seq2=[]

i=row_matrix
j=column_matrix

symbol_match='|'
symbol_mismatch='*'
#correct till here, do not check 

#string indexing and list indexing to be kept in mind
#back-tracing step 

def backtracing(i,j,seq1,seq2,first_matrix):
	while i>0 and j>0:
		matrix_current=first_matrix[i][j]
		matrix_diagonal=first_matrix[i-1][j-1]
		matrix1=first_matrix[i][j-1]
		matrix2=first_matrix[i-1][j]

		if matrix_current==matrix_diagonal + assignment_score(first_sequence[i-1],second_sequence[j-1]):
			seq1.append(first_sequence[i-1])
			seq2.append(second_sequence[j-1])
			i-=1
			j-=1
		elif matrix_current==matrix2 + gap_penalty:
			seq1.append(first_sequence[i-1])
			seq2.append('-')
			i-=1
		elif matrix_current==matrix1 + gap_penalty:
			seq1.append('-')
			seq2.append(second_sequence[j-1])
			j-=1

#since the sequence is appended from the last cell, has to be reversed 
	return(seq1[::-1],seq2[::-1]) 

a=main(first_sequence,second_sequence)
b=backtracing(i,j,seq1,seq2,a)

#getting the alignment and the pattern  
#list to get the pattern 
final=[]
alignment_score=0
#if the nucleotides match then the symbols are set
#alignment scores are obtained by calculating the match and the mismatch

for x,y in zip(b[0],b[1]):
	if x==y:
		final.append(symbol_match)
		alignment_score+=1
	elif x=='-' or y=='-':
		final.append(' ')
		alignment_score-=1
	else:
		final.append(symbol_mismatch)
		alignment_score-=1

print(''.join(b[0]))
print(''.join(final))
print(''.join(b[1]))
print("Alignment score: "+ str(alignment_score))



		

		
