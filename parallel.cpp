			/* Include Header Files */
#include<iostream>
#include<omp.h>
#include<time.h>
#include<cmath>
#include<fstream>
#include<string.h>
#include<string>
#include<map>
#include<vector>

using namespace std;

#define THREADS 10

int mat_mul(vector<vector<double> > M_mat, vector<double> &rank,int n)
{
	long  a, b;
	vector<double> new_rank(n);
	bool flag = false;

	#pragma omp parallel num_threads(THREADS) shared(new_rank)
	{
		#pragma omp for collapse(2) 
		for(int i=0;i<n;i++)
		{
			for(int j =0;j<n;j++)
			{
				#pragma omp atomic
				new_rank[i] += M_mat[i][j]*rank[j];
			}
		}
	}
		
	for(int i=0;i<n;i++)
	{	
		flag = false;		
		a = rank[i] * 1000000;
		b = new_rank[i] * 1000000;
		if(a == b)
			continue;
		else
		{
			flag = true;
			break;
		}
	}	

	if( ! flag)
		return 1;

	#pragma omp parallel for simd shared(rank)
	for(int i=0;i<n;i++)
	{	
		rank[i] = new_rank[i];
	}

	return 0;
}

void calc_rank(long n)
{

	vector<long>outdegree(n,0); // to store out degree of each node
	vector<vector<double> > A_mat(n, vector<double>(n,0)); // to store 1/L
	
	double rank_sum = 0.0; // Final rank sum of all page must be 1.0
	double damping_factor = 0.85;   // randomly given damping factor
	double init_rank = 1/(double)n; // initialise initial rank to  1/number of nodes
	
	vector<double>rank(n,init_rank);  // to store the final rank
	long in_node, out_node, node; // used to read in and out node from file
	long i, j, k; // iterator variables
	vector<vector<double> > M_mat(n, vector<double>(n,init_rank)); // to store damping factor equation
	
	
								/* read the file name*/

	string str = to_string(n);
	char file_name[30];
	strcpy(file_name,str.c_str());
	strcat(file_name,"_nodes.txt");
	
	ifstream file_read;
	file_read.open(file_name, ios::in);

	int no_lines = 0;

							/* Calculate outdegree of each node */

	while(file_read>>out_node>>in_node)
	{
		no_lines++;
		outdegree[out_node]++;
		
	}
	file_read.clear();
	file_read.seekg(0);

						/* Initialise each rank on the basis of 1/outdegree */
	long in = 0;
	while(file_read>>out_node>>in_node)
	{
		if(outdegree[out_node])
			A_mat[in_node][out_node] += 1/double(outdegree[out_node]);
		in++;

	}
	file_read.close();

							/* loops are independent so use collapse */

	#pragma omp parallel for collapse(2) num_threads(THREADS)
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			M_mat[i][j] = (damping_factor * A_mat[i][j]) + ((1 - damping_factor) * init_rank);
		}
	}

	k = 0;
	while(1)
	{	
		k++;	
		cout<<"Running Iteration "<<k<<"...\n";
		int ret_val = mat_mul(M_mat,rank,n);
		if(ret_val == 1)
			break;
	}

					/* calculate total rank sum */

	rank_sum = 0;
	cout << "\n********************************************************\n";
	cout << "Page No.\tRank\n";
	#pragma omp parallel for reduction(+:rank_sum)
	for(i = 0 ; i < n ; i++)
	{
		rank_sum += rank[i];
		cout<<"   "<< i << "\t\t" <<rank[i] <<"\n";
	}

						/*  PRINT STATEMENTS */


	cout << "********************************************************\n";
	cout << "NUMBER OF PAGES: " << n << endl;
	cout << "NUMBER OF LINKS: " << no_lines << endl;
	cout<<"TOTAL RANKSUM OF EACH PAGE: "<< rank_sum<<"\n";
	cout<<"TOTAL NUMBER OF ITERATIONS REQUIRED TO CONVERGE: "<< k;
}

int main()
{	
	cout << "\033[2J\033[1;1H";
	long n;
	cout << "Enter the number of pages(no. of graph nodes): ";
	cin >> n;
	cout << "\033[2J\033[1;1H";
	cout << "*********************** RESULTS ************************\n\n";
	double  start_time;
	start_time = omp_get_wtime();
	calc_rank(n); 
	double total_time = omp_get_wtime() - start_time;
	cout<<"\nTOTAL TIME TAKEN WITH PARALLEL ALGORITHM: " << total_time << endl;
	cout<<"NUMBER OF THREADS "<<THREADS<< endl << endl ;
	return 0;
}
