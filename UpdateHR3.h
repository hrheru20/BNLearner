#ifndef UPDATEHR3_H
#define UPDATEHR3_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include<iostream>
#include<string>
#include<fstream>
#include<cmath>
#include<vector>
#include<deque>
#include<list>



#if __GNUC__ < 3 && __GNUC__ >= 2 && __GNUC_MINOR__ >= 95
#  include <hash_set>
#  include <functional>
#  define gnu_namespace std
#elif __GNUC__ >= 3
#  include <ext/hash_set>
#  if __GNUC_MINOR__ == 0
#    include <functional>
#    define gnu_namespace std
#  else
#    include <ext/functional>
#    define gnu_namespace __gnu_cxx
#  endif
#else
#  include <hash_set.h>
#  include <functional.h>
#  define gnu_namespace std
#endif

using namespace gnu_namespace;




using namespace std;

#define MAX_NUM_OF_VARS (32)
#define EPSILON (0.000001)



//HR:
//class for total_order
class total_order{
   	public:
		int torderBP[MAX_NUM_OF_VARS];

		int torder[MAX_NUM_OF_VARS];

		int length;

		double scores;


       	total_order(int torder_length){
       		length = torder_length;
       		scores = 0;

       	}

       	void print(){
       		cout << "total order:" << endl;
       		cout << "length = " << length << endl;
       		cout << "scores = " << scores << endl;

       		cout << "torder : " << endl;
       		for(int i = 0; i < length; i++){
       			cout << "" << torder[i] << "  ";
       		}
       		cout << endl;

       		cout << "torderBP : " << endl;
       		for(int i = 0; i < length; i++){
       			cout << "" << torderBP[i] << "  ";
       		}
       		cout << endl;
       		cout << endl;
       	}


      	friend int operator<(const total_order & to1, const total_order & to2);

      	friend int operator==(const total_order & to1, const total_order & to2);

      	friend int operator>(const total_order & to1, const total_order & to2);

      	friend bool is_same_order(const total_order & to1, const total_order & to2);


};



int operator<(const total_order & to1, const total_order & to2){
 	if(fabs(to1.scores - to2.scores) < EPSILON){
     	return false;
   	}
    return(to1.scores < to2.scores);
}


int operator>(const total_order & to1, const total_order & to2){
 	if(fabs(to1.scores - to2.scores) < EPSILON){
     	return false;
   	}
    return(to1.scores > to2.scores);
}


int operator==(const total_order & to1, const total_order & to2){
	return(fabs(to1.scores - to2.scores) < EPSILON);
}


bool is_same_order(const total_order & to1, const total_order & to2){
	if(!(to1 == to2) || !(to1.length == to2.length)){
		return false;
	}

	bool is_same = true;
	for(int i = 0; i < to1.length; i++){
		if(to1.torder[i] != to2.torder[i]){
			is_same = false;
			break;
		}
	}
    return is_same;

}


bool descend_order_comparison_order(const total_order & to1, const total_order & to2){
	return (to1 > to2);
}




//HR:
//class for dag
class dag{
	public:
	int parents[MAX_NUM_OF_VARS];

	int length;

	double scores;


   	dag(int num_of_var){
   		length = num_of_var;
   		scores = 0;
   		for(int i = 0; i < length; i++){
   			parents[i] = -1;
   		}
   	}

   	dag(){
   		dag(MAX_NUM_OF_VARS);
   	}

   	void print() const{
   		cout << "dag:" << endl;
   		cout << "length = " << length << endl;
   		printf("scores = %18.8f\n", scores);

   		cout << "dag : " << endl;
   		for(int i = 0; i < length; i++){
   			cout << "parents[" << i << "]:  " << parents[i] << endl;
   		}
   		cout << endl;

   	}


  	friend int operator<(const dag & dag1, const dag & dag2);

  	friend int operator==(const dag & dag1, const dag & dag2);

  	friend int operator>(const dag & dag1, const dag & dag2);

  	friend bool is_same_dag(const dag & dag1, const dag & dag2);

  	friend bool haveClosescores(const dag & dag1, const dag & dag2);

};



int operator<(const dag & dag1, const dag & dag2){
 	if(fabs(dag1.scores - dag2.scores) < EPSILON){
     	return false;
   	}
    return(dag1.scores < dag2.scores);
}


int operator>(const dag & dag1, const dag & dag2){
 	if(fabs(dag1.scores - dag2.scores) < EPSILON){
     	return false;
   	}
    return(dag1.scores > dag2.scores);
}


bool descend_order_comparison(const dag & dag1, const dag & dag2){
	return (dag1 > dag2);
}



bool strict_descend_order_comparison(const dag & dag1, const dag & dag2){
	if(dag1 > dag2){
		return true;
	}
	else if(dag1 < dag2){
		return false;
	}
	else{
		for(int i = 0; i < dag1.length; i++){
			if(dag1.parents[i] > dag2.parents[i]){
				return true;
			}
			else if(dag1.parents[i] < dag2.parents[i]){
				return false;
			}
		}

		return false;
	}

}//bool strict_descend_order_comparison()




int operator==(const dag & dag1, const dag & dag2){
	return(fabs(dag1.scores - dag2.scores) < EPSILON);
}


bool is_same_dag(const dag & dag1, const dag & dag2){

	if(!haveClosescores(dag1, dag2) || !(dag1.length == dag2.length)){
		return false;
	}

	bool is_same = true;
	for(int i = 0; i < dag1.length; i++){
		if(dag1.parents[i] != dag2.parents[i]){
			is_same = false;
			break;
		}
	}
    return is_same;

}


bool haveClosescores(const dag & dag1, const dag & dag2){
	return(fabs(dag1.scores - dag2.scores) < 100 * EPSILON);
}


void print_dag_list(list<dag> & dag_list1){
	cout << "\nThe printed dag list has size() = " << dag_list1.size() << endl;
	int dag_index = 0;
	for (list<dag>::iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
		cout << "\ndag_index = " << dag_index << endl;
		(*it1).print();
		dag_index++;
	}

}



void get_ancesters(const list<int> incid_list[], const int Num_vars,
		const int s, bool ** ance_mat1){

	int color[Num_vars];
	for(int u = 0; u < Num_vars; u++){
		color[u] = 0;
	}

	color[s] = 1;

	deque<int> Queue1;
	Queue1.push_back(s);

	while(!Queue1.empty()){
		int u = Queue1.front();
		Queue1.pop_front();


		for (list<int>::const_iterator it1 = incid_list[u].begin();
				it1 != incid_list[u].end(); it1++){
			int v = *it1;
			if(color[v] == 0){
				color[v] = 1;
				Queue1.push_back(v);
			}
		}
		color[u] = 2;

	}

	for(int u = 0; u < Num_vars; u++){
		if(color[u] == 2 && u != s){
			ance_mat1[s][u] = true;
		}
	}


}//void get_ancesters()




void get_ancestor_matrix(const dag & dag1, bool ** ance_mat1){
	const int Num_vars = dag1.length;

	list<int> incid_list[Num_vars];

	for(int v = 0; v < Num_vars; v++){
		int parents = dag1.parents[v];

		for(int u = 0; u < Num_vars; u++){
			if( ((parents >> u) & 1) == 1 ) {
				incid_list[v].push_back(u);
			}
		}
	}


	for(int v = 0; v < Num_vars; v++){
		get_ancesters(incid_list, Num_vars, v, ance_mat1);
	}



}//void get_ancestor_matrix(const dag & dag1, bool ** ance_mat1)




void square_matrix_multiply(const double * const * mat1, const double * const * mat2,
		int n, double * const * result_mat){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){

			result_mat[i][j] = 0;

			for(int k = 0; k < n; k++){
				result_mat[i][j] += mat1[i][k] * mat2[k][j];
			}
		}
	}
}//void square_matrix_multiply()



void square_matrix_add(const double * const * mat1, const double * const * mat2,
		int n, double * const * result_mat){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			result_mat[i][j] = mat1[i][j] + mat2[i][j];
		}
	}
}//void square_matrix_add()



void square_matrix_right_add(double * const * result_mat, const double * const * mat1, int n){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			result_mat[i][j] += mat1[i][j];
		}
	}
}//void square_matrix_right_add()





//HR: get the adjacency matrix for dag1
void get_adj_matrix(const dag & dag1, double * const * adj_mat, int n){

	assert(dag1.length == n);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			adj_mat[i][j] = 0;
		}
	}

	for(int v = 0; v < dag1.length; v++){
		int parents = dag1.parents[v];
		for(int u = 0; u < dag1.length; u++){
			if( ((parents >> u) & 1) == 1 ) {
				adj_mat[u][v] = 1;
			}
		}
	}
}//void get_adj_matrix()



void print_square_matrix(const double * const * mat1, int n){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			cout << "mat1[" << i << "][" << j << "] = " << mat1[i][j] << endl;
		}
	}
}//void print_square_matrix()



void square_matrix_copy(const double * const * mat1, int n, double * const * result_mat){
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			result_mat[i][j] = mat1[i][j];
		}
	}
}//void square_matrix_copy()




void get_reachability_limited_leng(const double * const * adj_mat, int n, double * const * reachability_limited_leng_mat, int limited_leng){

	assert(limited_leng >= 1);

	double *** power_mat_arr1;
	power_mat_arr1 = new double ** [limited_leng];
	for(int i = 0; i < limited_leng; i++){
		power_mat_arr1[i] = new double * [n];
		for(int j = 0; j < n; j++){
			power_mat_arr1[i][j] = new double [n];
		}
	}

	square_matrix_copy(adj_mat, n, power_mat_arr1[0]);
	for(int i = 1; i < limited_leng; i++){
		square_matrix_multiply(power_mat_arr1[i-1], adj_mat, n, power_mat_arr1[i]);
	}

	double ** sum_mat1;
	sum_mat1 = new double *[n];
	for (int i = 0; i < n; i++) {
		sum_mat1[i] = new double[n];
	}

	square_matrix_copy(adj_mat, n, sum_mat1);
	for(int i = 1; i < limited_leng; i++){
		square_matrix_right_add(sum_mat1, power_mat_arr1[i], n);
	}


	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			reachability_limited_leng_mat[i][j] = 0;
		}
	}

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			//reachability_limited_leng_mat[i][j] == 1 iff sum_mat1[i][j] > 0
			if(sum_mat1[i][j] > 0){
				reachability_limited_leng_mat[i][j] = 1;
			}
		}
	}


	for (int i = 0; i < n; i++) {
		delete[] sum_mat1[i];
	}
	delete[] sum_mat1;


	for (int i = 0; i < limited_leng; i++) {
		for(int j = 0; j < n; j++){
			delete[] power_mat_arr1[i][j];
		}
		delete[] power_mat_arr1[i];
	}
	delete[] power_mat_arr1;


}//void get_reachability_limited_leng()




//HR:
//class for family
class family{
	public:
	int child;

	int parents;

   	void print(){
   		cout << "family:" << endl;
   		cout << "child = " << child << endl;
   		cout << "parents = " << parents << endl;
   		cout << endl;

   	}


};




bool is_same_family(const family & family1, const family & family2){
	bool is_same = (family1.child == family2.child) &&
					(family1.parents == family2.parents)
					;

    return is_same;

}





//write a DAG from memory to the file in the hard disk
void write_a_DAG_to_file(const dag & dag1, FILE* fout){

	for(int i = 0; i < MAX_NUM_OF_VARS; i++){
		fwrite(&(dag1.parents[i]), sizeof(int), 1, fout);
	}

	fwrite(&(dag1.length), sizeof(int), 1, fout);
	fwrite(&(dag1.scores), sizeof(double), 1, fout);

}


//write all the DAGs in dag_list1 to the file in the hard disk
void write_DAGs_to_file(list<dag> & dag_list1, FILE* fout){

	int dag_list_size = dag_list1.size();

	fwrite(&dag_list_size, sizeof(int), 1, fout);

	list<dag>::iterator it1;
	for (it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
		write_a_DAG_to_file((*it1), fout);
	}
}



//read a DAG from the file in the hard disk to memory
void read_a_DAG_from_file(dag & dag1, FILE* fin){
	for(int i = 0; i < MAX_NUM_OF_VARS; i++){
		fread(&(dag1.parents[i]), sizeof(int), 1, fin);
	}

	fread(&(dag1.length), sizeof(int), 1, fin);
	fread(&(dag1.scores), sizeof(double), 1, fin);

}



//read all the stored DAGs from the file in the hard disk and stored them to dag_list1
void read_DAGs_from_file(list<dag> & dag_list1, FILE* fin){
	int dag_list_size;

	fread(&dag_list_size, sizeof(int), 1, fin);

	for (int i = 0; i < dag_list_size; i++){
		dag dag1(MAX_NUM_OF_VARS);
		read_a_DAG_from_file(dag1, fin);
		dag_list1.push_back(dag1);

	}

}





struct hash_dag{
    size_t operator()(const dag & dag1) const{
    	hash<int> hash_value;
    	return hash_value((int) round(- dag1.scores * 10000));
    }
};

struct eq_dag{
    bool operator()(const dag & dag1, const dag & dag2) const {
    	return is_same_dag(dag1, dag2);
    }
};





typedef hash_set<dag, hash_dag, eq_dag> DAGHashSet;




#endif



