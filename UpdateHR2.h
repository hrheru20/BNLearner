#ifndef UPDATEHR2_H
#define UPDATEHR2_H

#include<stdio.h>
#include<stdlib.h>

#include<iostream>
#include<vector>

#include<deque>
#include<string>
#include<fstream>
#include<math.h>

#include"Arguments.h"
#include"Model.h"

using namespace std;



struct ADTreenode{
	int count;
	
	int startVaryNodeInd;
	
	struct ADVarynode ** ADVaryNodePs;
};

typedef struct ADTreenode ADTreeNode;

struct ADVarynode{

	int index;
	
	int mcv;
	
	struct ADTreenode ** ADTreeNodePs;
	
};

typedef struct ADVarynode ADVaryNode;


typedef struct {

	vector<int> contiTableCounts;
	
	vector<int> contiTableAttrs;
	
	vector<int> condiAttrs;
	vector<int> condiVals;
	
} CondiContiTable;





ADTreeNode * MakeADTree(int, vector<int> &, vector< vector<int> > &, vector<int> &);
ADVaryNode * MakeVaryNode(int, vector<int> &, vector< vector<int> > &, vector<int> &);
void PrintADTreeNode(ADTreeNode *, int, const int, vector<int> &);
void PrintADVaryNode(ADVaryNode *, int, const int, vector<int> &);
void FreeADTreeNode(ADTreeNode *, const int, vector<int> &);
void FreeADVaryNode(ADVaryNode *, const int, vector<int> &);


CondiContiTable MakeContab(deque<int>, ADTreeNode *, 
		vector<int> &, vector<int> &, const vector<int> &);
		
void MinusContab(CondiContiTable &, CondiContiTable *, int, int);

void ConcatContab(CondiContiTable *, int, int, CondiContiTable &);

		



ADTreeNode * MakeADTree(int a_i, vector<int> & dmIndex, vector< vector<int> > & dm, vector<int> & arities){
	//M = num of attr
	int M = dm[0].size();
	
	ADTreeNode * ADTreeNodeP1 = (ADTreeNode *) malloc(sizeof(ADTreeNode));
	ADTreeNodeP1->count = dmIndex.size();
	ADTreeNodeP1->startVaryNodeInd = a_i;
	//base
	if(a_i >= M){
		ADTreeNodeP1->ADVaryNodePs = NULL;
	}
	else{
		ADTreeNodeP1->ADVaryNodePs = (ADVaryNode **) malloc(sizeof(ADVaryNode*) * 
	                              (M - 1 - a_i + 1));
	}
	for(int a_j = a_i; a_j <= M-1; a_j++){
		ADTreeNodeP1->ADVaryNodePs[a_j - a_i] = MakeVaryNode(a_j, dmIndex, dm, arities);
	}
	return ADTreeNodeP1;
	                              
}//end ADTreeNode * MakeADTree()



ADVaryNode * MakeVaryNode(int a_i, vector<int> & dmIndex, vector< vector<int> > & dm, vector<int> & arities){
	ADVaryNode * ADVaryNodeP1 = (ADVaryNode *) malloc(sizeof(ADVaryNode));
	ADVaryNodeP1->index = a_i;

	int n_i = arities[a_i];

	ADVaryNodeP1->ADTreeNodePs = (ADTreeNode **) malloc(sizeof(ADTreeNode *) * n_i);
	vector<int> dmIndexSub[n_i];
	for(unsigned int j = 0; j <= dmIndex.size() - 1; j++){
		int v = dm[dmIndex[j]][a_i];
		dmIndexSub[v].push_back(dmIndex[j]);
	}
	ADVaryNodeP1->mcv = 0;
	int maxCount = dmIndexSub[0].size();
	for(int j = 1; j <= n_i - 1; j++){
		if((int) dmIndexSub[j].size() > maxCount){
			maxCount = dmIndexSub[j].size();
			ADVaryNodeP1->mcv = j;
		}
	}
	for(int j = 0; j <= n_i - 1; j++){
		if(dmIndexSub[j].size() == 0 || j == ADVaryNodeP1->mcv){
			ADVaryNodeP1->ADTreeNodePs[j] = NULL;
		}
		else{
			ADVaryNodeP1->ADTreeNodePs[j] = MakeADTree(a_i+1, dmIndexSub[j], dm, arities);
		}
	}

	return ADVaryNodeP1;
	
	
}//end ADVaryNode * MakeVaryNode()



void PrintADTreeNode(ADTreeNode * ADTreeNodeP1, int level, const int M, vector<int> & arities){
	string strSpace = "";
	for(int i = 0; i < level; i++){
		strSpace += "  ";
	}
	cout << strSpace << "ADTreeNode" << endl;
	cout << strSpace << "Count = " << ADTreeNodeP1->count << endl;
	cout << strSpace << "startVaryNodeInd = " << ADTreeNodeP1->startVaryNodeInd << endl;
	if(ADTreeNodeP1->ADVaryNodePs == NULL){
		return;
	}
	else{
		for(int i = 0; i < M - 1 - ADTreeNodeP1->startVaryNodeInd + 1; i++){
			PrintADVaryNode(ADTreeNodeP1->ADVaryNodePs[i], level+1, M, arities);
		}
	}
} //end void PrintADTreeNode()



void PrintADVaryNode(ADVaryNode * ADVaryNodeP1, int level, const int M, vector<int> & arities){
	string strSpace = "";
	for(int i = 0; i < level; i++){
		strSpace += "  ";
	}
	cout << strSpace << "ADVaryNode" << endl;
	cout << strSpace << "index = " << ADVaryNodeP1->index << endl;
	cout << strSpace << "mcv = " << ADVaryNodeP1->mcv << endl;
	int n_i = arities[ADVaryNodeP1->index];
	for(int j = 0; j <= n_i - 1; j++){
		if(j == ADVaryNodeP1->mcv){
			cout << strSpace << "  mcv" << endl;		
		}
		else if(ADVaryNodeP1->ADTreeNodePs[j] == NULL){
			cout << strSpace << "  NULL" << endl;
		}
		else{
			PrintADTreeNode(ADVaryNodeP1->ADTreeNodePs[j], level+1, M, arities);
		}
	}

}//end void PrintADVaryNode()



void FreeADTreeNode(ADTreeNode * ADTreeNodeP1, const int M, vector<int> & arities){

	if(ADTreeNodeP1->ADVaryNodePs != NULL){
		for(int i = 0; i < M - 1 - ADTreeNodeP1->startVaryNodeInd + 1; i++){
			FreeADVaryNode(ADTreeNodeP1->ADVaryNodePs[i], M, arities);
		}
		free(ADTreeNodeP1->ADVaryNodePs);
	}

	free(ADTreeNodeP1);
}  //end void FreeADTreeNode()



void FreeADVaryNode(ADVaryNode * ADVaryNodeP1, const int M, vector<int> & arities){

	int n_i = arities[ADVaryNodeP1->index]; //may be changed to global var
	for(int j = 0; j <= n_i - 1; j++){
		if(j == ADVaryNodeP1->mcv){
			;		
		}
		else if(ADVaryNodeP1->ADTreeNodePs[j] == NULL){
			;
		}
		else{
			FreeADTreeNode(ADVaryNodeP1->ADTreeNodePs[j], M, arities);
		}
	}
	free(ADVaryNodeP1->ADTreeNodePs);
	free(ADVaryNodeP1);

}//end void FreeADVaryNode()




void print_vec(vector<int> & vec1, ofstream & of1){
	for(int i = 0; i < (int) vec1.size(); i++){
		of1 << vec1[i] << " ";
	}
	of1 << endl;
}

void print_deque(deque<int> & v1, ofstream & of1){
	for(int i = 0; i < (int) v1.size(); i++){
		of1 << v1[i] << " ";
	}
	of1 << endl;
}

void print_CondiContiTable(CondiContiTable & ccTable1, ofstream & of1){
	of1 << "ccTable1.contiTableCounts: " << endl;
	print_vec(ccTable1.contiTableCounts, of1);
	
	of1 << "ccTable1.contiTableAttrs: " << endl;
	print_vec(ccTable1.contiTableAttrs, of1);
	
	of1 << "ccTable1.condiAttrs: " << endl;
	print_vec(ccTable1.condiAttrs, of1);
	
	of1 << "ccTable1.condiVals: " << endl;
	print_vec(ccTable1.condiVals, of1);
	
}



void print_vec(vector<int> & vec1){
	for(int i = 0; i < (int) vec1.size(); i++){
		cout << vec1[i] << " ";
	}
	cout << endl;
}


void print_deque(deque<int> & v1){
	for(int i = 0; i < (int) v1.size(); i++){
		cout << v1[i] << " ";
	}
	cout << endl;
}

void print_CondiContiTable(CondiContiTable & ccTable1){
	cout << "ccTable1.contiTableCounts: " << endl;
	print_vec(ccTable1.contiTableCounts);
	
	cout << "ccTable1.contiTableAttrs: " << endl;
	print_vec(ccTable1.contiTableAttrs);
	
	cout << "ccTable1.condiAttrs: " << endl;
	print_vec(ccTable1.condiAttrs);
	
	cout << "ccTable1.condiVals: " << endl;
	print_vec(ccTable1.condiVals);
	
}



CondiContiTable MakeContab(deque<int> inquiryAttrs, ADTreeNode * ADNP, 
		vector<int> & givenAttrs, vector<int> & givenVals, const vector<int> & arities){
			

	
	//base 1
	if(ADNP == NULL){
		CondiContiTable bottomTable;
		
		//base 1.1
		//if(inquiryAttrs.size() == 0){
		if(inquiryAttrs.empty()){
			bottomTable.contiTableCounts.push_back(0);
			//bottomTable.contiTableAttrs is empty
			bottomTable.condiAttrs = givenAttrs;
			bottomTable.condiVals = givenVals;
			
			//cout << "//base 1.1: if(ADNP == NULL) and if(inquiryAttrs.empty())" << endl;
			//print_CondiContiTable(bottomTable);
			
		}
		//base 1.2
		//if(inquiryAttrs.size() > 0)
		else{
			int rowNum = 1;
			for(unsigned int i = 0; i < inquiryAttrs.size(); i++){

				int lastAttr = inquiryAttrs[inquiryAttrs.size() - 1 - i ];
				rowNum *= arities[lastAttr];
				
				bottomTable.contiTableAttrs.push_back(lastAttr);
			}		
			bottomTable.condiAttrs = givenAttrs;
			bottomTable.condiVals = givenVals;
			
			for(int i = 0; i < rowNum; i++){
				bottomTable.contiTableCounts.push_back(0);
			}
			
			//cout << "//base 1.2: if(ADNP == NULL) and if(!inquiryAttrs.empty())" << endl;
			//print_CondiContiTable(bottomTable);
		}
		
		
		
		return bottomTable;
		
	}// end if(ADNP == NULL)
	
	//base 2
	else if (inquiryAttrs.empty()){
		CondiContiTable bottomTable;
		bottomTable.contiTableCounts.push_back(ADNP->count);
		//bottomTable.contiTableAttrs is empty
		bottomTable.condiAttrs = givenAttrs;
		bottomTable.condiVals = givenVals;
		
		//cout << "//base 2: if(ADNP != NULL) and if(inquiryAttrs.empty())" << endl;
		//print_CondiContiTable(bottomTable);
			
		return bottomTable;
	}
	else{
		int a_i_1 = inquiryAttrs.front();
		int VNIndex = a_i_1 - ADNP->startVaryNodeInd;
		ADVaryNode * VNP = ADNP->ADVaryNodePs[VNIndex];
		int mcv = VNP->mcv;
		int n_i_1 = arities[a_i_1];
		CondiContiTable CTs[n_i_1];
		
		inquiryAttrs.pop_front();
		
		for(int k = 0; k < n_i_1; k++){
			if(k != mcv){
				ADTreeNode * ADNP_k = VNP -> ADTreeNodePs[k];

				vector<int> newGivenAttrs = givenAttrs;
				newGivenAttrs.push_back(a_i_1);
				vector<int> newGivenVals = givenVals;
				newGivenVals.push_back(k);
				CTs[k] = MakeContab(inquiryAttrs, ADNP_k, newGivenAttrs, newGivenVals, arities);
				
			}
		}
		CondiContiTable sumCTs = MakeContab(inquiryAttrs, ADNP, givenAttrs, givenVals, arities);
		MinusContab(sumCTs, CTs, n_i_1, mcv);
		
		CondiContiTable result;
		ConcatContab(CTs, n_i_1, a_i_1, result);
		
		//cout << "//non-base: " << endl;
		//print_CondiContiTable(result);
		
		return result;
	}
	
}// end CondiContiTable MakeContab()


void MinusContab(CondiContiTable & sumCTs, CondiContiTable * CTs, int n_i_1, int mcv){
	CTs[mcv].contiTableCounts = sumCTs.contiTableCounts;
	CTs[mcv].contiTableAttrs = sumCTs.contiTableAttrs;
	for(unsigned int i = 0; i < CTs[mcv].contiTableCounts.size(); i++){
		for(int k = 0; k < n_i_1; k++){
			if(k != mcv){
				CTs[mcv].contiTableCounts[i] -= CTs[k].contiTableCounts[i];
			}
		}
	}
	//assume arities of each attr >=2
	if(mcv != 0){
		CTs[mcv].condiAttrs = CTs[0].condiAttrs;
		CTs[mcv].condiVals = CTs[0].condiVals;
	}
	else{
		CTs[mcv].condiAttrs = CTs[1].condiAttrs;
		CTs[mcv].condiVals = CTs[1].condiVals;
	}
	CTs[mcv].condiVals.pop_back();
	CTs[mcv].condiVals.push_back(mcv);
	
} //end void MinusContab()


void ConcatContab(CondiContiTable * CTs, int n_i_1, int a_i_1, CondiContiTable & concatedTable){
	for(int i = 0; i < n_i_1; i++){

		for(unsigned int j = 0; j < CTs[i].contiTableCounts.size(); j++){
			concatedTable.contiTableCounts.push_back(CTs[i].contiTableCounts[j]);
		}
	}
	
	if(a_i_1 != CTs[0].condiAttrs.back()){
		cout << "Error: a_i_1 != CTs[0].condiAttrs.back()" << endl;
		cout << "a_i_1 = " << a_i_1 << endl;
		cout << "CTs[0].condiAttrs.back() = " << CTs[0].condiAttrs.back() << endl;
	}
	
	concatedTable.contiTableAttrs = CTs[0].contiTableAttrs;
	concatedTable.contiTableAttrs.push_back(a_i_1);
	
	concatedTable.condiAttrs = CTs[0].condiAttrs;
	concatedTable.condiAttrs.pop_back();
	
	concatedTable.condiVals = CTs[0].condiVals;
	concatedTable.condiVals.pop_back();
	
} //end void ContatContab()


		

#endif

