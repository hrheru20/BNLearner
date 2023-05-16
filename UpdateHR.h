#ifndef UPDATEHR_H
#define UPDATEHR_H

#include<stdio.h>
#include<stdlib.h>
#include"math.h"

#include"Arguments.h"
#include"Model.h"
#include"Engine.h"



#define LOGADD_NEW(logA, logB) if(logA == MARK) logA = logB; else if (logA == LOG_ZERO) logA = logB; else if (logB == LOG_ZERO) logA = logA; else logA = logA + log(1+exp(logB-logA));
    	

#define LOGMINUS_NEW(logA, logB) if(logA == MARK) logA = logB; else if (logA == LOG_ZERO) logA = logB; else if (logB == LOG_ZERO) logA = logA; else logA = logA + log(1-exp(logB-logA));


#define epsilon 1.0e-5

bool logSpecValEquals(double logX, double logSpecVal){
	return (logX >= logSpecVal - epsilon && logX <= logSpecVal + epsilon);
}


void logAdd_New(double & logA, double & logB){
	    //if(logA == MARK){
	    if(logSpecValEquals(logA, MARK)){
    		logA = logB;
    	}
    	//else if (logA == LOG_ZERO){
    	else if (logSpecValEquals(logA, LOG_ZERO)){
    		logA = logB;
    	}
    	//else if(logB == MARK){
    	else if(logSpecValEquals(logB, MARK)){
    		logA = logA;
    	}
    	//else if(logB == LOG_ZERO){
    	else if(logSpecValEquals(logB, LOG_ZERO)){
    		logA = logA;
    	}
    	//Neither logA nor logB is MARK/LOG_ZERO
    	//logA < 0 && logB < 0
    	else{
    		//1.1
			//Gaurd logB-logA <= 100 to avoid that exp(logB-logA) go to inf
    		if(logB - logA <= 100 ){
    			logA = logA + log( 1 + exp(logB-logA) );

    		} 
    		//1.2		
    		//if logB - logA > 100, then logB > 100 + logA so that B is much bigger than A (B > e^100 * A)
    		//    So (A + B) almost = B 
    		else{

	    		logA = logB + log( exp(logA - logB) + 1);
    		}
    		
    	}
}// end of void logAdd_New(double & logA, double & logB)


void logMinus_New(double & logA, double & logB){
	    //if(logA == MARK){
	    if(logSpecValEquals(logA, MARK)){
    		cerr << "Error: logMinus_New(double & logA, double & logB): logA == MARK\n"; 
    		exit(1);
    	}
    	//else if (logA == LOG_ZERO && logB != LOG_ZERO){
    	else if (logSpecValEquals(logA, LOG_ZERO) && ! logSpecValEquals(logB, LOG_ZERO)){
    		cerr << "Error: logMinus_New(double & logA, double & logB): logA == LOG_ZERO && logB != LOG_ZERO\n"; 
    		exit(1);
    	}
    	//else if(logB == MARK){
    	else if(logSpecValEquals(logB, MARK)){
    		logA = logA;
    	}
    	//else if(logB == LOG_ZERO){
    	else if(logSpecValEquals(logB, LOG_ZERO)){
    		logA = logA;
    	}
    	//Nither logA nor logB is MARK/LOG_ZERO
    	//logA >= logB
    	else{
    		if(logA == logB){
    			logA = LOG_ZERO;
    		}
    		//logA > logB
    		else{
    			if(logA > logB){
    				logA = logA + log(1 - exp(logB-logA));
    			}
    			//if logA < logB

    			else{
    				cerr << "Warning: logMinus_New(double & logA, double & logB): logA < logB" << endl;
    				logA = -( logB + log(1 - exp(logA-logB)) );
    			}

    			
    		}
    		
    	}
} // end of void logMinus_New(double & logA, double & logB)



void logMinus_New_NonnegRes(double & logA, double & logB){
	    //if(logA == MARK){
	    if(logSpecValEquals(logA, MARK)){
    		cerr << "Error: logMinus_New_NonnegRes(double & logA, double & logB): logA == MARK\n"; 
    		exit(1);
    	}
    	//else if(logB == MARK){
    	else if(logSpecValEquals(logB, MARK)){
    		cerr << "Error: logMinus_New_NonnegRes(double & logA, double & logB): logB == MARK\n"; 
    		exit(1);
    	}
    	//else if(logB == LOG_ZERO){
    	else if(logSpecValEquals(logB, LOG_ZERO)){
    		logA = logA;
    	}
    	//Nither logA nor logB is MARK/LOG_ZERO
    	//logA >= logB
    	else{
    		if(logA == logB){
    			logA = LOG_ZERO;
    		}
    		//logA > logB
    		else{
    			if(logA > logB){
    				logA = logA + log(1 - exp(logB-logA));
    			}
    			//if logA < logB, return LOG_ZERO as the bound
    			else{
    				logA = LOG_ZERO;
    			}
	
    		}
    		
    	}
} // end of void logMinus_New_NonnegRes(double & logA, double & logB)



void logAddComp(double & logA, double & logB){

    	//if(logA == MARK){
	    if(logSpecValEquals(logA, MARK)){
    		logA = logB;
    	}
    	//else if (logA == LOG_ZERO){
    	else if (logSpecValEquals(logA, LOG_ZERO)){
    		logA = logB;
    	}
    	//else if(logB == MARK){
    	else if(logSpecValEquals(logB, MARK)){
    		logA = logA;
    	}
    	//else if(logB == LOG_ZERO){
    	else if(logSpecValEquals(logB, LOG_ZERO)){
    		logA = logA;
    	}
    	
    	
    	//Nither logA nor logB is MARK/LOG_ZERO
    	else{
    		//1
    		if(logA < 0 && logB < 0 ){
    			//1.1
				//Gaurd logB-logA <= 100 to avoid that exp(logB-logA) go to inf
	    		if(logB - logA <= 100 ){
	    			logA = logA + log( 1 + exp(logB-logA) );
	
	    		} 
	    		//1.2		
	    		//if logB - logA > 100, then logB > 100 + logA so that B is much bigger than A (B > e^100 * A)
	    		//    So (A + B) almost = B 
	    		else{

		    		logA = logB + log( exp(logA - logB) + 1);
	    		}
    		}
    		
    		//2
    		else if(logA > 0 && logB > 0){
    			//always make exp(.): . > 0
    			if(-logB + logA <= 100){
    				//2.1
    				logA = -( -logA + log( 1 + exp(-logB+logA) ) );
    			} 			
    			else{
    				//2.2
    				logA = -( -logB + log( 1 + exp(-logA+logB) ) );
    			}
    			
    		}
    		//3
    		else if(logA < 0 && logB > 0){
    			//3.1
    			if(logA + logB > 0){   	
    				if(logA + logB <= 100){
    					//3.1.2
    					logA = -logB + log( exp(logA + logB) - 1 );
    				}	
    				else{
    					//3.1.1
    					logA = logA + log( 1 - exp(-logB-logA) );
    				}	
	
    			}
    			//3.2
    			else if(logA + logB < 0){
    				if(-logA - logB <= 100){
    					//3.2.2
    					logA = -( logA + log( exp(-logB-logA)  - 1) );
    				}
    				else{
    					//3.2.1
    					logA = - ( -logB + log(1 - exp(logA+logB) ) );
    				}

    			}
    			//3.3
    			else{
    				logA = LOG_ZERO;
    			}
    		}
    		//4
    		else if(logA > 0 && logB < 0){
    			//4.1
    			if(logA + logB > 0){
    				if(logA + logB <= 100){
    					//4.1.2
    					logA = -logA + log( exp(logB+logA) - 1);
    				}
    				else{
    					//4.1.1
    					logA = logB + log( 1 - exp(-logA-logB) );
    				}
    				
    			}
    			//4.2
    			else if(logA + logB < 0){
    				if(-logA-logB <= 100){
    					//4.2.2
    					logA = -( logB + log( exp(-logA-logB) - 1) );
    				}
    				else{
    					//4.2.1
    					logA = -( -logA + log(1 - exp(logB+logA)) );
    				}
	
    			}
    			//4.3
    			else{
    				logA = LOG_ZERO;
    			}			
    		}
		
    	}
	
} // end of void logAddComp(double & logA, double & logB)

    	
    	

//forward
double * hFunc;
double * hFuncNume;


//backward
double * AA_S;

double * RRFunc;
double * RRFuncNume;


//mixIndegree
double * Eta_U;
double ** KFunc;

//In_Out
double * RRFuncInOut;



//ofstream totalOfMixEta("debugEngineMixForEta_U.txt");




void printIntArr(int * TSet, int leng){
		for(int ind = 0; ind < leng; ind++){
			cout << "Set[" << ind << "]=" << TSet[ind] << endl;
		}	
}



int binIntLen(int S){
	int SSetLen = 0;

	while(S > 0){ 
		if(S & 1){
			SSetLen++;
			
		}
		S >>= 1;
	}
	return SSetLen;
	

}



void binIntToDecIntArr(int S, int ones, int *SSet){
	int posFromLeft = 0;
	int ind = 0;
	while(S > 0){ 
		if(S & 1){
			SSet[ind++] = posFromLeft;
			
		}
		S >>= 1;
		posFromLeft++;
	}
	//printIntArr(SSet, ones);
	

	
}



//HR: Compute F(S, T)
double compute_f(int nh, int S, int T, int sizeofT, int *TSet, double * h, double ** alpha){
	

	double result_f = MARK;
	double result_h = h[S-T];
	//LOGADD(result_f, result_h);
	result_f = result_h;
	if (result_f == LOG_ZERO){
		return result_f;
	}

	
	for(int index = 0; index < sizeofT; index++){
		int j = TSet[index];
		double result_A = alpha[j][S-T];
		//result_f += result_A;
		if (result_A == LOG_ZERO){
			return LOG_ZERO;
		}
		else{
			result_f += result_A;
		}
	}
	return result_f;
}




void rec_compute_sum_F_S_T(
	int d,
	int nh, 
	int S, 
	int ones, 
	int * SSet, 
	int T,
	int sizeofT,
	int * TSet,
	double * h, 
	double ** alpha, 
	double & sumOdd, 
	double & sumEven)
{

	if(d < nh){	
		rec_compute_sum_F_S_T(d+1, nh, S, ones, SSet, T, sizeofT, TSet, h, alpha, sumOdd, sumEven);
		
		if(((S>>d) & 1) == 1){
			TSet[sizeofT] = d;
			rec_compute_sum_F_S_T(d+1, nh, S, ones, SSet, (T | (1<<d)), (sizeofT + 1), TSet, h, alpha, sumOdd, sumEven);
		}
	}
	else{
		if(sizeofT == 0){
			return;
			
		}
		else{	
			double result_F = compute_f(nh, S, T, sizeofT, TSet, h, alpha);
			//sizeofT is odd
			if(sizeofT % 2 == 1){				
				logAdd_New(sumOdd, result_F);
				//LOGADD_NEW(sumOdd, result_F);
			}
			else{
				logAdd_New(sumEven, result_F);
				//LOGADD_NEW(sumEven, result_F);
			}
			
		}
		
	}
} // end rec_compute_sum_F_S_T()

	



void compute_sum_F_S_T(int nh, int S, int ones, int * SSet, double * h, double ** alpha, double & sumOdd, double & sumEven){
	int * TSet = new int[ones];
	rec_compute_sum_F_S_T(0, nh, S, ones, SSet, 0, 0, TSet, h, alpha, sumOdd, sumEven);
	delete[] TSet;
	
} // end compute_sum_F_S_T()




void sub_h_Fast(int d, int nh, int S, int ones, int * SSet, double * h, double ** alpha){
	if(d < nh){
		
		sub_h_Fast(d+1, nh, S, ones, SSet, h, alpha);
		
		SSet[ones] = d;
		sub_h_Fast(d+1, nh, S|(1<<d), (ones+1), SSet, h, alpha);
	}

	else{
		if(S == 0){
			h[0] = 0; //log(1) = 0;
			//cout << "h[0] = 0" << endl;
			//cout << endl;
			return;

		}
		
		double sumOdd = MARK;
		double sumEven = MARK;
		
		
		compute_sum_F_S_T(nh, S, ones, SSet, h, alpha, sumOdd, sumEven);
		
		if(sumEven != MARK){
			logMinus_New(sumOdd, sumEven);
			//LOGMINUS_NEW(sumOdd, sumEven);
		}
						
		h[S] = sumOdd;
		//cout << "h[" << S << "] = "  << h[S] << endl;
		//cout << endl;

	}
} // end sub_h_Fast()


//HR:
//	Compute h
void compute_h_Fast(int nh, double * h, double ** alpha){
	int * SSet = new int[nh];
	
	sub_h_Fast(0, nh, 0, 0, SSet, h, alpha);
	
	delete [] SSet;
	
} // end compute_h_Fast()









void compute_sum_F_S_T_Iter(int nh, int S, double * h, double ** alpha, double & sumOdd, double & sumEven){
	
	for(int T = 1; T <= S; T++){

		if((S|T) == S){
			
			int TSetLen = binIntLen(T);
			int * TSet = new int[TSetLen];
			binIntToDecIntArr(T, TSetLen, TSet);
			
			double result_F = compute_f(nh, S, T, TSetLen, TSet, h, alpha);
			//sizeofT is odd
			if(TSetLen % 2 == 1){				
				logAdd_New(sumOdd, result_F);
				//LOGADD_NEW(sumOdd, result_F);
			}
			else{
				logAdd_New(sumEven, result_F);
				//LOGADD_NEW(sumEven, result_F);
			}
			
			delete[] TSet;
	
		}
	
	}
	
}


//HR:
//	Compute h by iteration instead of recursion
void compute_h_Iter(int nh, double * h, double ** alpha){
	
	//
	h[0] = 0;

	//h[S]: S: from 1 to (1<<nh)-1
	for(int S = 1; S < (1<<nh); S++){
		
		double sumOdd = MARK;
		double sumEven = MARK;
		
		compute_sum_F_S_T_Iter(nh, S, h, alpha, sumOdd, sumEven);
		
		//if sumEven has been assigned some value instead of MARK
		if(sumEven != MARK){
			logMinus_New(sumOdd, sumEven);
			//LOGMINUS_NEW(sumOdd, sumEven);
		}
						
		h[S] = sumOdd;
	}
	
}




void rec_compute_AA_S(
	int d,
	int nh,
	int S,
	int ones, //leng of SSet
    int * SSet,
    int T,
    int sizeofT,
    int * TSet,
    double * AA_S,
    double ** alpha
){


	if(d < nh){	
		rec_compute_AA_S(d+1, nh, S, ones, SSet, T, sizeofT, TSet, AA_S, alpha);
		
		if(((S>>d) & 1) == 1){
			TSet[sizeofT] = d;
			rec_compute_AA_S(d+1, nh, S, ones, SSet, (T | (1<<d)), (sizeofT + 1), TSet, AA_S, alpha);
		}
	}
	else{
		if(T == 0){
			//cout << "\nT == 0" << endl;
			
			AA_S[T] = 0;
			
			//cout << "AA_S[" << T << "]=" << AA_S[T] << endl;
			//cout << endl;
		}
		else{

			if(alpha[TSet[0]][(1<<nh) - 1 - S] == LOG_ZERO){
				AA_S[T] = LOG_ZERO;
			}
			else{
				AA_S[T] = alpha[TSet[0]][(1<<nh) - 1 - S] + AA_S[T - (1<<TSet[0])];
			}

			//cout << "AA_S[" << T << "]=" << AA_S[T] << endl;
			//cout << endl;
		}
		
	}


}


//	Given S \subset V,
//  Compute all T \subset S: AA_S(T) at one time
void compute_AA_S(
	int nh,
	int S,
	int ones, //leng of SSet
    int * SSet,
    double * AA_S,
    double ** alpha
){
	int * TSet = new int[ones];
	rec_compute_AA_S(0, nh, S, ones, SSet, 0, 0, TSet, AA_S, alpha);
	delete[] TSet;
	
}


void compute_AA_S_Iter(
	int nh,
	int S,
    double * AA_S,
    double ** alpha
){

    AA_S[0] = 0;
    
	for(int T = 1; T <= S; T++){
		//if T is the subset of S
		if((S|T) == S){
			
			int TSetLen = binIntLen(T);
			int * TSet = new int[TSetLen];
			binIntToDecIntArr(T, TSetLen, TSet);
			
			AA_S[T] = alpha[TSet[TSetLen-1]][(1<<nh) - 1 - S] + AA_S[T - (1<<TSet[TSetLen-1])];
			
			delete[] TSet;
	
		}
	
	}
	
}






//HR: Compute RF(S, T)
double compute_RF(int nh, int S, int T, int sizeofT, int *TSet, double * RR, double ** alpha){

	
	double result_RF = MARK;
	double result_RR = RR[S-T];
	//LOGADD(result_f, result_h);
	result_RF = result_RR;
	if (result_RF == LOG_ZERO){
		return result_RF;
	}

	if(AA_S[T] == LOG_ZERO){
		return LOG_ZERO;
	}
	else{
		result_RF += AA_S[T];
	}
	
	return result_RF;
}




//HR: Compute: \sum_{k=1}^{|S|} (-1)^{k+1} \sum{T \subset S, |T|=k} RF(S,T)
void rec_compute_sum_RF_S_T(
	int d,
	int nh, 
	int S, 
	int ones, 
	int * SSet, 
	int T,
	int sizeofT,
	int * TSet,
	double * RR, 
	double ** alpha, 
	double & sumOdd, 
	double & sumEven)
{

	if(d < nh){	
		rec_compute_sum_RF_S_T(d+1, nh, S, ones, SSet, T, sizeofT, TSet, RR, alpha, sumOdd, sumEven);
		
		if(((S>>d) & 1) == 1){
			TSet[sizeofT] = d;
			rec_compute_sum_RF_S_T(d+1, nh, S, ones, SSet, (T | (1<<d)), (sizeofT + 1), TSet, RR, alpha, sumOdd, sumEven);
		}
	}
	else{
		if(sizeofT == 0){
			return;
			
		}
		else{	
			double result_RF = compute_RF(nh, S, T, sizeofT, TSet, RR, alpha);
			//sizeofT is odd
			if(sizeofT % 2 == 1){				
				logAdd_New(sumOdd, result_RF);
				//LOGADD_NEW(sumOdd, result_RF);
			}
			else{
				logAdd_New(sumEven, result_RF);
				//LOGADD_NEW(sumEven, result_RF);
			}
			
		}
		
	}
} // end rec_compute_sum_RF_S_T()

	


//HR: Compute: \sum_{k=1}^{|S|} (-1)^{k+1} \sum{T \subset S, |T|=k} RF(S,T)

void compute_sum_RF_S_T(int nh, int S, int ones, int * SSet, double * RR, double ** alpha, double & sumOdd, double & sumEven){
	
	int * TSet = new int[ones];
	rec_compute_sum_RF_S_T(0, nh, S, ones, SSet, 0, 0, TSet, RR, alpha, sumOdd, sumEven);
	delete[] TSet;
	
} // end compute_sum_RF_S_T()








void sub_RR_Fast(int d, int nh, int S, int ones, int * SSet, double * RR, double ** alpha){
	if(d < nh){
		sub_RR_Fast(d+1, nh, S, ones, SSet, RR, alpha);
		
		SSet[ones] = d;
		sub_RR_Fast(d+1, nh, S|(1<<d), (ones+1), SSet, RR, alpha);
	}

	else{

		if(S == 0){
			RR[0] = 0; //log(1) = 0;
			//cout << "R[0] = 0" << endl;
			//cout << endl;
			return;

		}
		
		AA_S = new double [S+1];
		compute_AA_S(nh, S, ones, SSet, AA_S, alpha);
		
		double sumOdd = MARK;
		double sumEven = MARK;
		
		
		compute_sum_RF_S_T(nh, S, ones, SSet, RR, alpha, sumOdd, sumEven);
		
		
		if(sumEven != MARK){
			logMinus_New(sumOdd, sumEven);
			//LOGMINUS_NEW(sumOdd, sumEven);
		}
	
		RR[S] = sumOdd;


		delete[] AA_S;

	}
} // end sub_RR_Fast()







//HR:
//	Compute RR
void compute_RR_Fast(int nh, double * RR, double ** alpha){
	int * SSet = new int[nh];
	
	sub_RR_Fast(0, nh, 0, 0, SSet, RR, alpha);
	
	delete[] SSet;
	
} // end compute_RR_Fast()


//HR: Compute: \sum_{k=1}^{|S|} (-1)^{k+1} \sum{T \subset S, |T|=k} RF(S,T)
void rec_compute_sum_RF_S_T_TSize(
	int d,
	int nh, 
	int S, 
	int ones, 
	int * SSet, 
	int T,
	int sizeofT,
	int * TSet,
	double * RR, 
	double ** alpha, 
	double & sumOdd, 
	double & sumEven,
	int TSize)
{

	if(d < nh){	
		rec_compute_sum_RF_S_T_TSize(d+1, nh, S, ones, SSet, T, sizeofT, TSet, RR, alpha, sumOdd, sumEven, TSize);
		
		if(((S>>d) & 1) == 1){
			TSet[sizeofT] = d;
			rec_compute_sum_RF_S_T_TSize(d+1, nh, S, ones, SSet, (T | (1<<d)), (sizeofT + 1), TSet, RR, alpha, sumOdd, sumEven, TSize);
		}
	}
	else{
		if(sizeofT == 0){
			return;
			
		}
		else if(sizeofT <= TSize){
			double result_RF = compute_RF(nh, S, T, sizeofT, TSet, RR, alpha);
			//sizeofT is odd
			if(sizeofT % 2 == 1){				
				logAdd_New(sumOdd, result_RF);
				//LOGADD_NEW(sumOdd, result_RF);
			}
			else{
				logAdd_New(sumEven, result_RF);
				//LOGADD_NEW(sumEven, result_RF);
			}
			
		}
		
	}
} // end rec_compute_sum_RF_S_T_TSize()



void compute_sum_RF_S_T_TSize(int nh, int S, int ones, int * SSet, double * RR, double ** alpha, double & sumOdd, double & sumEven, int TSize){
	
	int * TSet = new int[ones];
	rec_compute_sum_RF_S_T_TSize(0, nh, S, ones, SSet, 0, 0, TSet, RR, alpha, sumOdd, sumEven, TSize);
	delete[] TSet;
	
	
}// end compute_sum_RF_S_T_TSize()



void sub_RR_Fast_TSize(int d, int nh, int S, int ones, int * SSet, double * RR, double ** alpha, int TSize){
	if(d < nh){
		sub_RR_Fast_TSize(d+1, nh, S, ones, SSet, RR, alpha, TSize);
		
		SSet[ones] = d;
		sub_RR_Fast_TSize(d+1, nh, S|(1<<d), (ones+1), SSet, RR, alpha, TSize);
	}
	else{
		if(S == 0){
			RR[0] = 0; //log(1) = 0;
			//cout << "R[0] = 0" << endl;
			//cout << endl;
			return;

		}
		
		AA_S = new double [S+1];
		compute_AA_S(nh, S, ones, SSet, AA_S, alpha);	
		
		double sumOdd = MARK;
		double sumEven = MARK;
		
		
		compute_sum_RF_S_T_TSize(nh, S, ones, SSet, RR, alpha, sumOdd, sumEven, TSize);

		
		//if sumEven has been assigned some value instead of MARK
		if(sumEven != MARK){
			logMinus_New_NonnegRes(sumOdd, sumEven);
			//logMinus_New(sumOdd, sumEven);
			//LOGMINUS_NEW(sumOdd, sumEven);
		}
	
		RR[S] = sumOdd;

		delete[] AA_S;

	}
}// end sub_RR_Fast_TSize()



void compute_RR_Fast_TSize(int nh, double * RR, double ** alpha, int TSize){
	int * SSet = new int[nh];
	
	sub_RR_Fast_TSize(0, nh, 0, 0, SSet, RR, alpha, TSize);

	delete[] SSet;
	
} // end compute_RR_Fast_TSize()





void rec_compute_sum_F_S_T_TSize(
	int d,
	int nh, 
	int S, 
	int ones, 
	int * SSet, 
	int T,
	int sizeofT,
	int * TSet,
	double * h, 
	double ** alpha, 
	double & sumOdd, 
	double & sumEven,
	int TSize)
{

	if(d < nh){	
		rec_compute_sum_F_S_T_TSize(d+1, nh, S, ones, SSet, T, sizeofT, TSet, h, alpha, sumOdd, sumEven, TSize);
		
		if(((S>>d) & 1) == 1){
			TSet[sizeofT] = d;
			rec_compute_sum_F_S_T_TSize(d+1, nh, S, ones, SSet, (T | (1<<d)), (sizeofT + 1), TSet, h, alpha, sumOdd, sumEven, TSize);
		}
	}
	else{
		if(sizeofT == 0){

			return;
			
		}

		else if(sizeofT <= TSize){
			double result_F = compute_f(nh, S, T, sizeofT, TSet, h, alpha);
			//sizeofT is odd
			if(sizeofT % 2 == 1){		
			//if(sizeofT % 2 == 1 && sizeofT <= TSize){			
				logAdd_New(sumOdd, result_F);
				//LOGADD_NEW(sumOdd, result_F);
			}
			else{
			//else if(sizeofT % 2 == 0 && sizeofT <= TSize){
				logAdd_New(sumEven, result_F);
				//LOGADD_NEW(sumEven, result_F);
			}
			
		}
		
	}
} // end rec_compute_sum_F_S_T_TSize()




void compute_sum_F_S_T_TSize(int nh, int S, int ones, int * SSet, double * h, double ** alpha, double & sumOdd, double & sumEven, int TSize){
	int * TSet = new int[ones];
	rec_compute_sum_F_S_T_TSize(0, nh, S, ones, SSet, 0, 0, TSet, h, alpha, sumOdd, sumEven, TSize);
	delete[] TSet;
	
	
}//end void compute_sum_F_S_T_TSize()



void sub_h_Fast_TSize(int d, int nh, int S, int ones, int * SSet, double * h, double ** alpha, int TSize){
	if(d < nh){
		
		sub_h_Fast_TSize(d+1, nh, S, ones, SSet, h, alpha, TSize);
		
		SSet[ones] = d;
		sub_h_Fast_TSize(d+1, nh, S|(1<<d), (ones+1), SSet, h, alpha, TSize);
	}
	else{
		if(S == 0){
			h[0] = 0; //log(1) = 0;
			//cout << "h[0] = 0" << endl;
			//cout << endl;
			return;

		}
		
		double sumOdd = MARK;
		double sumEven = MARK;
		
		
		compute_sum_F_S_T_TSize(nh, S, ones, SSet, h, alpha, sumOdd, sumEven, TSize);
		
		//if sumEven has been assigned some value instead of MARK
		if(sumEven != MARK){
			logMinus_New_NonnegRes(sumOdd, sumEven);
			//logMinus_New(sumOdd, sumEven);
			//LOGMINUS_NEW(sumOdd, sumEven);
		}
						
		h[S] = sumOdd;
		//cout << "h[" << S << "] = "  << h[S] << endl;
		//cout << endl;

	}
} // end sub_h_Fast_TSize()



void compute_h_Fast_TSize(int nh, double * h, double ** alpha, int TSize){
	int * SSet = new int[nh];
	
	sub_h_Fast_TSize(0, nh, 0, 0, SSet, h, alpha, TSize);
	
	delete [] SSet;
	
} // end compute_h_Fast_TSize()







//HR: Compute: \sum_{k=1}^{|S|} (-1)^{k+1} \sum{T \subset S, |T|=k} RF(S,T)
void compute_sum_RF_S_T_Iter(int nh, int S, double * RR, double ** alpha, double & sumOdd, double & sumEven){
	
	for(int T = 1; T <= S; T++){
		//if T is the subset of S
		if((S|T) == S){
			
			int TSetLen = binIntLen(T);
			int * TSet = new int[TSetLen];
			binIntToDecIntArr(T, TSetLen, TSet);

			double result_RF = compute_RF(nh, S, T, TSetLen, TSet, RR, alpha);
			//sizeofT is odd
			if(TSetLen % 2 == 1){				
				logAdd_New(sumOdd, result_RF);
				//LOGADD_NEW(sumOdd, result_RF);
			}
			else{
				logAdd_New(sumEven, result_RF);
				//LOGADD_NEW(sumEven, result_RF);
			}
			
			delete[] TSet;
		}

	}

}

//HR:
//	Compute RR

void compute_RR_Iter(int nh, double * RR, double ** alpha){
	
	
	RR[0] = 0; //log(1) = 0;
	
	for(int S = 1; S < (1<<nh); S++){
		AA_S = new double [S+1];
		compute_AA_S_Iter(nh, S, AA_S, alpha);		
		
		double sumOdd = MARK;
		double sumEven = MARK;
				
		compute_sum_RF_S_T_Iter(nh, S, RR, alpha, sumOdd, sumEven);
		
		//if sumEven has been assigned some value instead of MARK
		if(sumEven != MARK){
			logMinus_New(sumOdd, sumEven);
			//LOGMINUS_NEW(sumOdd, sumEven);
		}
	
		RR[S] = sumOdd;
		//cout << "RR[" << S << "] = "  << h[S] << endl;
		//cout << endl;
		
		//HR
		delete[] AA_S;
		
	}
	
}


	
	
void rec_compute_Eta_U(
	int d,
	int nh,
	int TMax,
	int TMaxSetLen, //leng of SSet
    int * TMaxSet,
    int U,
    int T,
    int sizeofT,
    int * TSet,
    double * Eta_U,
    double ** alpha
){

	if(d < nh){	
		rec_compute_Eta_U(d+1, nh, TMax, TMaxSetLen, TMaxSet, U, T, sizeofT, TSet, Eta_U, alpha);
		
		if(((TMax>>d) & 1) == 1){
			TSet[sizeofT] = d;
			rec_compute_Eta_U(d+1, nh, TMax, TMaxSetLen, TMaxSet, U, (T | (1<<d)), (sizeofT + 1), TSet, Eta_U, alpha);
		}
	}
	else{
		if(T == 0){
			//cout << "\nT == 0" << endl;
			
			Eta_U[T] = 0;
			
			//cout << "AA_S[" << T << "]=" << AA_S[T] << endl;
			//cout << endl;
		}
		else{
	
			Eta_U[T] = alpha[TSet[0]][U] + Eta_U[T - (1<<TSet[0])];
			
		}
		
	}

}



void compute_Eta_U(
	int nh,
	int TMax,
	int TMaxSetLen, 
    int * TMaxSet,
    int U,
    double * Eta_U,
    double ** alpha
){
	int * TSet = new int[TMaxSetLen];
	rec_compute_Eta_U(0, nh, TMax, TMaxSetLen, TMaxSet, U, 0, 0, TSet, Eta_U, alpha);
	delete[] TSet;
	
}



void rec_compute_K_v(
	int d,
	int nh,
	int TMax,
	int TMaxSetLen, //leng of SSet
    int * TMaxSet,
    int U,
    int T,
    int sizeofT,
    int * TSet,
    double * Eta_U,
    double * K_v,
    double ** alpha,
    double & sumOdd,
    double & sumEven
){
	if(d < nh){	
		rec_compute_K_v(d+1, nh, TMax, TMaxSetLen, TMaxSet, U, T, sizeofT, TSet, Eta_U, K_v, alpha, sumOdd, sumEven);
		
		if(((TMax>>d) & 1) == 1){
			TSet[sizeofT] = d;
			rec_compute_K_v(d+1, nh, TMax, TMaxSetLen, TMaxSet, U, (T | (1<<d)), (sizeofT + 1), TSet, Eta_U, K_v, alpha, sumOdd, sumEven);
		}
	}
	else{
			double abs_result = RRFunc[TMax - T] + Eta_U[T];
			//sizeofT is odd
			if(sizeofT % 2 == 1){				
				logAdd_New(sumOdd, abs_result);
				//LOGADD_NEW(sumOdd, abs_result);
			}
			else{
				logAdd_New(sumEven, abs_result);
				//LOGADD_NEW(sumEven, abs_result);
			}
		
	}	
	
}



void compute_K_v(
	int nh,
	double * K_v,
	int v, 
	int U,
	int USetLen, 
    int * USet,
    double * Eta_U,
    double ** alpha
)


{
	int TMax = (1 << nh) - 1 - (1 << v) - U;
	int TMaxSetLen = nh - 1 - USetLen;
	int * TMaxSet = new int[TMaxSetLen];
	binIntToDecIntArr(TMax, TMaxSetLen, TMaxSet);
	
	Eta_U = new double [TMax+1];
	compute_Eta_U(nh, TMax, TMaxSetLen, TMaxSet, U, Eta_U, alpha);
	
//	totalOfMixEta << "\nEta_U U == " << U << ":" << endl;
//	for(int index = 0; index < TMax+1; index++){
//		totalOfMixEta << "Eta_U[" << index << "] = " << Eta_U[index] << endl;
//	}

	
	double result_K_v;
	
	double sumOdd = MARK;
	double sumEven = MARK;
	
	int * TSet = new int[TMaxSetLen];
	rec_compute_K_v(0, nh, TMax, TMaxSetLen, TMaxSet, U, 0, 0, TSet, Eta_U, K_v, alpha, sumOdd, sumEven);
	
	

	if(sumOdd == MARK){
		result_K_v = sumEven;	
	}
	else{

		if(sumEven > sumOdd){
			logMinus_New(sumEven, sumOdd);

			result_K_v = sumEven;
		}
		else if(sumEven == sumOdd){
			result_K_v = LOG_ZERO;
		}
		else{
			logMinus_New(sumOdd, sumEven);

			result_K_v = - sumOdd;

		}
	}
	K_v[U] = result_K_v;	
	
	delete[] TSet;
	
	delete[] Eta_U;			
	

}



void init_all_K_v(int d, int S, double * K_v, double value, int ones, int nh){

		//cerr<<" d: "<<d<<endl; 
	
	if (d < nh){
		init_all_K_v(d+1, S, K_v, value, ones, nh);
			
		init_all_K_v(d+1, S | (1 << d), K_v, value, ones+1, nh);
	}
	else{
		K_v[S] = value;
		
	}
}
	

void rec_compute_all_K_v(
	int d,
	int nh,
	int UMax,
	int UMaxSetLen, //leng of SSet
    int * UMaxSet,
    int U,
    int USetLen,
    int * USet,
    double * K_v,
    int v,
    double ** alpha
){
	if(d < nh){	
		rec_compute_all_K_v(d+1, nh, UMax, UMaxSetLen, UMaxSet, U, USetLen, USet, K_v, v, alpha);
		
		if(((UMax>>d) & 1) == 1){
			USet[USetLen] = d;
			rec_compute_all_K_v(d+1, nh, UMax, UMaxSetLen, UMaxSet, (U | (1<<d)), (USetLen + 1), USet, K_v, v, alpha);
		}
	}
	else{			
		compute_K_v(nh, K_v, v, U, USetLen, USet, Eta_U, alpha);
	
	}	
	
}



void compute_all_K_v(
	int nh,
	double * K_v,
	int v, 
    double ** alpha
){
	
	init_all_K_v(0, 0, K_v, LOG_ZERO, 0, nh);
	
	int UMax = (1 << nh) - 1 - (1 << v);
	int UMaxSetLen = nh - 1;
	int * UMaxSet = new int[UMaxSetLen];
	binIntToDecIntArr(UMax, UMaxSetLen, UMaxSet);
	
	int * USet = new int[UMaxSetLen];
	rec_compute_all_K_v(0, nh, UMax, UMaxSetLen, UMaxSet, 0, 0, USet, K_v, v, alpha);
	
	delete[] USet;
	
}



		
		

#endif


