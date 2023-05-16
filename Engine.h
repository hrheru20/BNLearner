#ifndef ENGINE_H
#define ENGINE_H

#include<stdio.h>
#include<stdlib.h>
#include <assert.h>


#include <sys/time.h>
#include <time.h>

#include <list>
#include <algorithm>
using namespace std;


#include"Arguments.h"
#include"Model.h"



#define MARK 9.0e99


#define LOGADD(la,lb) if(la>8e99||lb-la>100) la=lb;else la+=log(1+exp(lb-la));

#define LOG_ZERO -1.0e101


#include"UpdateHR.h"
#include"UpdateHR3.h"


using namespace std;



//ofstream totalOf2("debugEngine.txt");
//ofstream totalOf3("debugEngineUpdated.txt");
//ofstream totalOfAB("debugEngineAB.txt");
//ofstream totalOfMix("debugEngineMix.txt");




//----------------
// Fast Mobius transform routines.
void sub_fumt(int j, int d, int S, int overlap, double *t, double *s, int n, int k){
	if (d < n && overlap < k){
		sub_fumt(j, d+1, S, overlap, t, s, n, k);

		S |= (1 << d);
		overlap += (d >= j+1);
		sub_fumt(j, d+1, S, overlap, t, s, n, k);
	}

	else{
		s[S] = MARK;

		//jinS: j is in S, i.e. S_j == 1, i.e. T_j == 1
		int jinS = ((S >> j) & 1);
		if (overlap + jinS <= k) s[S] = t[S];

		if (jinS) logAdd_New(s[S], t[S - (1 << j)]);
	}
}

// Fast upward Mobius transform:
void fumt(double *t, double *s, int n, int k){
	double *tmp, *sprev = new double[1 << n];
	for (int T = 0; T < (1<<n); T ++)
		sprev[T] = t[T];
	for (int j = 0; j < n; j ++){
		sub_fumt(j, 0, 0, 0, sprev, s, n, k);
		tmp = sprev;
		sprev = s;
		s = tmp;
	}
	if (n % 2 == 0){
		for (int S = 0; S < (1<<n); S ++){
			s[S] = sprev[S];
		}
	}
	else{
		tmp = sprev;
		sprev = s;
		s = tmp;
	}
	delete [] sprev;
}
void test_fumt(){
	int n = 4;
	double t[16], s[16];
	for (int S = 0; S < 16; S ++){
		t[S] = log(S+1);
	}
	fumt(t, s, n, 2);
	for (int S = 0; S < 16; S ++){
		cerr<<" s["<<S<<"] = "<<exp(s[S])<<endl;
	}
}



void sub_fdmt(int j, int d, int T, int overlap, double *s, double *t, int n, int k){
	//cerr<<" T = "<<T<<":"; print_set(cerr, T);
	//cerr<<" overlap = "<<overlap<<endl;
	if (d < n && overlap <= k){
		sub_fdmt(j, d+1, T, overlap, s, t, n, k);


		T |= (1 << d);
		overlap += (d <= j);
		sub_fdmt(j, d+1, T, overlap, s, t, n, k);
	}

	// d == n || overlap > k
	else if (d == n){
		t[T] = s[T];

		int jinT = ((T >> j) & 1);

		if (jinT == 0)
			logAdd_New(t[T], s[T + (1 << j)]);
		if (overlap > k)
			t[T] = MARK;
	}

}


// Fast downward Mobius transform:
void fdmt(double *s, double *t, int n, int k){
	double *tmp, *tprev = new double[1 << n];
	for (int S = 0; S < (1<<n); S ++)
		tprev[S] = s[S];

	for (int j = 0; j < n; j ++){
		//for each T_{0:j}
		sub_fdmt(j, 0, 0, 0, tprev, t, n, k);
		tmp = t; t = tprev; tprev = tmp;
	}

	if (n % 2 == 0){
		for (int T = 0; T < (1<<n); T ++)
			t[T] = tprev[T];
	}

	else{
		tmp = t; t = tprev; tprev = tmp;
	}
	delete [] tprev;
}



void sub_fdmtLogAddComp(int j, int d, int T, int overlap, double *s, double *t, int n, int k){
	//cerr<<" T = "<<T<<":"; print_set(cerr, T);
	//cerr<<" overlap = "<<overlap<<endl;
	if (d < n && overlap <= k){
		sub_fdmtLogAddComp(j, d+1, T, overlap, s, t, n, k);
		T |= (1 << d);
		overlap += (d <= j);
		sub_fdmtLogAddComp(j, d+1, T, overlap, s, t, n, k);
	}
	else if (d == n){
		t[T] = s[T];
		int jinT = ((T >> j) & 1);
		if (jinT == 0) logAddComp(t[T], s[T + (1 << j)]);
		if (overlap > k) t[T] = MARK;
	}
}

void fdmtLogAddComp(double *s, double *t, int n, int k){
	double *tmp, *tprev = new double[1 << n];
	for (int S = 0; S < (1<<n); S ++) tprev[S] = s[S];
	for (int j = 0; j < n; j ++){

		sub_fdmtLogAddComp(j, 0, 0, 0, tprev, t, n, k);
		tmp = t; t = tprev; tprev = tmp;
	}
	if (n % 2 == 0){ for (int T = 0; T < (1<<n); T ++) t[T] = tprev[T];}
	else{ tmp = t; t = tprev; tprev = tmp;}
	delete [] tprev;
}



void sub_fdmtNoLog(int j, int d, int T, int overlap, double *s, double *t, int n, int k){
	//cerr<<" T = "<<T<<":"; print_set(cerr, T);
	//cerr<<" overlap = "<<overlap<<endl;
	if (d < n && overlap <= k){
		sub_fdmtNoLog(j, d+1, T, overlap, s, t, n, k);
		T |= (1 << d);
		overlap += (d <= j);
		sub_fdmtNoLog(j, d+1, T, overlap, s, t, n, k);
	}
	else if (d == n){
		t[T] = s[T];
		int jinT = ((T >> j) & 1);
		if (jinT == 0) t[T] += s[T + (1 << j)];
		if (overlap > k) t[T] = MARK;
	}
}


void fdmtNoLog(double *s, double *t, int n, int k){
	double *tmp, *tprev = new double[1 << n];
	for (int S = 0; S < (1<<n); S ++) tprev[S] = s[S];
	for (int j = 0; j < n; j ++){

		sub_fdmtNoLog(j, 0, 0, 0, tprev, t, n, k);
		tmp = t; t = tprev; tprev = tmp;
	}
	if (n % 2 == 0){ for (int T = 0; T < (1<<n); T ++) t[T] = tprev[T];}
	else{ tmp = t; t = tprev; tprev = tmp;}
	delete [] tprev;
}


void test_fdmt(){
	int n = 4;
	double t[16], s[16];
	for (int S = 0; S < 16; S ++){
		s[S] = log(S+1);
	}
	fdmt(s, t, n, 3);
	for (int T = 0; T < 16; T ++){
		cerr<<" t["<<T<<"] = "<<exp(t[T])<<endl;
	}
}

//----------------



// NOTE: Engine only calls a veru limited set of interface functions
// implemented in Model.
class Engine {
public:
	Engine(){}
	~Engine(){}

	void init(Model *m){
		model = m;

		//HR: restore
		//test();

	}
	void test(){
		vector<int> T;
		T.push_back(1); T.push_back(2);

		cerr<<" log p(x_0 | x_1,2): "<<model->log_lcp(0, T)<<endl;

		for (int h = 0; h < model->num_layers(); h ++){
			int *Vh, *Vu, *Vl;
			int nh, nu, nl;
			model->layer(h, &Vh, &nh);
			model->upper_layers(h, &Vu, &nu);
			model->lower_layers(h, &Vl, &nl);
			cerr<<"Layer "<<h<<":";
			cerr<<endl<<" layer: ";
			print_nodes(cerr, Vh, nh);
			cerr<<endl<<" upper: ";
			print_nodes(cerr, Vu, nu);
			cerr<<endl<<" lower: ";
			print_nodes(cerr, Vl, nl);
			cerr<<endl;
		}

		//test_fumt();
		test_fdmt();
		exit(1);

	}

	// Computes all edge probabilities.
	void compute_edge_probabilities(){
//		int l = model->num_layers();
//		totalOf2 << "compute_edge_probabilities()" << endl;
//		totalOf2 << "l=" << l << endl;

		// Separate computations for each layer.
//		for (int h = 0; h < l; h ++){
//			if (model->edges_within(h))
//				compute_edge_probabilities(h);
//		}

		//HR:
		//Arguments::option == 0 means the rebel method with new Dirichlet hyper-param
		if(Arguments::option == 0){
			compute_edge_probabilities(0);
		}

		//Arguments::option == 3 means the mix method in the UAI09 paper
		else if(Arguments::option == 3){
			compute_edge_probabilities_mixIndegree(0);

		}

		//Arguments::option == 8: DOS
		//Arguments::option == 9: DDS
		else if(Arguments::option == 8 || Arguments::option == 9){
			compute_edge_probabilities_sample1(0);
		}


		//Arguments::option == 11: IW-DDS
		else if(Arguments::option == 11){
			compute_edge_probabilities_sample3(0);
		}

		//Arguments::option == 12:  IW-DDS with huge N_o
		else if(Arguments::option == 12){
			compute_edge_probabilities_sample4(0);
		}

	}


	// Computes probabilities for edges pointing to the h-th layer.
	void compute_edge_probabilities(int h){
		cerr<<"\n REBEL Method:" << endl;
//		cerr<<" Compute edge probabilities for Layer "<<h<<":"<<endl;

		// Step 1: Compute beta[][];
		// Step 2: Compute alpha[][];
		int *Vh, *Vu; int nh, nu;
		model->layer(h, &Vh, &nh);
		model->upper_layers(h-1, &Vu, &nu);

//		totalOf2 << "nh=" << nh << endl;
//		for(int i = 0; i < nh; i++){
//			totalOf2 << "Vh[" << i << "] = " << Vh[i] << endl;
//		}
//		totalOf2 << "nu=" << nu << endl;
//		for(int i = 0; i < nu; i++){
//			totalOf2 << "Vu[" << i << "] = " << Vu[i] << endl;
//		}

		beta = new double*[nh];
		alpha = new double*[nh];

		k = model->max_indegree();

//		cerr<<" . "<<nh<<" nodes"<<endl;

		struct timeval startTime;
		struct timeval endTime;
		struct timezone myTimeZone;


		double betaTime;
		double alphaTime;
		double forwardTime;
		double backwardTime;
		double gammaTime;

		double step2BTime;

		gettimeofday(&startTime, &myTimeZone);


		if(Arguments::ADtree == 1){
			model->makeADTree();
		}

		for (int j = 0; j < nh; j ++){

			beta[j] = new double[1 << nh];
			//int i = Vh[j];
			compute_beta(j, Vh, nh, Vu, nu);

		}

		if(Arguments::ADtree == 1){
			model->freeADTree();
		}

		gettimeofday(&endTime, &myTimeZone);
		betaTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

		cerr<<" . Tables beta have been computed."<<endl;

		gettimeofday(&startTime, &myTimeZone);
		for (int j = 0; j < nh; j ++){

			alpha[j] = new double[1 << nh];
			compute_alpha(j, Vh, nh);


		}
		gettimeofday(&endTime, &myTimeZone);
		alphaTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

		cerr<<" . Tables alpha have been computed."<<endl;


		// Step 3: Compute g_forward[];
		gf = new double[1 << nh];

		gettimeofday(&startTime, &myTimeZone);

		compute_g_forward(nh);

		gettimeofday(&endTime, &myTimeZone);
		forwardTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

		//HR: restore
		//cerr<<" . gf[all] = "<<gf[(1<<nh)-1]<<endl;
		cerr<<" . Forward step has been computed."<<endl;


		// Step 4: Compute g_backward[];
		gb = new double[1 << nh];

		gettimeofday(&startTime, &myTimeZone);

		compute_g_backward(nh);

		gettimeofday(&endTime, &myTimeZone);
		backwardTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

		//HR: restore
		//cerr<<" . gb[all] = "<<gb[(1<<nh)-1]<<endl;
		cerr<<" . Backward step has been computed."<<endl;


		// Step 5: Loop over end-point nodes j.

		gammaTime = 0.0;
		step2BTime = 0.0;

		double *gamma = new double[1 << nh];
		double *gfb = new double[1 << nh];

		double ppDirectedFromTo[nh][nh];
		double ppUndirected[nh][nh];
		for(int i = 0; i < nh; i++){
			for(int j = 0; j < nh; j++){
				ppDirectedFromTo[i][j] = 0;
				ppUndirected[i][j] = 0;

			}
		}


		for (int j = 0; j < nh; j ++){
			// Compute gfb[].

			gettimeofday(&startTime, &myTimeZone);

			for (int S = 0; S < (1<<nh); S ++){

				//if j \in the binary repre of S
				if ((S>>j) & 1){
					//HR: can regard as log(0) = -inf; set gfb[S] to be -MARK so that FTDMT can be used directly to V instead of V - {j}
					gfb[S] = -MARK;
				}
				else{
					int cS = (1 << nh) - 1 - S - (1 << j);
					gfb[S] = gf[S] + gb[cS];
					//cerr<<" j = "<<j<<": ";
					//print_set(cerr, S);
					//cerr<<"; ";
					//print_set(cerr, cS);
					//cerr<<": "<<gfb[S];
					//cerr<<endl;
				}
			}
			// Compute gamma_j() := gamma[].
			fdmt(gfb, gamma, nh, k);

			gettimeofday(&endTime, &myTimeZone);
			gammaTime += (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;


			//cerr<<" gamma[0]: "<<gamma[0]<<endl;
			//cerr<<" gfb[0]: "<<gfb[0]<<endl;

			cerr<<" . Incoming edges for node "<<j<<":";

			gettimeofday(&startTime, &myTimeZone);

			for (int i = 0; i < nh; i ++){
				if (i == j) continue;
				cerr<<" "<<i;

				//HR: compute p(x,e), where e = i -> j
				double log_prob = eval_edge(i,  j, gamma, beta[j], nh, k);

				model->print_edge_prob(
					cout, Vh[i], Vh[j], exp(log_prob - gb[(1<<nh)-1])); //HR: - gb[(1<<nh)-1] mean /L(v)

				ppDirectedFromTo[i][j] = exp(log_prob - gb[(1<<nh)-1]);


			}
			gettimeofday(&endTime, &myTimeZone);
			step2BTime += (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;


			cerr<<endl;

		}


		cout << "\n\n\nTime Statistics" << endl;
		cout << "betaTime = " << betaTime <<  endl;
		cout << "alphaTime = " << alphaTime <<  endl;
		cout << "forwardTime = " << forwardTime <<  endl;
		cout << "backwardTime = " << backwardTime <<  endl;
		cout << "gammaTime = " << gammaTime <<  endl;

		cout << "step2BTime = " << step2BTime <<  endl;

		double totalTime = betaTime + alphaTime + forwardTime
						+ backwardTime + gammaTime + step2BTime;
		cout << "totalTime = " << totalTime <<  endl;


		cout << "\nExact_ppDirectedFromTo \t Exact_ppDirectedFromTo" << endl;
		for(int i = 0; i < nh; i++){
			for(int j = 0; j < nh; j++){
				cout << "ppDirectedFromTo[" << j << "][" << i << "] \t " << ppDirectedFromTo[j][i] << endl;

			}
		}
		cout << "\n\nPlain output:" << endl;
		cout << "Exact_ppDirectedFromTo" << endl;
		for(int i = 0; i < nh; i++){
			for(int j = 0; j < nh; j++){
				cout << ppDirectedFromTo[j][i] << endl;

			}
		}



//		cout << "\nExact_ppUndirected \t Exact_ppUndirected" << endl;
//		for(int i = 0; i < nh; i++){
//			for(int j = i + 1; j < nh; j++){
//				ppUndirected[i][j] = ppDirectedFromTo[i][j] + ppDirectedFromTo[j][i];
//				cout << "ppUndirected[" << i << "][" << j << "] \t " << ppUndirected[i][j] << endl;
//
//
//			}
//		}



		delete [] gamma;
		delete [] gfb;

		delete [] gf; delete [] gb;
		for (int j = 0; j < nh; j ++){
			delete [] beta[j]; delete [] alpha[j];
		}
		delete [] beta; delete [] alpha;
		cerr<<" Posterior of each edge now computed."<<endl;

	}//void compute_edge_probabilities(int h)




	//	  Used to implement the idea of DOS and DDS
	void compute_edge_probabilities_sample1(int h){
		//cerr<<"\nOrder-sampling Based Method (without bias correction) (DOS/DDS):" << endl;
		//Arguments::option == 8: DOS
		//Arguments::option == 9: DDS
		if(Arguments::option == 8){
			cerr<<"\n DOS method:" << endl;
		}
		else if(Arguments::option == 9){
			cerr<<"\n DDS method:" << endl;
		}
		//cerr<<" Compute edge probabilities for Layer "<<h<<":"<<endl;

		// Step 1: Compute beta[][];
		// Step 2: Compute alpha[][];
		int *Vh, *Vu; int nh, nu;
		model->layer(h, &Vh, &nh);
		model->upper_layers(h-1, &Vu, &nu);

//		totalOf2 << "nh=" << nh << endl;
//		for(int i = 0; i < nh; i++){
//			totalOf2 << "Vh[" << i << "] = " << Vh[i] << endl;
//		}
//		totalOf2 << "nu=" << nu << endl;
//		for(int i = 0; i < nu; i++){
//			totalOf2 << "Vu[" << i << "] = " << Vu[i] << endl;
//		}

		beta = new double*[nh];
		alpha = new double*[nh];

		k = model->max_indegree();

//		cerr<<" . "<<nh<<" nodes"<<endl;

		struct timeval startTime;
		struct timeval endTime;
		struct timezone myTimeZone;


		double betaTime;
		double alphaTime;
		double forwardTime;


		gettimeofday(&startTime, &myTimeZone);


		//If Arguments::ADtree == 1, user wants to use ADtree, then make the whole tree first
		if(Arguments::ADtree == 1){
			model->makeADTree();
		}

		for (int j = 0; j < nh; j ++){

			beta[j] = new double[1 << nh];
			//int i = Vh[j];
			compute_beta(j, Vh, nh, Vu, nu);

		}

		//If Arguments::ADtree == 1, then delete the whole tree
		if(Arguments::ADtree == 1){
			model->freeADTree();
		}

		gettimeofday(&endTime, &myTimeZone);
		betaTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//		cerr<<" . Tables beta are now ready."<<endl;

		gettimeofday(&startTime, &myTimeZone);
		for (int j = 0; j < nh; j ++){

			alpha[j] = new double[1 << nh];
			compute_alpha(j, Vh, nh);


		}
		gettimeofday(&endTime, &myTimeZone);
		alphaTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//		cerr<<" . Tables alpha are now ready."<<endl;


		// Step 3: Compute g_forward[];
		gf = new double[1 << nh];

		gettimeofday(&startTime, &myTimeZone);

		compute_g_forward(nh);

		gettimeofday(&endTime, &myTimeZone);
		forwardTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

		//HR: restore
		//cerr<<" . gf[all] = "<<gf[(1<<nh)-1]<<endl;

		cerr<<" \tDP step is done."<<endl;


		//HR: UAI_12 starts here
		double sampleOrderTime;


		//HR: Step 4.1
		//HR: sample orders
		//HR: Reset the random number generator with the system clock.
		srand(time(NULL));
		//set the same value for the testing
		//srand(100);


		gettimeofday(&startTime, &myTimeZone);

		const int num_of_orders = Arguments::num_of_sampled_orders;

		const int num_of_dags_per_order = 1;


		list<total_order> t_order_list1;

		total_order * t_order_temp_p;
		for(int order_ind = 0; order_ind < num_of_orders; order_ind++){
			t_order_temp_p = sample_one_order(nh);
			t_order_list1.push_back(* t_order_temp_p);

			delete t_order_temp_p;
		}
//		cout << "t_order_list1.size() = " << t_order_list1.size() << endl;

		gettimeofday(&endTime, &myTimeZone);
		sampleOrderTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

		cerr<<" \tOrder-Sampling is done."<<endl;


		//HR: new Step 4.2
		//	compute the posterior of each edge using the modular feature formula
		if(Arguments::option == 8){

			double computeEdgePosteriorTime;
			gettimeofday(&startTime, &myTimeZone);


			//the estimated posterior of each edge (n*(n-1) edges) is the average of each P(edge | < , D)
			double ** edge_posterior;
			edge_posterior = new double * [nh];
			for (int i = 0; i < nh; i++){
				edge_posterior[i] = new double[nh];
			}


			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					edge_posterior[i][j] = 0;
				}
			}

			for (list<total_order>::iterator it1 = t_order_list1.begin(); it1 != t_order_list1.end(); it1++){
				 add_edge_poster_given_order(edge_posterior, (*it1), nh);
			}

			//finally / num_of_orders
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					edge_posterior[i][j] = edge_posterior[i][j] / ((double) t_order_list1.size());
				}
			}

			gettimeofday(&endTime, &myTimeZone);
			computeEdgePosteriorTime = (endTime.tv_sec - startTime.tv_sec) +
			           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//			cerr<<" \tComputing the posterior of each edge is done."<<endl;


			cout << "\nedge_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					cout << "edge_posterior[" << j << "][" << i << "] = \t" << edge_posterior[j][i] << endl;
				}
			}

			cout << "\n\nPlain version:" << endl;
			cout << "edge_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					cout << edge_posterior[j][i] << endl;
				}
			}

			//write to the file
			ofstream of_poster("./results/edge_poster");;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					of_poster << edge_posterior[j][i] << endl;
				}
			}


			//free edge_posterior
			for(int i = 0; i < nh; i++){
				delete [] edge_posterior[i];
			}
			delete [] edge_posterior;


			cout << "\n\n\nTime Statistics" << endl;
			cout << "betaTime = " << betaTime <<  endl;
			cout << "alphaTime = " << alphaTime <<  endl;
			cout << "forwardTime = " << forwardTime <<  endl;

			cout << "sampleOrderTime = " << sampleOrderTime <<  endl;

			cout << "computeEdgePosteriorTime = " << computeEdgePosteriorTime <<  endl;


			double totalTime = betaTime + alphaTime + forwardTime
					+ sampleOrderTime + computeEdgePosteriorTime;

			cout << "totalTime = " << totalTime <<  endl;


			string str1;
			str1.assign("./results/total_time_");
			//+6 to delete the prefix "cases/"
			str1.append((Arguments::datafile + 6));
			FILE* file_time = fopen(str1.c_str(), "a");
			fprintf(file_time,"%f \n", totalTime);
			fclose(file_time);


			double DPTime = betaTime + alphaTime + forwardTime;
			cout << "DPTime = " << DPTime <<  endl;

			string str1_2;
			str1_2.assign("./results/DP_time_");
			//+6 to delete the prefix "cases/"
			str1_2.append((Arguments::datafile + 6));
			FILE* file_time_2 = fopen(str1_2.c_str(), "a");
			fprintf(file_time_2,"%f \n", DPTime);
			fclose(file_time_2);


//			string str1_3;
//			str1_3.assign("./results/beta_time_");
//			//+6 to delete the prefix "cases/"
//			str1_3.append((Arguments::datafile + 6));
//			FILE* file_time_3 = fopen(str1_3.c_str(), "a");
//			fprintf(file_time_3,"%f \n", betaTime);
//			fclose(file_time_3);


			string str1_4;
			str1_4.assign("./results/sampleOrder_time_");
			//+6 to delete the prefix "cases/"
			str1_4.append((Arguments::datafile + 6));
			FILE* file_time_4 = fopen(str1_4.c_str(), "a");
			fprintf(file_time_4,"%f \n", sampleOrderTime);
			fclose(file_time_4);


			//Finally
			//free the rest of the memory

			delete [] gf;

			for (int j = 0; j < nh; j ++){
				delete [] beta[j];
			}
			delete [] beta;
			for (int j = 0; j < nh; j ++){
				delete [] alpha[j];
			}
			delete [] alpha;
			cerr<<"\n Posterior of each edge now computed."<<endl;


		}//if(Arguments::option == 8)

		//HR: Step 4.2
		//	  compute the posterior of each edge by sampling DAGs
		//    DDS
		else if(Arguments::option == 9){
			//free gf[] before DAG sampling step to save the memory
			delete [] gf;

			double sampleDAGTime;

			gettimeofday(&startTime, &myTimeZone);

			//if(nh > 17){
				t_order_list1.sort();
//				cout << "*******************" << endl;
//				cout << "After sort():" << endl;
//				cout << "t_order_list1.size() = " << t_order_list1.size() << endl;
			//}


			//Only remember parents (parent set) instead of the whole family to save the space cost
			vector<int> * * * U_i_fam_vector_ptr;
			U_i_fam_vector_ptr = new vector<int> * * [nh];
			for(int i = 0; i < nh; i++){
				U_i_fam_vector_ptr[i] = new vector<int> * [(1<<nh)];
			}

			//initialize U_i_raw_fam_vector_ptr
			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					U_i_fam_vector_ptr[ch][fam] = NULL;
				}
			}

			double * * * U_i_fam_prob_cumul_vector_ptr;
			U_i_fam_prob_cumul_vector_ptr = new double * * [nh];
			for(int i = 0; i < nh; i++){
				U_i_fam_prob_cumul_vector_ptr[i] = new double * [(1<<nh)];
			}


			//initialize U_i_fam_prob_cumul_vector_ptr
			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					U_i_fam_prob_cumul_vector_ptr[ch][fam] = NULL;
				}
			}

			//the number of non-NULL U_i_fam_vector
			int U_i_fam_vector_setup_count = 0;

			//the number of cached families
			int family_setup_count = 0;

			//change int type into short type to save the memory
			short * * U_i_fam_vector_used_count;
			U_i_fam_vector_used_count = new short * [nh];
			for(int i = 0; i < nh; i++){
				U_i_fam_vector_used_count[i] = new short [(1<<nh)];
			}
			for (int ch = 0; ch < nh; ch++) {
				for (int fam = 0; fam < (1 << nh); fam++) {
					U_i_fam_vector_used_count[ch][fam] = 0;
				}
			}
			const int max_recorded_used_count = 10;


//			cout << "After initialize U_i_fam_vector_used_count[ch][fam] " << endl;

			list<dag> dag_list1;


			int order_num = 0;
			int numOfRecycle = 0;
//			cout << "order_num = " << order_num << endl;

			for (list<total_order>::iterator it1 = t_order_list1.begin(); it1 != t_order_list1.end(); it1++){
				order_num++;

				for(int ind = 0; ind < (*it1).length; ind++ ){
					int sigma_i = (*it1).torder[ind];
					int U_sigma_i = (*it1).torderBP[ind];
					if(U_i_fam_vector_ptr[sigma_i][U_sigma_i] == NULL){

						int numOfFam = getNumbOfFam(k, ind);

						U_i_fam_vector_ptr[sigma_i][U_sigma_i] = get_raw_family_vector(sigma_i, U_sigma_i, nh, k, numOfFam);
						assert(numOfFam == (int) (*U_i_fam_vector_ptr[sigma_i][U_sigma_i]).size());

						U_i_fam_prob_cumul_vector_ptr[sigma_i][U_sigma_i] = get_prob_cumul_array(U_i_fam_vector_ptr[sigma_i][U_sigma_i], alpha[sigma_i][U_sigma_i], sigma_i);
						U_i_fam_vector_setup_count++;

						if(U_i_fam_vector_used_count[sigma_i][U_sigma_i] < max_recorded_used_count){
							U_i_fam_vector_used_count[sigma_i][U_sigma_i]++;
						}

						family_setup_count += (*U_i_fam_vector_ptr[sigma_i][U_sigma_i]).size();


					}//if(U_i_fam_vector_ptr[sigma_i][U_sigma_i] == NULL)
					else {

						if(U_i_fam_vector_used_count[sigma_i][U_sigma_i] < max_recorded_used_count){
							U_i_fam_vector_used_count[sigma_i][U_sigma_i]++;
						}

					}

				}//for(int ind = 0; ind < (*it1).length; ind++ )

				if(order_num % 1000 == 0){
//					cout << "order_num = " << order_num << endl;
//					cout << "U_i_fam_vector_setup_count = " << U_i_fam_vector_setup_count << endl;
//					cout << "family_setup_count = " << family_setup_count << endl;

				}


				for(int dag_seq_per_order = 1; dag_seq_per_order <= num_of_dags_per_order; dag_seq_per_order++){

					dag samp_dag = sample_one_dag((*it1), U_i_fam_vector_ptr, U_i_fam_prob_cumul_vector_ptr);

					dag_list1.push_back(samp_dag);

				}

				const double Max_family_setup_num = 3.0e7; //for 2 G memory

				const double Min_released_family_setup_num = 1.9e7; //for 2 G memory

				//need to delete some cache if it occupies too much memory (when family_setup_count >= Max_family_setup_num)
				if(family_setup_count > Max_family_setup_num ){
//					cout << "\norder_num = " << order_num << endl;
					numOfRecycle++;

					recycleMemory(
							U_i_fam_vector_ptr,
							U_i_fam_prob_cumul_vector_ptr,
							U_i_fam_vector_setup_count,
							family_setup_count,
							U_i_fam_vector_used_count,
							Max_family_setup_num,
							Min_released_family_setup_num,
							nh);

				}//if(family_setup_count > Max_family_setup_num )


			}//for (vector<total_order>::iterator it1 = t_order_list1.begin(); ... )



			gettimeofday(&endTime, &myTimeZone);
			sampleDAGTime = (endTime.tv_sec - startTime.tv_sec) +
					(endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

			cerr<<" \tDAG-Sampling is done."<<endl;

//			cout << "numOfRecycle = " << numOfRecycle << endl;


			//get the estimated posterior of each edge (n*(n-1) edges) by weighting the number of dags in dag_list1
			double computeEdgePosteriorTime;
			gettimeofday(&startTime, &myTimeZone);

			double ** edge_posterior;
			edge_posterior = new double * [nh];
			for (int i = 0; i < nh; i++){
				edge_posterior[i] = new double[nh];
			}

//			cout << "dag_list1.size() = " << dag_list1.size() << endl;
			get_edge_poster_by_num_of_dags(edge_posterior, nh, dag_list1);

			gettimeofday(&endTime, &myTimeZone);
			computeEdgePosteriorTime = (endTime.tv_sec - startTime.tv_sec) +
			           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//			cerr<<" Computing the posterior of each edge are now done."<<endl;


			cout << "\nedge_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					cout << "edge_posterior[" << j << "][" << i << "] = \t" << edge_posterior[j][i] << endl;
				}
			}

			cout << "\n\nPlain version:" << endl;
			cout << "edge_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					cout << edge_posterior[j][i] << endl;
				}
			}

			//write to the file
			ofstream of_poster("./results/edge_poster");;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					of_poster << edge_posterior[j][i] << endl;
				}
			}



			//free the memory
			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					if(U_i_fam_vector_ptr[ch][fam] != NULL){
						(*U_i_fam_vector_ptr[ch][fam]).clear();
						delete U_i_fam_vector_ptr[ch][fam];
					}

				}
			}

			for(int i = 0; i < nh; i++){
				delete [] U_i_fam_vector_ptr[i];
			}
			delete [] U_i_fam_vector_ptr;


			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					if(U_i_fam_prob_cumul_vector_ptr[ch][fam] != NULL){
						delete [] U_i_fam_prob_cumul_vector_ptr[ch][fam];
					}

				}
			}

			for(int i = 0; i < nh; i++){
				delete [] U_i_fam_prob_cumul_vector_ptr[i];
			}
			delete [] U_i_fam_prob_cumul_vector_ptr;


			for(int i = 0; i < nh; i++){
				delete [] edge_posterior[i];
			}
			delete [] edge_posterior;


			cout << "\n\n\nTime Statistics" << endl;
			cout << "betaTime = " << betaTime <<  endl;
			cout << "alphaTime = " << alphaTime <<  endl;
			cout << "forwardTime = " << forwardTime <<  endl;

			cout << "sampleOrderTime = " << sampleOrderTime <<  endl;
			cout << "sampleDAGTime = " << sampleDAGTime <<  endl;
			cout << "computeEdgePosteriorTime = " << computeEdgePosteriorTime <<  endl;

			double totalTime = betaTime + alphaTime + forwardTime
					+ sampleOrderTime + sampleDAGTime + computeEdgePosteriorTime;

			cout << "totalTime = " << totalTime <<  endl;


			string str1;
			str1.assign("./results/total_time_");

			str1.append((Arguments::datafile + 6));
			FILE* file_time = fopen(str1.c_str(), "a");
			fprintf(file_time,"%f \n", totalTime);
			fclose(file_time);


			double DPTime = betaTime + alphaTime + forwardTime;
			cout << "DPTime = " << DPTime <<  endl;

			string str1_2;
			str1_2.assign("./results/DP_time_");

			str1_2.append((Arguments::datafile + 6));
			FILE* file_time_2 = fopen(str1_2.c_str(), "a");
			fprintf(file_time_2,"%f \n", DPTime);
			fclose(file_time_2);


//			string str1_3;
//			str1_3.assign("./results/beta_time_");
//
//			str1_3.append((Arguments::datafile + 6));
//			FILE* file_time_3 = fopen(str1_3.c_str(), "a");
//			fprintf(file_time_3,"%f \n", betaTime);
//			fclose(file_time_3);


			string str1_4;
			str1_4.assign("./results/sampleOrder_time_");

			str1_4.append((Arguments::datafile + 6));
			FILE* file_time_4 = fopen(str1_4.c_str(), "a");
			fprintf(file_time_4,"%f \n", sampleOrderTime);
			fclose(file_time_4);


			string str1_5;
			str1_5.assign("./results/sampleDAG_time_");

			str1_5.append((Arguments::datafile + 6));
			FILE* file_time_5 = fopen(str1_5.c_str(), "a");
			fprintf(file_time_5,"%f \n", sampleDAGTime);
			fclose(file_time_5);


			//Finally
			//free the rest of the memory

			for (int j = 0; j < nh; j ++){
				delete [] alpha[j];
			}
			delete [] alpha;
			for (int j = 0; j < nh; j ++){
				delete [] beta[j];
			}
			delete [] beta;
			cerr<<"\n Posterior of each edge now computed."<<endl;


		}//else if(Arguments::option == 9)


	}//void compute_edge_probabilities_sample1(int h)




	//	  used to implement IW-DDS
	void compute_edge_probabilities_sample3(int h){
		cerr << "\n IW-DDS method:" << endl;

//		cerr<<" Compute edge probabilities for Layer "<<h<<":"<<endl;

		// Step 1: Compute beta[][];
		// Step 2: Compute alpha[][];
		int *Vh, *Vu; int nh, nu;
		model->layer(h, &Vh, &nh);
		model->upper_layers(h-1, &Vu, &nu);

//		totalOf2 << "nh=" << nh << endl;
//		for(int i = 0; i < nh; i++){
//			totalOf2 << "Vh[" << i << "] = " << Vh[i] << endl;
//		}
//		totalOf2 << "nu=" << nu << endl;
//		for(int i = 0; i < nu; i++){
//			totalOf2 << "Vu[" << i << "] = " << Vu[i] << endl;
//		}

		beta = new double*[nh];
		alpha = new double*[nh];

		k = model->max_indegree();

//		cerr<<" . "<<nh<<" nodes"<<endl;

		struct timeval startTime;
		struct timeval endTime;
		struct timezone myTimeZone;

		double betaTime;
		double alphaTime;
		double forwardTime;


		gettimeofday(&startTime, &myTimeZone);


		//If Arguments::ADtree == 1, user wants to use ADtree, then make the whole tree first
		if(Arguments::ADtree == 1){
			model->makeADTree();
		}

		for (int j = 0; j < nh; j ++){

			beta[j] = new double[1 << nh];
			compute_beta(j, Vh, nh, Vu, nu);

		}

		//If Arguments::ADtree == 1, then delete the whole tree
		if(Arguments::ADtree == 1){
			model->freeADTree();
		}

		gettimeofday(&endTime, &myTimeZone);
		betaTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//		cerr<<" . Tables beta are now ready."<<endl;

		gettimeofday(&startTime, &myTimeZone);
		for (int j = 0; j < nh; j ++){

			alpha[j] = new double[1 << nh];
			compute_alpha(j, Vh, nh);


		}
		gettimeofday(&endTime, &myTimeZone);
		alphaTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//		cerr<<" . Tables alpha are now ready."<<endl;


		// Step 3: Compute g_forward[];
		gf = new double[1 << nh];

		gettimeofday(&startTime, &myTimeZone);

		compute_g_forward(nh);

		gettimeofday(&endTime, &myTimeZone);
		forwardTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

		//cerr<<" . gf[all] = "<<gf[(1<<nh)-1]<<endl;
		cerr<<" \tDP step is done."<<endl;


		double sampleOrderTime;


		//HR: Step 4.1
		//HR: sample orders
		//HR: Reset the random number generator with the system clock.
		srand(time(NULL));
		//set the same value for the testing
		//srand(100);


		gettimeofday(&startTime, &myTimeZone);

		const int num_of_orders = Arguments::num_of_sampled_orders;


		list<total_order> t_order_list1;
		list<total_order> t_order_list2;
		total_order * t_order_temp_p;
		for(int order_ind = 0; order_ind < num_of_orders; order_ind++){
			t_order_temp_p = sample_one_order(nh);
			t_order_list1.push_back(* t_order_temp_p);
			t_order_list2.push_back(* t_order_temp_p);

			delete t_order_temp_p;
		}
//		cout << "t_order_list1.size() = " << t_order_list1.size() << endl;


		t_order_list2.sort();
//		cout << "*******************" << endl;
//		cout << "After sort():" << endl;
//		cout << "t_order_list2:" << endl;
//		cout << "t_order_list2.size() = " << t_order_list2.size() << endl;


		unique_HR(t_order_list2);
//		cout << "*******************" << endl;
//		cout << "See after uniqueHR(t_order_list2):" << endl;
//		cout << "t_order_list2:" << endl;
//		cout << "t_order_list2.size() = " << t_order_list2.size() << endl;

//		double sum_order_scores  = get_sum_order_scores(t_order_list2);
//		cout << "exp(sum_order_scores) = " << exp(sum_order_scores) << endl;

		t_order_list2.clear();


		gettimeofday(&endTime, &myTimeZone);
		sampleOrderTime = (endTime.tv_sec - startTime.tv_sec)
				+ (endTime.tv_usec - startTime.tv_usec) / 1000000.0;

		cerr << " \tOrder-Sampling is done." << endl;


		//HR: new Step 4.2
		//	compute the posterior of each edge by sampling k=1 DAGs per order
//		if(Arguments::option == 11){
		    //free gf[] before DAG sampling step to save the memory
//			const double gf_all = gf[(1<<nh)-1];
			delete [] gf;

			double sampleDAGTime;

			gettimeofday(&startTime, &myTimeZone);

//			if(nh > 17){
				t_order_list1.sort(); //
//				cout << "*******************" << endl;
//				cout << "After sort():" << endl;
//				cout << "t_order_list1.size() = " << t_order_list1.size() << endl;

//			}


			//Only remember parents (parent set) instead of the whole family to save the space cost
			vector<int> * * * U_i_fam_vector_ptr;
			U_i_fam_vector_ptr = new vector<int> * * [nh];
			for(int i = 0; i < nh; i++){
				U_i_fam_vector_ptr[i] = new vector<int> * [(1<<nh)];
			}

			//initialize U_i_raw_fam_vector_ptr
			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					U_i_fam_vector_ptr[ch][fam] = NULL;
				}
			}


			double * * * U_i_fam_prob_cumul_vector_ptr;
			U_i_fam_prob_cumul_vector_ptr = new double * * [nh];
			for(int i = 0; i < nh; i++){
				U_i_fam_prob_cumul_vector_ptr[i] = new double * [(1<<nh)];
			}

			//initialize U_i_fam_prob_cumul_vector_ptr
			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					U_i_fam_prob_cumul_vector_ptr[ch][fam] = NULL;
				}
			}
			//


			//the number of non-NULL U_i_fam_vector
			int U_i_fam_vector_setup_count = 0;

			//the number of cached families
			int family_setup_count = 0;


			//change int type into short type to save the memory
			short * * U_i_fam_vector_used_count;
			U_i_fam_vector_used_count = new short * [nh];
			for(int i = 0; i < nh; i++){
				U_i_fam_vector_used_count[i] = new short [(1<<nh)];
			}
			for (int ch = 0; ch < nh; ch++) {
				for (int fam = 0; fam < (1 << nh); fam++) {
					U_i_fam_vector_used_count[ch][fam] = 0;
				}
			}
			const int max_recorded_used_count = 10;


			list<dag> dag_list1;

			int total_num_sampled_dags = 0;

			//HR: Use the on-fly version to avoid using too much memory
			//for each sampled total order <
			int order_num = 0;
			int numOfRecycle = 0;

			for (list<total_order>::iterator it1 = t_order_list1.begin(); it1 != t_order_list1.end(); it1++){
				order_num++;

				for(int ind = 0; ind < (*it1).length; ind++ ){
					int sigma_i = (*it1).torder[ind];
					int U_sigma_i = (*it1).torderBP[ind];
					if(U_i_fam_vector_ptr[sigma_i][U_sigma_i] == NULL){

						int numOfFam = getNumbOfFam(k, ind);

						U_i_fam_vector_ptr[sigma_i][U_sigma_i] = get_raw_family_vector(sigma_i, U_sigma_i, nh, k, numOfFam);
						assert(numOfFam == (int) (*U_i_fam_vector_ptr[sigma_i][U_sigma_i]).size());

						U_i_fam_prob_cumul_vector_ptr[sigma_i][U_sigma_i] = get_prob_cumul_array(U_i_fam_vector_ptr[sigma_i][U_sigma_i], alpha[sigma_i][U_sigma_i], sigma_i);

						U_i_fam_vector_setup_count++;

						if(U_i_fam_vector_used_count[sigma_i][U_sigma_i] < max_recorded_used_count){
							U_i_fam_vector_used_count[sigma_i][U_sigma_i]++;
						}

						family_setup_count += (*U_i_fam_vector_ptr[sigma_i][U_sigma_i]).size();


					}//if(U_i_fam_vector_ptr[sigma_i][U_sigma_i] == NULL)
					else {

						if(U_i_fam_vector_used_count[sigma_i][U_sigma_i] < max_recorded_used_count){
							U_i_fam_vector_used_count[sigma_i][U_sigma_i]++;
						}

					}

				}//for(int ind = 0; ind < (*it1).length; ind++ )




				const int req_num_of_sampled_dags_per_order = 1;

				int num_sampled_dags_per_order = 0;

				//order_dag_list1 stored all the unique sampled DAGs per order in the descending order
				list<dag> order_dag_list1;


				while (num_sampled_dags_per_order < req_num_of_sampled_dags_per_order) {
					//dag samp_dag = sample_one_dag((*it1), U_i_fam_vector_ptr);
					dag samp_dag = sample_one_dag((*it1), U_i_fam_vector_ptr, U_i_fam_prob_cumul_vector_ptr);
					order_dag_list1.push_back(samp_dag);
					num_sampled_dags_per_order++;

				}//while()


				dag_list1.splice(dag_list1.begin(), order_dag_list1);


				total_num_sampled_dags += num_sampled_dags_per_order;

				if(order_num % 20000 == 1){
//					cout << "\norder_num = " << order_num << endl;
//					cout << "total_num_sampled_dags = " << total_num_sampled_dags << endl;
//
//					cout << "dag_list1.size() = " << dag_list1.size() << endl;
//					cout << "U_i_fam_vector_setup_count = " << U_i_fam_vector_setup_count << endl;
//					cout << "family_setup_count = " << family_setup_count << endl;
				}


				if(order_num % 200000 == 1){
//					cout << "\norder_num = " << order_num << endl;
//					cout << "dag_list1.sort(descend_order_comparison) " << endl;
//					cout << "unique_HR(dag_list1) " << endl;
					dag_list1.sort(descend_order_comparison);
					unique_HR(dag_list1);

				}



				const double Max_family_setup_num = 3.0e7; //for 2 G memory

				const double Min_released_family_setup_num = 1.9e7; //for 2 G memory

				//need to delete some cache if it occupies too much memory (when family_setup_count >= Max_family_setup_num)
				if(family_setup_count > Max_family_setup_num ){
//					cout << "\norder_num = " << order_num << endl;
					numOfRecycle++;

					recycleMemory(
							U_i_fam_vector_ptr,
							U_i_fam_prob_cumul_vector_ptr,
							U_i_fam_vector_setup_count,
							family_setup_count,
							U_i_fam_vector_used_count,
							Max_family_setup_num,
							Min_released_family_setup_num,
							nh);


				}//if(family_setup_count > Max_family_setup_num )


			}//for (vector<total_order>::iterator it1 = t_order_list1.begin(); ... )

//			cout << "\ntotal_num_sampled_dags = " << total_num_sampled_dags << endl;

			dag_list1.sort(descend_order_comparison);
			unique_HR(dag_list1);


			gettimeofday(&endTime, &myTimeZone);
			sampleDAGTime = (endTime.tv_sec - startTime.tv_sec) +
					(endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

			cerr<<" \tDAG-Sampling and Bias-Correction are done."<<endl;

//			cout << "numOfRecycle = " << numOfRecycle << endl;

//			cout <<"\nFinally:" << endl;
//			cout << "dag_list1.size() = " << dag_list1.size() << endl;


			//get the estimated posterior of each edge (n*(n-1) edges) by weighting the number of dags in dag_list1
			double computeEdgePosteriorTime;
			gettimeofday(&startTime, &myTimeZone);

			double ** edge_posterior;
			edge_posterior = new double * [nh];
			for (int i = 0; i < nh; i++){
				edge_posterior[i] = new double[nh];
			}


			double sum_posterior_dag_list = LOG_ZERO;
			double temp_dag_scores;
			for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){

				temp_dag_scores = (*it1).scores;
				logAdd_New(sum_posterior_dag_list, temp_dag_scores);
			}

//			printf("sum_posterior_dag_list = %18.6f\n", sum_posterior_dag_list);
//
//			printf("P(D) = gf[all] = %18.6f\n", gf_all);
//
//			printf("sum_i P(G_i|D) = %18.6f\n", exp(sum_posterior_dag_list - gf_all));



			//get_edge_poster_by_num_of_dags(edge_posterior, nh, dag_list1);
			get_edge_poster_by_posterir_of_dags(edge_posterior, nh, dag_list1);

			gettimeofday(&endTime, &myTimeZone);
			computeEdgePosteriorTime = (endTime.tv_sec - startTime.tv_sec) +
			           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//			cerr<<"Compute the posterior of each edge are now done."<<endl;


			//note the order of the double loop, consistent with the column order of Matlab
			cout << "\nedge_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					cout << "edge_posterior[" << j << "][" << i << "] = \t" << edge_posterior[j][i] << endl;
				}
			}

			cout << "\n\nPlain version:" << endl;
			cout << "edge_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					cout << edge_posterior[j][i] << endl;
				}
			}

			ofstream of_edge_poster("./results/edge_poster");;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					of_edge_poster << edge_posterior[j][i] << endl;
				}
			}



			//get the estimated posterior of each path (n*(n-1) paths) by weighting the number of dags in dag_list1
			double computePathPosteriorTime;
			gettimeofday(&startTime, &myTimeZone);

			double ** path_posterior;
			path_posterior = new double * [nh];
			for (int i = 0; i < nh; i++){
				path_posterior[i] = new double[nh];
			}

			get_path_poster_by_posterir_of_dags(path_posterior, nh, dag_list1);

			gettimeofday(&endTime, &myTimeZone);
			computePathPosteriorTime = (endTime.tv_sec - startTime.tv_sec) +
			           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//			cerr<<"Compute the posterior of each directed path are now done."<<endl;


			cout << "\npath_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					cout << "path_posterior[" << j << "][" << i << "] = \t" << path_posterior[j][i] << endl;
				}
			}

			cout << "\n\nPlain version:" << endl;
			cout << "path_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					cout << path_posterior[j][i] << endl;
				}
			}

			ofstream of_path_poster("./results/path_poster");;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					of_path_poster << path_posterior[j][i] << endl;
				}
			}


			//get the estimated posterior of each limited length path (n*(n-1) paths) by weighting the number of dags in dag_list1
			double computeLimitedLengPathPosteriorTime;
			gettimeofday(&startTime, &myTimeZone);

			const int limited_leng = 2;

			double ** limited_leng_path_posterior;
			limited_leng_path_posterior= new double * [nh];
			for (int i = 0; i < nh; i++){
				limited_leng_path_posterior[i] = new double[nh];
			}

			get_limited_leng_path_poster_by_posterir_of_dags(limited_leng_path_posterior, nh, dag_list1, limited_leng);


			gettimeofday(&endTime, &myTimeZone);
			computeLimitedLengPathPosteriorTime = (endTime.tv_sec - startTime.tv_sec) +
			           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//			cerr<<"Compute the posterior of each directed limited length path are now done."<<endl;

			cout << "\nlimited_leng_path_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					cout << "limited_leng_path_posterior[" << j << "][" << i << "] = \t" << limited_leng_path_posterior[j][i] << endl;
				}
			}

			cout << "\n\nPlain version:" << endl;
			cout << "limited_leng_path_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					cout << limited_leng_path_posterior[j][i] << endl;
				}
			}

			ofstream of_limited_leng_path_posterior("./results/limited_leng_path_poster");;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					of_limited_leng_path_posterior << limited_leng_path_posterior[j][i] << endl;
				}
			}


			//get the estimated posterior of each combined path feature i ~~> j ~~> k (n*(n-1)*(n-2) paths) by weighting the posterior of dags in dag_list1
			double computeCombinedTwoPathPosteriorTime;
			gettimeofday(&startTime, &myTimeZone);

			double *** combined_two_path_posterior;

			combined_two_path_posterior = new double ** [nh];
			for (int i = 0; i < nh; i++){
				combined_two_path_posterior[i] = new double * [nh];
				for(int j = 0; j < nh; j++){
					combined_two_path_posterior[i][j] = new double [nh];
				}
			}

			get_combined_two_path_poster_by_posterir_of_dags(combined_two_path_posterior, nh, dag_list1);

			gettimeofday(&endTime, &myTimeZone);
			computeCombinedTwoPathPosteriorTime = (endTime.tv_sec - startTime.tv_sec) +
			           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//			cerr<<"Compute the posterior of each combination of two directed path are now done."<<endl;


			cout << "\ncombined_two_path_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					for(int k = 0; k < nh; k++){
						cout << "combined_two_path_posterior[" << k << "][" << j << "][" << i << "] = \t" << combined_two_path_posterior[k][j][i] << endl;
					}
				}
			}

			cout << "\n\nPlain version:" << endl;
			cout << "combined_two_path_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					for(int k = 0; k < nh; k++){
						cout << combined_two_path_posterior[k][j][i] << endl;
					}

				}
			}

			//combined_two_path_posterior[i][j][k] records the posterior of directed path i ~~> j ~~> k
			ofstream of_combined_two_path_poster("./results/combined_two_path_poster");;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					for(int k = 0; k < nh; k++){
						of_combined_two_path_poster << combined_two_path_posterior[k][j][i] << endl;
					}

				}
			}



			//get the estimated posterior of each combined path feature f4  (n*(n-1)*(n-2) paths) by weighting the posterior of dags in dag_list1
			//	  combined_two_f4_path_posterior[i][j][k] records the posterior of directed path i ~~> j and i ~~> k and j ~= k
			//	  the same as Matlab: C(i, j, k) = 1 iff R(i, j) = 1 and R(i, k) = 1 and j ~= k
			double computeCombinedTwo_f4_PathPosteriorTime;
			gettimeofday(&startTime, &myTimeZone);

			double *** combined_two_f4_path_posterior;

			combined_two_f4_path_posterior = new double ** [nh];
			for (int i = 0; i < nh; i++){
				combined_two_f4_path_posterior[i] = new double * [nh];
				for(int j = 0; j < nh; j++){
					combined_two_f4_path_posterior[i][j] = new double [nh];
				}
			}

			get_combined_two_f4_path_poster_by_posterir_of_dags(combined_two_f4_path_posterior, nh, dag_list1);

			gettimeofday(&endTime, &myTimeZone);
			computeCombinedTwo_f4_PathPosteriorTime = (endTime.tv_sec - startTime.tv_sec) +
			           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//			cerr<<"Compute the posterior of each f4 are now done."<<endl;


			cout << "\ncombined_two_f4_path_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					for(int k = 0; k < nh; k++){
						cout << "combined_two_f4_path_posterior[" << k << "][" << j << "][" << i << "] = \t" << combined_two_f4_path_posterior[k][j][i] << endl;
					}
				}
			}

			cout << "\n\nPlain version:" << endl;
			cout << "combined_two_f4_path_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					for(int k = 0; k < nh; k++){
						cout << combined_two_f4_path_posterior[k][j][i] << endl;
					}

				}
			}

			ofstream of_combined_two_f4_path_poster("./results/combined_two_f4_path_poster");;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					for(int k = 0; k < nh; k++){
						of_combined_two_f4_path_poster << combined_two_f4_path_posterior[k][j][i] << endl;
					}

				}
			}


			//get the estimated posterior of each combined path feature f5  (n*(n-1)*(n-2) paths) by weighting the posterior of dags in dag_list1
			//	  combined_two_f5_path_posterior[i][j][k] records the posterior of directed path i ~~> j and i not ~~> k and i ~= k
			//    	  the same as Matlab: C(i, j, k) = 1 iff R(i, j) = 1 and R(i, k) = 0 and i ~= k
			double computeCombinedTwo_f5_PathPosteriorTime;
			gettimeofday(&startTime, &myTimeZone);

			double *** combined_two_f5_path_posterior;

			combined_two_f5_path_posterior = new double ** [nh];
			for (int i = 0; i < nh; i++){
				combined_two_f5_path_posterior[i] = new double * [nh];
				for(int j = 0; j < nh; j++){
					combined_two_f5_path_posterior[i][j] = new double [nh];
				}
			}

			get_combined_two_f5_path_poster_by_posterir_of_dags(combined_two_f5_path_posterior, nh, dag_list1);

			gettimeofday(&endTime, &myTimeZone);
			computeCombinedTwo_f5_PathPosteriorTime = (endTime.tv_sec - startTime.tv_sec) +
			           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//			cerr<<"Compute the posterior of each f5 are now done."<<endl;


			cout << "\ncombined_two_f5_path_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					for(int k = 0; k < nh; k++){
						cout << "combined_two_f5_path_posterior[" << k << "][" << j << "][" << i << "] = \t" << combined_two_f5_path_posterior[k][j][i] << endl;
					}
				}
			}

			cout << "\n\nPlain version:" << endl;
			cout << "combined_two_f5_path_posterior:" << endl;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					for(int k = 0; k < nh; k++){
						cout << combined_two_f5_path_posterior[k][j][i] << endl;
					}

				}
			}


			ofstream of_combined_two_f5_path_poster("./results/combined_two_f5_path_poster");;
			for (int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					for(int k = 0; k < nh; k++){
						of_combined_two_f5_path_poster << combined_two_f5_path_posterior[k][j][i] << endl;
					}

				}
			}


			//write the estimated log posterior of each dag, i.e., log P(G|D), to the file
			ofstream of_dag_poster("./results/dag_poster");
			const int req_num_output_dags = 500;
			int num_output_dags = 0;
			for (list<dag>::const_iterator it1 = dag_list1.begin();
					it1 != dag_list1.end() && num_output_dags < req_num_output_dags;
					it1++){
				of_dag_poster << ((*it1).scores - sum_posterior_dag_list) << endl;
				num_output_dags++;
			}


			//free the memory
			//HR: add
			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					if(U_i_fam_vector_ptr[ch][fam] != NULL){
						(*U_i_fam_vector_ptr[ch][fam]).clear();
						delete U_i_fam_vector_ptr[ch][fam];
					}

				}
			}
			for(int i = 0; i < nh; i++){
				delete [] U_i_fam_vector_ptr[i];
			}
			delete [] U_i_fam_vector_ptr;


			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					if(U_i_fam_prob_cumul_vector_ptr[ch][fam] != NULL){
						delete [] U_i_fam_prob_cumul_vector_ptr[ch][fam];
					}
				}
			}
			for(int i = 0; i < nh; i++){
				delete [] U_i_fam_prob_cumul_vector_ptr[i];
			}
			delete [] U_i_fam_prob_cumul_vector_ptr;



			for(int i = 0; i < nh; i++){
				delete [] U_i_fam_vector_used_count[i];
			}
			delete [] U_i_fam_vector_used_count;



			for(int i = 0; i < nh; i++){
				delete [] edge_posterior[i];
			}
			delete [] edge_posterior;


			for(int i = 0; i < nh; i++){
				delete [] path_posterior[i];
			}
			delete [] path_posterior;


			for(int i = 0; i < nh; i++){
				delete [] limited_leng_path_posterior[i];
			}
			delete [] limited_leng_path_posterior;


			for(int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					delete [] combined_two_path_posterior[i][j];
				}
				delete [] combined_two_path_posterior[i];
			}
			delete [] combined_two_path_posterior;


			for(int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					delete [] combined_two_f4_path_posterior[i][j];
				}
				delete [] combined_two_f4_path_posterior[i];
			}
			delete [] combined_two_f4_path_posterior;



			for(int i = 0; i < nh; i++){
				for(int j = 0; j < nh; j++){
					delete [] combined_two_f5_path_posterior[i][j];
				}
				delete [] combined_two_f5_path_posterior[i];
			}
			delete [] combined_two_f5_path_posterior;




			cout << "\n\n\nTime Statistics" << endl;
			cout << "betaTime = " << betaTime <<  endl;
			cout << "alphaTime = " << alphaTime <<  endl;
			cout << "forwardTime = " << forwardTime <<  endl;


			cout << "sampleOrderTime = " << sampleOrderTime <<  endl;
			cout << "sampleDAGTime = " << sampleDAGTime <<  endl;
			cout << "computeEdgePosteriorTime = " << computeEdgePosteriorTime <<  endl;
			cout << "computePathPosteriorTime = " << computePathPosteriorTime <<  endl;
			cout << "computeLimitedLengPathPosteriorTime = " << computeLimitedLengPathPosteriorTime <<  endl;
			cout << "computeCombinedTwoPathPosteriorTime = " << computeCombinedTwoPathPosteriorTime <<  endl;
			cout << "computeCombinedTwo_f4_PathPosteriorTime = " << computeCombinedTwo_f4_PathPosteriorTime <<  endl;
			cout << "computeCombinedTwo_f5_PathPosteriorTime = " << computeCombinedTwo_f5_PathPosteriorTime <<  endl;



			double totalTime = betaTime + alphaTime + forwardTime
					+ sampleOrderTime + sampleDAGTime + computeEdgePosteriorTime;
			cout << "totalTime = " << totalTime <<  endl;

			string str1;
			str1.assign("./results/total_time_");

			str1.append((Arguments::datafile + 6));
			FILE* file_time = fopen(str1.c_str(), "a");
			fprintf(file_time,"%f \n", totalTime);
			fclose(file_time);


			double DPTime = betaTime + alphaTime + forwardTime;
			cout << "DPTime = " << DPTime <<  endl;

			string str1_2;
			str1_2.assign("./results/DP_time_");

			str1_2.append((Arguments::datafile + 6));
			FILE* file_time_2 = fopen(str1_2.c_str(), "a");
			fprintf(file_time_2,"%f \n", DPTime);
			fclose(file_time_2);


//			string str1_3;
//			str1_3.assign("./results/beta_time_");
//
//			str1_3.append((Arguments::datafile + 6));
//			FILE* file_time_3 = fopen(str1_3.c_str(), "a");
//			fprintf(file_time_3,"%f \n", betaTime);
//			fclose(file_time_3);


			string str1_4;
			str1_4.assign("./results/sampleOrder_time_");

			str1_4.append((Arguments::datafile + 6));
			FILE* file_time_4 = fopen(str1_4.c_str(), "a");
			fprintf(file_time_4,"%f \n", sampleOrderTime);
			fclose(file_time_4);


			string str1_5;
			str1_5.assign("./results/sampleDAG_time_");

			str1_5.append((Arguments::datafile + 6));
			FILE* file_time_5 = fopen(str1_5.c_str(), "a");
			fprintf(file_time_5,"%f \n", sampleDAGTime);
			fclose(file_time_5);



			string str2;
			str2.assign("./results/sum_prob_dags_");

			str2.append((Arguments::datafile + 6));
			FILE* file_sum_prob = fopen(str2.c_str(), "a");

			fprintf(file_sum_prob,"%18.6f \n", sum_posterior_dag_list);
			fclose(file_sum_prob);


//		}//if(Arguments::option == 11)



		//Finally
		//free the rest of the memory

		for (int j = 0; j < nh; j ++){
			delete [] beta[j];
			delete [] alpha[j];
		}
		delete [] beta;
		delete [] alpha;
		cerr<<"\n Posterior of each edge and posterior of each of five non-modular features now computed."<<endl;

	}//void compute_edge_probabilities_sample3(int h)





	//	  IW-DDS with huge N_o
	void compute_edge_probabilities_sample4(int h){
		cerr << "\n IW-DDS method with huge N_o:" << endl;

//		cerr<<" Compute edge probabilities for Layer "<<h<<":"<<endl;

		// Step 1: Compute beta[][];
		// Step 2: Compute alpha[][];
		int *Vh, *Vu; int nh, nu;
		model->layer(h, &Vh, &nh);
		model->upper_layers(h-1, &Vu, &nu);

//		totalOf2 << "nh=" << nh << endl;
//		for(int i = 0; i < nh; i++){
//			totalOf2 << "Vh[" << i << "] = " << Vh[i] << endl;
//		}
//		totalOf2 << "nu=" << nu << endl;
//		for(int i = 0; i < nu; i++){
//			totalOf2 << "Vu[" << i << "] = " << Vu[i] << endl;
//		}

		beta = new double*[nh];
		alpha = new double*[nh];

		k = model->max_indegree();

//		cerr<<" . "<<nh<<" nodes"<<endl;

		struct timeval startTime;
		struct timeval endTime;
		struct timezone myTimeZone;


		double betaTime;
		double alphaTime;
		double forwardTime;


		gettimeofday(&startTime, &myTimeZone);


		//If Arguments::ADtree == 1, user wants to use ADtree, then make the whole tree first
		if(Arguments::ADtree == 1){
			model->makeADTree();
		}

		for (int j = 0; j < nh; j ++){

			beta[j] = new double[1 << nh];
			//int i = Vh[j];
			compute_beta(j, Vh, nh, Vu, nu);

		}

		//If Arguments::ADtree == 1, then delete the whole tree
		if(Arguments::ADtree == 1){
			model->freeADTree();
		}

		gettimeofday(&endTime, &myTimeZone);
		betaTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//		cerr<<" . Tables beta are now ready."<<endl;

		gettimeofday(&startTime, &myTimeZone);
		for (int j = 0; j < nh; j ++){

			alpha[j] = new double[1 << nh];
			compute_alpha(j, Vh, nh);


		}
		gettimeofday(&endTime, &myTimeZone);
		alphaTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//		cerr<<" . Tables alpha are now ready."<<endl;


		// Step 3: Compute g_forward[];
		gf = new double[1 << nh];

		gettimeofday(&startTime, &myTimeZone);

		compute_g_forward(nh);

		gettimeofday(&endTime, &myTimeZone);
		forwardTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

		//HR: restore
		//cerr<<" . gf[all] = "<<gf[(1<<nh)-1]<<endl;
		cerr<<" \tDP step is done."<<endl;



		//HR: Step 4.1
		//HR: sample orders
		//HR: Reset the random number generator with the system clock.
		srand(time(NULL));
		//set the same value for the testing
		//srand(100);


		const int num_of_orders = Arguments::num_of_sampled_orders;

		int num_of_orders_per_run;

		if(nh <= 19){
			//num_of_orders_per_run = 20000;
			num_of_orders_per_run = 2000000;
		}
		else{
			//num_of_orders_per_run = 10000;
			num_of_orders_per_run = 1000000;
		}

		int num_of_runs;

		if(num_of_orders % num_of_orders_per_run == 0){
			num_of_runs = num_of_orders / num_of_orders_per_run;
		}
		else{
			num_of_runs = num_of_orders / num_of_orders_per_run + 1;
		}


		double sampleOrderTime = 0;
		double sampleDAGTime = 0;

		int numOfDAGsStoredInHD = 0;


		double * DAGScoresThresholdPtr = new double[num_of_runs - 1];


	    //free gf[] before DAG sampling step to save the memory
		//const double gf_all = gf[(1<<nh)-1];
		//delete [] gf;

		//each run
		for(int run_index = 0; run_index < num_of_runs; run_index++){
			cerr << "\n Run Index = " << run_index << endl;
			cout << "\n-------------------------" << endl;
			cout << "run_index = " << run_index << endl;


		gettimeofday(&startTime, &myTimeZone);


		list<total_order> t_order_list1;
		list<total_order> t_order_list2;
		total_order * t_order_temp_p;

		for(int order_ind = 0; order_ind < num_of_orders_per_run; order_ind++){

			t_order_temp_p = sample_one_order(nh);

			t_order_list1.push_back(* t_order_temp_p);
			t_order_list2.push_back(* t_order_temp_p);

			delete t_order_temp_p;
		}
//		cout << "t_order_list1.size() = " << t_order_list1.size() << endl;


		t_order_list2.sort();
//		cout << "*******************" << endl;
//		cout << "After sort():" << endl;
//		cout << "t_order_list2:" << endl;
//		cout << "t_order_list2.size() = " << t_order_list2.size() << endl;

		unique_HR(t_order_list2);
//		cout << "*******************" << endl;
//		cout << "See after uniqueHR(t_order_list2):" << endl;
//		cout << "t_order_list2:" << endl;
//		cout << "t_order_list2.size() = " << t_order_list2.size() << endl;

//		double sum_order_scores  = get_sum_order_scores(t_order_list2);
//		cout << "exp(sum_order_scores) = " << exp(sum_order_scores) << endl;

		t_order_list2.clear();


		gettimeofday(&endTime, &myTimeZone);
		sampleOrderTime += (endTime.tv_sec - startTime.tv_sec)
				+ (endTime.tv_usec - startTime.tv_usec) / 1000000.0;

		cerr << " \tOrder-Sampling per run is done." << endl;


		//HR: new Step 4.2
//		if(Arguments::option == 12){

			//double sampleDAGTime = 0;

			gettimeofday(&startTime, &myTimeZone);

//			if(nh > 17){
				t_order_list1.sort();
//				cout << "*******************" << endl;
//				cout << "After sort():" << endl;
//
//				cout << "t_order_list1.size() = " << t_order_list1.size() << endl;

//			}

			//Only remember parents (parent set) instead of the whole family to save the space cost
			vector<int> * * * U_i_fam_vector_ptr;
			U_i_fam_vector_ptr = new vector<int> * * [nh];
			for(int i = 0; i < nh; i++){
				U_i_fam_vector_ptr[i] = new vector<int> * [(1<<nh)];
			}

			//initialize U_i_raw_fam_vector_ptr
			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					U_i_fam_vector_ptr[ch][fam] = NULL;
				}
			}


			double * * * U_i_fam_prob_cumul_vector_ptr;
			U_i_fam_prob_cumul_vector_ptr = new double * * [nh];
			for(int i = 0; i < nh; i++){
				U_i_fam_prob_cumul_vector_ptr[i] = new double * [(1<<nh)];
			}

			//initialize U_i_fam_prob_cumul_vector_ptr
			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					U_i_fam_prob_cumul_vector_ptr[ch][fam] = NULL;
				}
			}
			//


			//the number of non-NULL U_i_fam_vector
			int U_i_fam_vector_setup_count = 0;

			//the number of cached families
			int family_setup_count = 0;


			//change int type into short type to save the memory
			short * * U_i_fam_vector_used_count;
			U_i_fam_vector_used_count = new short * [nh];
			for(int i = 0; i < nh; i++){
				U_i_fam_vector_used_count[i] = new short [(1<<nh)];
			}
			for (int ch = 0; ch < nh; ch++) {
				for (int fam = 0; fam < (1 << nh); fam++) {
					U_i_fam_vector_used_count[ch][fam] = 0;
				}
			}
			const int max_recorded_used_count = 10;


			list<dag> dag_list1;

			//DAGHashSet DAGHashSet1;

			int total_num_sampled_dags = 0;


			//HR: Use the on-fly version to avoid using too much memory
			//for each sampled total order <
			int order_num = 0;
			int numOfRecycle = 0;

			for (list<total_order>::iterator it1 = t_order_list1.begin(); it1 != t_order_list1.end(); it1++){
				order_num++;

				for(int ind = 0; ind < (*it1).length; ind++ ){
					int sigma_i = (*it1).torder[ind];
					int U_sigma_i = (*it1).torderBP[ind];
					if(U_i_fam_vector_ptr[sigma_i][U_sigma_i] == NULL){

						int numOfFam = getNumbOfFam(k, ind);

						U_i_fam_vector_ptr[sigma_i][U_sigma_i] = get_raw_family_vector(sigma_i, U_sigma_i, nh, k, numOfFam);
						assert(numOfFam == (int) (*U_i_fam_vector_ptr[sigma_i][U_sigma_i]).size());

						U_i_fam_prob_cumul_vector_ptr[sigma_i][U_sigma_i] = get_prob_cumul_array(U_i_fam_vector_ptr[sigma_i][U_sigma_i], alpha[sigma_i][U_sigma_i], sigma_i);

						U_i_fam_vector_setup_count++;

						if(U_i_fam_vector_used_count[sigma_i][U_sigma_i] < max_recorded_used_count){
							U_i_fam_vector_used_count[sigma_i][U_sigma_i]++;
						}


						family_setup_count += (*U_i_fam_vector_ptr[sigma_i][U_sigma_i]).size();


					}//if(U_i_fam_vector_ptr[sigma_i][U_sigma_i] == NULL)
					else {

						if(U_i_fam_vector_used_count[sigma_i][U_sigma_i] < max_recorded_used_count){
							U_i_fam_vector_used_count[sigma_i][U_sigma_i]++;
						}

					}

				}//for(int ind = 0; ind < (*it1).length; ind++ )


				int num_sampled_dags_per_order = 0;

				list<dag> order_dag_list1;

				dag samp_dag = sample_one_dag((*it1), U_i_fam_vector_ptr, U_i_fam_prob_cumul_vector_ptr);
				num_sampled_dags_per_order++;


				dag_list1.push_front(samp_dag);



				total_num_sampled_dags += num_sampled_dags_per_order;

				if(order_num % 20000 == 1){
//					cout << "\norder_num = " << order_num << endl;
//					cout << "total_num_sampled_dags = " << total_num_sampled_dags << endl;
//
//					cout << "dag_list1.size() = " << dag_list1.size() << endl;
//					cout << "U_i_fam_vector_setup_count = " << U_i_fam_vector_setup_count << endl;
//					cout << "family_setup_count = " << family_setup_count << endl;
				}


				//To reduce the number of DAGs in the dag_list1:
				//		sort in the descending order and then unique_HR
				if(order_num % 200000 == 1){
//					cout << "\norder_num = " << order_num << endl;
//					cout << "dag_list1.sort(descend_order_comparison) " << endl;
//					cout << "unique_HR(dag_list1) " << endl;

					dag_list1.sort(strict_descend_order_comparison);

					unique_HR(dag_list1);


				}


				const double Max_family_setup_num = 3.0e7; //for 2 G memory

				const double Min_released_family_setup_num = 1.9e7; //for 2 G memory

				//need to delete some cache if it occupies too much memory (when family_setup_count >= Max_family_setup_num)
				if(family_setup_count > Max_family_setup_num ){
//					cout << "\norder_num = " << order_num << endl;
					numOfRecycle++;

					recycleMemory(
							U_i_fam_vector_ptr,
							U_i_fam_prob_cumul_vector_ptr,
							U_i_fam_vector_setup_count,
							family_setup_count,
							U_i_fam_vector_used_count,
							Max_family_setup_num,
							Min_released_family_setup_num,
							nh);


				}//if(family_setup_count > Max_family_setup_num )


			}//for (vector<total_order>::iterator it1 = t_order_list1.begin(); ... )

//			cout << "\ntotal_num_sampled_dags = " << total_num_sampled_dags << endl;

			dag_list1.sort(strict_descend_order_comparison);

			unique_HR(dag_list1);

			dag_list1.sort(strict_descend_order_comparison);

//			cout << "dag_list1.size() = " << dag_list1.size() << endl;

			numOfDAGsStoredInHD += dag_list1.size();


			if(run_index == 0){
				int segmentLength = dag_list1.size() / num_of_runs;
//				cout << "run_index == 0" << endl;
//				cout << "segmentLength = " << segmentLength << endl;
				int numOfScannedDAGs = 0;
				for (list<dag>::const_iterator dag_list_it1 = dag_list1.begin(); dag_list_it1 != dag_list1.end(); dag_list_it1++){
					numOfScannedDAGs++;
					if((numOfScannedDAGs % segmentLength == 0) && (numOfScannedDAGs / segmentLength <= num_of_runs - 1)){
						DAGScoresThresholdPtr[numOfScannedDAGs / segmentLength - 1] = (*dag_list_it1).scores;
					}
				}

			}
			//



			//test write
			FILE* fout;

			string outStr;

			outStr.assign("dag_list_");
			char listNoStr1[40];
			//int v1 = 0;
			//sprintf(listNoStr1, "%d", v1);
			sprintf(listNoStr1, "%d", run_index);
			outStr.append(listNoStr1);

			fout = fopen(outStr.c_str(), "wb");
			write_DAGs_to_file(dag_list1, fout);
			fclose(fout);


			dag_list1.clear();
			t_order_list1.clear();

			//free the memory
			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					if(U_i_fam_vector_ptr[ch][fam] != NULL){
						(*U_i_fam_vector_ptr[ch][fam]).clear();
						delete U_i_fam_vector_ptr[ch][fam];
					}

				}
			}
			for(int i = 0; i < nh; i++){
				delete [] U_i_fam_vector_ptr[i];
			}
			delete [] U_i_fam_vector_ptr;


			for(int ch = 0; ch < nh; ch++){
				for(int fam = 0; fam < (1<<nh); fam++){
					if(U_i_fam_prob_cumul_vector_ptr[ch][fam] != NULL){
						delete [] U_i_fam_prob_cumul_vector_ptr[ch][fam];
					}
				}
			}
			for(int i = 0; i < nh; i++){
				delete [] U_i_fam_prob_cumul_vector_ptr[i];
			}
			delete [] U_i_fam_prob_cumul_vector_ptr;


			for(int i = 0; i < nh; i++){
				delete [] U_i_fam_vector_used_count[i];
			}
			delete [] U_i_fam_vector_used_count;



			gettimeofday(&endTime, &myTimeZone);
			sampleDAGTime += (endTime.tv_sec - startTime.tv_sec) +
					(endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

			cerr<<" \tDAG-Sampling per run is done."<<endl;

//			cout << "numOfRecycle = " << numOfRecycle << endl;



		}//for(int run_index = 0; run_index < num_of_runs; run_index++)


		for(int i = 0; i <= num_of_runs - 2; i++){
//			printf("DAGScoresThresholdPtr[%d] = %18.6f\n", i, DAGScoresThresholdPtr[i]);
		}


		//Now it is fine to free the following memory so that more DAGs can be stored in the memory
		delete [] gf;

		for (int j = 0; j < nh; j ++){
			delete [] beta[j];
			delete [] alpha[j];
		}
		delete [] beta;
		delete [] alpha;


		double computeEdgePosteriorTime;
		gettimeofday(&startTime, &myTimeZone);

		double ** edge_posterior;
		edge_posterior = new double *[nh];
		for (int i = 0; i < nh; i++) {
			edge_posterior[i] = new double[nh];
		}

		double sum_posterior_dag_list = LOG_ZERO;

//		cout << "\nnumOfDAGsStoredInHD = " << numOfDAGsStoredInHD << endl;


		computeEdgePosteriorInDirectNew(nh,
				num_of_runs,
				edge_posterior,
				sum_posterior_dag_list,
				DAGScoresThresholdPtr);

		gettimeofday(&endTime, &myTimeZone);
		computeEdgePosteriorTime = (endTime.tv_sec - startTime.tv_sec)
				+ (endTime.tv_usec - startTime.tv_usec) / 1000000.0;


//		printf("sum_posterior_dag_list = %18.6f\n", sum_posterior_dag_list);
//		cerr << "Compute the posterior of each edge are now done." << endl;

		cout << "\nedge_posterior:" << endl;
		for (int i = 0; i < nh; i++) {
			for (int j = 0; j < nh; j++) {
				cout << "edge_posterior[" << j << "][" << i << "] = \t"
						<< edge_posterior[j][i] << endl;
			}
		}

		cout << "\n\nPlain version:" << endl;
		cout << "edge_posterior:" << endl;
		for (int i = 0; i < nh; i++) {
			for (int j = 0; j < nh; j++) {
				cout << edge_posterior[j][i] << endl;
			}
		}

		//write the estimated posterior of each edge to the file
		ofstream of_edge_poster("./results/edge_poster");
		for (int i = 0; i < nh; i++) {
			for (int j = 0; j < nh; j++) {
				of_edge_poster << edge_posterior[j][i] << endl;
			}
		}
		//


			for(int i = 0; i < nh; i++){
				delete [] edge_posterior[i];
			}
			delete [] edge_posterior;


			delete [] DAGScoresThresholdPtr;


			cout << "\n\n\nTime Statistics" << endl;
			cout << "betaTime = " << betaTime <<  endl;
			cout << "alphaTime = " << alphaTime <<  endl;
			cout << "forwardTime = " << forwardTime <<  endl;

			cout << "sampleOrderTime = " << sampleOrderTime <<  endl;
			cout << "sampleDAGTime = " << sampleDAGTime <<  endl;
			cout << "computeEdgePosteriorTime = " << computeEdgePosteriorTime <<  endl;


			double totalTime = betaTime + alphaTime + forwardTime
					+ sampleOrderTime + sampleDAGTime + computeEdgePosteriorTime;
			cout << "totalTime = " << totalTime <<  endl;

			string str1;
			str1.assign("./results/total_time_");

			str1.append((Arguments::datafile + 6));
			FILE* file_time = fopen(str1.c_str(), "a");
			fprintf(file_time,"%f \n", totalTime);
			fclose(file_time);


			double DPTime = betaTime + alphaTime + forwardTime;
			cout << "DPTime = " << DPTime <<  endl;

			string str1_2;
			str1_2.assign("./results/DP_time_");

			str1_2.append((Arguments::datafile + 6));
			FILE* file_time_2 = fopen(str1_2.c_str(), "a");
			fprintf(file_time_2,"%f \n", DPTime);
			fclose(file_time_2);


//			string str1_3;
//			str1_3.assign("./results/beta_time_");
//
//			str1_3.append((Arguments::datafile + 6));
//			FILE* file_time_3 = fopen(str1_3.c_str(), "a");
//			fprintf(file_time_3,"%f \n", betaTime);
//			fclose(file_time_3);


			string str1_4;
			str1_4.assign("./results/sampleOrder_time_");

			str1_4.append((Arguments::datafile + 6));
			FILE* file_time_4 = fopen(str1_4.c_str(), "a");
			fprintf(file_time_4,"%f \n", sampleOrderTime);
			fclose(file_time_4);


			string str1_5;
			str1_5.assign("./results/sampleDAG_time_");

			str1_5.append((Arguments::datafile + 6));
			FILE* file_time_5 = fopen(str1_5.c_str(), "a");
			fprintf(file_time_5,"%f \n", sampleDAGTime);
			fclose(file_time_5);



			string str2;
			str2.assign("./results/sum_prob_dags_");

			str2.append((Arguments::datafile + 6));
			FILE* file_sum_prob = fopen(str2.c_str(), "a");
			//fprintf(file_sum_prob,"%f \n", exp(sum_posterior_dag_list));
			fprintf(file_sum_prob,"%18.6f \n", sum_posterior_dag_list);
			fclose(file_sum_prob);


//		}//if(Arguments::option == 12)



		cerr<<" Posterior of each edge now computed."<<endl;

	}//void compute_edge_probabilities_sample4(int h)






	//HR: new mix with max in-degree method,
	// Computes probabilities for edges pointing to the h-th layer:
	//HR: to vertex uu, from vertex vv; uu, vv \in {0, 1, ..., nh-1}
	void compute_edge_probabilities_mixIndegree(int h){
		cerr<<" Mix with max in-degree method: " << endl;
//		cerr<<" Compute edge probabilities for Layer "<<h<<":"<<endl;

		// Step 1: Compute beta[][];
		// Step 2: Compute alpha[][];
		int *Vh, *Vu; int nh, nu;
		model->layer(h, &Vh, &nh);
		model->upper_layers(h-1, &Vu, &nu);


//		totalOf3 << "nh=" << nh << endl;
//		for(int i = 0; i < nh; i++){
//			totalOf3 << "Vh[" << i << "] = " << Vh[i] << endl;
//		}
//		totalOf3 << "nu=" << nu << endl;
//		for(int i = 0; i < nu; i++){
//			totalOf3 << "Vu[" << i << "] = " << Vu[i] << endl;
//		}


		//Step 1(a) (b)

		beta = new double*[nh];
		alpha = new double*[nh];

		k = model->max_indegree();

//		cerr<<" . "<<nh<<" nodes"<<endl;

		struct timeval startTime;
		struct timeval endTime;
		struct timezone myTimeZone;


		double betaTime;
		double alphaTime;
		double RRFuncTime;
		double HFuncTime;
		double KFuncTime;
		double GammaTime;

		double step2Time;


		//Step 1(a)

		gettimeofday(&startTime, &myTimeZone);


		//If Arguments::ADtree == 1, user wants to use ADtree, then make the whole tree first
		if(Arguments::ADtree == 1){
			//cout << "Arguments::ADtree == 1"  << endl;
			model->makeADTree();
		}
		else{
			//cout << "Arguments::ADtree == 0"  << endl;
		}

		for (int j = 0; j < nh; j ++){

			beta[j] = new double[1 << nh];

			//int i = Vh[j];
			compute_beta(j, Vh, nh, Vu, nu);

		}



		//If Arguments::ADtree == 1, then delete the whole tree
		if(Arguments::ADtree == 1){
			model->freeADTree();
		}

		gettimeofday(&endTime, &myTimeZone);
		betaTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;



		cerr<<"  Tables beta have been computed." << endl;



		//step 1(b)
		gettimeofday(&startTime, &myTimeZone);
		for (int j = 0; j < nh; j ++){

			alpha[j] = new double[1 << nh];
			compute_alpha(j, Vh, nh);


		}
		gettimeofday(&endTime, &myTimeZone);
		alphaTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;


		cerr<<"  Tables alpha have been computed." << endl;



//		totalOfMix << "\nalpha" << endl;
//		for (int j = 0; j < nh; j ++){
//			for(int index=0; index < (1 << nh); index++){
//				totalOfMix<< "alpha[" << j << "][" << index << "] =" <<  alpha[j][index] << endl;
//			}
//			totalOfMix << endl;
//		}
//		totalOfMix << endl;


		//Step 1(c)
		//HR: backward method to compute: for all S \subset V RR(S)
		RRFunc = new double[1 << nh];

//		cerr<<" . To compute Table RR" << endl;

		gettimeofday(&startTime, &myTimeZone);

		compute_RR_Fast(nh, RRFunc, alpha);

		gettimeofday(&endTime, &myTimeZone);
		RRFuncTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;

//		cerr<<"  P(D) = RRFunc[all] = "<<RRFunc[(1<<nh)-1]<<endl;
//		fprintf(stderr, "P(D) = RRFunc[all] = %18.8f\n", RRFunc[(1<<nh)-1]);
		cerr<<"  Table RR has been computed." <<endl;



//		totalOfMix << "\nRRFunc:" << endl;
//		for(int index = 0; index < (1<<nh) ; index++){
//			totalOfMix << "RRFunc[" << index << "] = " << RRFunc[index] << endl;
//		}


		//Step 1(d)
		//HR: forward method to compute: for all S \subset V H(S)

		hFunc = new double[1 << nh];

//		cerr<<" . To compute Table H" << endl;

		gettimeofday(&startTime, &myTimeZone);

		compute_h_Fast(nh, hFunc, alpha);

		gettimeofday(&endTime, &myTimeZone);
		HFuncTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;


//		cerr<<"hFunc[all] = "<<hFunc[(1<<nh)-1]<<endl;
//		fprintf(stderr, "hFunc[all] = %18.8f\n", hFunc[(1<<nh)-1]);
		cerr<<"  Table H has been computed." <<endl;

		//Added to record p(D)
//		cout << "\nRRFunc[all] = " << RRFunc[(1<<nh)-1] << endl;
		fprintf(stdout, "p(D) = RRFunc[all] = %18.8f\n\n", RRFunc[(1<<nh)-1]);
//		cout << "hFunc[all] = " << hFunc[(1<<nh)-1] << endl;
//		fprintf(stdout, "p(D) = hFunc[all] = %18.8f\n\n", hFunc[(1<<nh)-1]);


//		totalOfMix << "\nhFunc:" << endl;
//		for(int index = 0; index < (1<<nh) ; index++){
//			totalOfMix << "hFunc[" << index << "] = " << hFunc[index] << endl;
//		}

		//Step 1 (e)
		KFunc = new double*[nh];
//		cerr<<" . To compute Table K" << endl;

		gettimeofday(&startTime, &myTimeZone);
		for(int i = 0; i < nh; i++){
			//cerr<<" For v = " << i << endl;
			KFunc[i] = new double [1<<nh];
			compute_all_K_v(nh, KFunc[i], i, alpha);

		}

		gettimeofday(&endTime, &myTimeZone);
		KFuncTime = (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;


		cerr<<"  Table K has been computed." <<endl;

//		totalOfMix << "\nKFunc:" << endl;
//		for(int i = 0; i < nh; i++){
//			for(int index = 0; index < (1<<nh) ; index++){
//				totalOfMix <<  "KFunc[" << i <<"][" << index << "] = "  << KFunc[i][index] << endl;
//			}
//			totalOfMix << endl;
//		}




		// Step 1 (f): For all i \in V, for all Pa_i \subset V - {i} with |Pa_i| <= k, compute Gamma_i(Pa_i)

		double * Gamma = new double[1 << nh];
		double * Gfb = new double[1 << nh];
//		ofstream ofLog("debugEngineLog.txt");



		double ppDirectedFromTo[nh][nh];
		double ppUndirected[nh][nh];
		for(int i = 0; i < nh; i++){
			for(int j = 0; j < nh; j++){
				ppDirectedFromTo[i][j] = 0;
				ppUndirected[i][j] = 0;

			}
		}



		GammaTime = 0.0;
		step2Time = 0.0;

		for (int v = 0; v < nh; v++){
//			cerr<<" For v = " << v <<": " << endl;

			//start counting the GammaTime for each v
			gettimeofday(&startTime, &myTimeZone);

			// Compute Gfb[].
			for (int U = 0; U < (1<<nh); U++){
				//if v \in the binary repre of U
				if ((U>>v) & 1){
					Gfb[U] = LOG_ZERO;
				}
				else{

					if(KFunc[v][U] == LOG_ZERO){
						Gfb[U] = LOG_ZERO;
					}

					else if(KFunc[v][U] <= 0){
						Gfb[U] = hFunc[U] + KFunc[v][U];
					}
					else{

//						ofLog << "KFunc[v][U] > 0" << endl;
//						ofLog << "KFunc[" << v << "][" << U << "] = " << KFunc[v][U] << endl;
//						ofLog << "Gfb > 0" << endl;
//
						Gfb[U] = -(hFunc[U] - KFunc[v][U]);
//
//						ofLog << "Gfb[" << U << " ] = " << Gfb[U] << endl;
					}


					//cerr<<" j = "<<j<<": ";
					//print_set(cerr, S);
					//cerr<<"; ";
					//print_set(cerr, cS);
					//cerr<<": "<<gfb[S];
					//cerr<<endl;
				}
			}


			//cerr<< "Table Gfb has been computed" << endl;



//			ofLog <<"For v = " << v <<": " << endl;
//			ofLog << "\nGfb: " << endl;
//			for(int index = 0; index < (1 << nh); index++){
//				ofLog << "Gfb[" << index << " ] =" << Gfb[index] << endl;
//			}
//			ofLog << endl;

			// Compute gamma_j() := gamma[].


			fdmtLogAddComp(Gfb, Gamma, nh, k);

			gettimeofday(&endTime, &myTimeZone);
			 //end counting GammaTime for each v
			GammaTime += (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;





//			ofLog <<"For v = " << v <<": " << endl;
//			ofLog << "\nGamma: " << endl;
//			for(int index = 0; index < (1 << nh); index++){
//				ofLog << "Gamma[" << index << " ] =" << Gamma[index] << endl;
//			}
//			ofLog << endl;

//			cerr<< "   Table Gamma has been computed." << endl;


			//cerr<< " Gfb[0]: " <<Gfb[0]<<endl;
			//cerr<< " Gamma[0]: " << Gamma[0]<<endl;


			//Step 2
			//HR: v is the to vertex, i is the from vertex
//			cerr<<" . Incoming edges for node "<< v <<":";
			// Loop over start-point nodes i.


			//start counting step2Time for each v
			gettimeofday(&startTime, &myTimeZone);

			for (int i = 0; i < nh; i++){
				if (i == v) continue;
//				cerr<<" "<< i;
				double log_prob_nume = eval_edge_mixIndegree(i,  v, Gamma, beta[v], nh, k);

				model->print_edge_prob(
					cout, Vh[i], Vh[v], exp(log_prob_nume - RRFunc[(1<<nh)-1])); //HR: - gb[(1<<nh)-1] mean /L(v)

				//model->print_edge_prob(
					//ofResult, Vh[i], Vh[v], exp(log_prob_nume - RRFunc[(1<<nh)-1])); //HR: - gb[(1<<nh)-1] mean /L(v)

				ppDirectedFromTo[i][v] = exp(log_prob_nume - RRFunc[(1<<nh)-1]);

			}


			gettimeofday(&endTime, &myTimeZone);
			//end counting step2Time for each v
			step2Time += (endTime.tv_sec - startTime.tv_sec) +
		           (endTime.tv_usec - startTime.tv_usec) / 1000000.0 ;


//			cerr<<endl;

		} //for (int v = 0; v < nh; v++)


		cout << "\n\n\nTime Statistics" << endl;
		cout << "betaTime = " << betaTime <<  endl;
		cout << "alphaTime = " << alphaTime <<  endl;
		cout << "RRFuncTime = " << RRFuncTime <<  endl;
		cout << "HFuncTime = " << HFuncTime <<  endl;
		cout << "KFuncTime = " << KFuncTime <<  endl;

		cout << "GammaTime = " << GammaTime <<  endl;
		cout << "step2Time = " << step2Time <<  endl;
		double totalTime = betaTime + alphaTime + RRFuncTime
						+ HFuncTime + KFuncTime + GammaTime + step2Time;
		cout << "totalTime = " << totalTime <<  endl;





		//HR:
		//note the order of the double loop, consistent with the column order of Matlab
		cout << "\nExact_ppDirectedFromTo \t Exact_ppDirectedFromTo" << endl;
		for(int i = 0; i < nh; i++){
			for(int j = 0; j < nh; j++){
				cout << "ppDirectedFromTo[" << j << "][" << i << "] \t " << ppDirectedFromTo[j][i] << endl;


			}
		}


		cout << "\n\nPlain output:" << endl;
		cout << "Exact_ppDirectedFromTo" << endl;
		for(int i = 0; i < nh; i++){
			for(int j = 0; j < nh; j++){
				cout << ppDirectedFromTo[j][i] << endl;

			}
		}

//		cout << "\nExact_ppUndirected \t Exact_ppUndirected" << endl;
//		for(int i = 0; i < nh; i++){
//			for(int j = i + 1; j < nh; j++){
//				ppUndirected[i][j] = ppDirectedFromTo[i][j] + ppDirectedFromTo[j][i];
//				cout << "ppUndirected[" << i << "][" << j << "] \t " << ppUndirected[i][j] << endl;
//
//
//			}
//		}


		delete [] Gamma;
		delete [] Gfb;

		delete [] hFunc;

		delete [] RRFunc;

		delete [] Eta_U;

		for (int j = 0; j < nh; j ++){
			delete [] beta[j];
			delete [] alpha[j];

			delete [] KFunc[j];

		}
		delete [] beta;
		delete [] alpha;

		delete [] KFunc;

		//cerr<<" Edge probabilities now computed. (test)"<<endl; //
		cerr<<" Posterior of each edge now computed."<<endl;

		return;



	}
	//HR: end mixIndegree




	void compute_beta(int j, int *Vh, int nh, int *Vu, int nu){

		//cerr<<"Compute beta, j = "<<j<<", i = "<< Vh[j];
		//cerr<<", k = "<<k<<", nh = "<<nh<<":"<<endl;

//		totalOf2 <<"Compute beta, j = "<< j << endl;
		//totalOf2 << "call: sub_init()" << endl;


		//sub_init(0, 0, beta[j], MARK, 0, nh);
		init_beta_new(beta[j], LOG_ZERO, nh);

		//HR:
		//totalOf2 << "return: sub_init()" << endl;


		//HR: restore
		//cerr<<" Initialization done."<<endl;

		int *S = new int[k];

		//totalOf2 << "call: sub_beta()" << endl;
		sub_beta(-1, nh, beta[j], 0, S, Vh[j], Vu, nu, 0);
		//totalOf2 << "return: sub_beta()" << endl;
		delete [] S;


	}// end void compute_beta()





	void compute_betaNume(int j, int *Vh, int nh, int *Vu, int nu, int uu, int vv){

		//cerr<<"Compute betaNume, j = "<<j<<", i = "<< Vh[j];
		//cerr<<", k = "<<k<<", nh = "<<nh<<":"<<endl;
		//cerr<<"uu = "<< uu << ", vv = "<< vv <<endl;

		sub_init(0, 0, betaNume[j], MARK, 0, nh);

		//HR: restore
		//cerr<<" Initialization done."<<endl;

		int *S = new int[k];

		//check whether uu == j
		if(uu != j){
			sub_beta(-1, nh, betaNume[j], 0, S, Vh[j], Vu, nu, 0);
		}
		else{
			//cerr << "uu == j == " << uu << ": so call sub_betaNume()" << endl;
			//totalOf3 << "uu == j == " << uu << ": so call sub_betaNume()" << endl;
			sub_betaNume(-1, nh, betaNume[j], 0, S, Vh[j], Vu, nu, 0, uu, vv);
		}



		delete [] S;
	}




	void compute_betaInOut(int j, int *Vh, int nh, int *Vu, int nu){

		//cerr << "call: sub_init(0, 0, betaInOut[j], MARK, 0, nh)" << endl;
		sub_init(0, 0, betaInOut[j], MARK, 0, nh);
		//cerr << "return: sub_init(0, 0, betaInOut[j], MARK, 0, nh)" << endl;

		//HR: restore
		//cerr<<" Initialization done."<<endl;

		int *S = new int[k];

		//cerr << "call: sub_betaInOut()" << endl;

		sub_betaInOut(-1, nh, betaInOut[j], 0, S, Vh[j], Vu, nu, 0);
		//sub_beta(-1, nh, betaInOut[j], 0, S, Vh[j], Vu, nu, 0);

		//cerr << "return: sub_betaInOut()" << endl;

		delete [] S;

	}//end of compute_betaInOut()





	void sub_beta(int jprev, int nh, double *b, int d, int *S, int i, int *Vu, int nu, int T){

		//cerr<<" S:"; print_nodes(cerr, S, d);
		//cerr<<"; T:"; print_nodes(cerr, T, Vu, nu); cerr<<endl;


		double lb;
		//if(Arguments::option < 1){
		if(Arguments::option < 1){
			//rebel method

			if(Arguments::ADtree == 0){

				if(Arguments::rebel_para_choice == 1){

					//Original rebel
					//K2 score with prior  1 / nchoosek(n-1, |Pa_i|)  for \rho_i(Pa_i)
					lb = model->log_prior(i, S, d) + model->log_lcp(i, S, d);
				}
				else if(Arguments::rebel_para_choice == 2){
					//rebelOld
					//BDeu score with prior  1 / nchoosek(n-1, |Pa_i|)  for \rho_i(Pa_i)
					lb = model->log_prior(i, S, d) + model->log_lcpHR(i, S, d);
				}
				else if(Arguments::rebel_para_choice == 3){
					//rebelNew
					//BDeu score with prior  1  for \rho_i(Pa_i)
					lb = 0 + model->log_lcpHR(i, S, d); 	//rebleNew
				}

			}
			else{

				if(Arguments::rebel_para_choice == 1){
					//Original rebel
					lb = model->log_prior(i, S, d) + model->log_lcp(i, S, d);
				}
				else if(Arguments::rebel_para_choice == 2){
					//rebelOld
					lb = model->log_prior(i, S, d) + model->log_lcpHR(i, S, d);
				}
				else if(Arguments::rebel_para_choice == 3){
					//rebelNew
					lb = 0 + model->log_lcpHR(i, S, d); 	//rebleNew
				}
			}

		}
		else if(Arguments::option == 8 || Arguments::option == 9){

			if(Arguments::ADtree == 0){
				//K2 score with prior  1 / nchoosek(n-1, |Pa_i|)  for \rho_i(Pa_i)
				lb = model->log_prior(i, S, d) + model->log_lcp(i, S, d);

			}
			else{
				//K2 score with prior  1 / nchoosek(n-1, |Pa_i|)  for \rho_i(Pa_i)
				lb = model->log_prior(i, S, d) + model->log_lcp(i, S, d);
			}

		}
		//All the other cases except opt 0, 8, 9
		else if(Arguments::option >= 1) {

			//BDeu score with prior  1  for \rho_i(Pa_i)
			if(Arguments::ADtree == 0){
				lb = 0 + model->log_lcpHR(i, S, d);
			}
			else{
				lb = 0 + model->log_lcpHR_ADtree(i, S, d);
			}
			//totalOf2 << "end: lb = 0 + model->log_lcpHR(i, S, d)" << endl;
		}


		logAdd_New(b[T], lb);


		if (d < k){
			for (int j = jprev + 1; j < nu; j ++){

				if (Vu[j] != i){

					S[d] = Vu[j];

					if (j < nh) {

						sub_beta(j, nh, b, d+1, S, i, Vu, nu, T | (1 << j));
					}

					else{

						sub_beta(j, nh, b, d+1, S, i, Vu, nu, T);
					}
				}
			}
		}
		//totalOf2 <<" End: sub_beta()" << endl;
	}



	void sub_betaNume(int jprev, int nh, double *b, int d, int *S, int i, int *Vu, int nu, int T, int uu, int vv){

		//cerr<<" S:"; print_nodes(cerr, S, d);
		//cerr<<"; T:"; print_nodes(cerr, T, Vu, nu); cerr<<endl;

		int vvInS = 0;


		for(int ind = 0; ind < d; ind++){
			if(vv == S[ind]){
				vvInS = 1;
				break;
			}
		}

		if(vvInS == 1){
//			totalOf3 << "vvInS == 1: double lb = model->log_prior(i, S, d) + model->log_lcp(i, S, d);" << endl;

			double lb;
			if(Arguments::option < 1){

				if(Arguments::ADtree == 0){
					lb = model->log_prior(i, S, d) + model->log_lcpHR(i, S, d);
				}
				else{
					lb = model->log_prior(i, S, d) + model->log_lcpHR_ADtree(i, S, d);
				}
			}
			else if(Arguments::option >= 1) {
				//HR: for 1 prior for q(G_i) and new Dirichlet hyper-param

				if(Arguments::ADtree == 0){
					lb = 0 + model->log_lcpHR(i, S, d);
				}
				else{
					lb = 0 + model->log_lcpHR_ADtree(i, S, d);
				}

			}


			logAdd_New(b[T], lb);
		}
		else{

			double tempLogB = LOG_ZERO;
			logAdd_New(b[T], tempLogB);
		}



		if (d < k){
			for (int j = jprev + 1; j < nu; j ++){

				if (Vu[j] != i){

					S[d] = Vu[j];

					if (j < nh) {

						sub_betaNume(j, nh, b, d+1, S, i, Vu, nu, T | (1 << j), uu, vv);
					}

					else{

						sub_betaNume(j, nh, b, d+1, S, i, Vu, nu, T, uu, vv);
					}
				}
			}
		}
	}





	void getInOutFeatures(int no_vars){

		for(int i = 0; i < no_vars; i++){
			inFeatureEdges[i] = 0;
		}
		for(int i = 0; i < no_vars; i++){
			outFeatureEdges[i] = 0;
		}

		//string str;
		//str.assign(Arguments::inFeasFileName);
		//FILE * fpIn = fopen(str.c_str(), "r");

		FILE * fpIn = fopen(Arguments::inFeasFileName, "r");

		int fromNode;
		int toNode;
		//char c;
		//int i=0;
		while(!feof(fpIn)){

			fromNode = -1;
			toNode = -1;


			fscanf(fpIn,"%d->%d\n",&fromNode, &toNode);

			//cout << "fromNode = "<< fromNode << endl;
			//cout << "toNode = "<< toNode << endl;

			//To guard against toNode == -1
			if(toNode == -1){
				break;
			}
			else{
				inFeatureEdges[toNode] += (1 << fromNode);
			}

		}

		fclose(fpIn);

//		for(int i = 0; i < no_vars; i++){
//			cout << "inFeatureEdges[" << i << "] " << inFeatureEdges[i] << endl;
//		}


		FILE *fp2 = fopen(Arguments::outFeasFileName, "r");

		while(!feof(fp2)){

			fromNode = -1;
			toNode = -1;

			fscanf(fp2,"%d->%d\n",&fromNode, &toNode);

			//cout << "fromNode = "<< fromNode << endl;
			//cout << "toNode = "<< toNode << endl;

			if(toNode == -1){
				break;
			}
			else{
				outFeatureEdges[toNode] += (1 << fromNode);
			}


		}
		fclose(fp2);

//		for(int i = 0; i < no_vars; i++){
//			cout << "outFeatureEdges[" << i << "] " << outFeatureEdges[i] << endl;
//		}

	} //end of getInOutFeatures




	void sub_betaInOut(int jprev, int nh, double *b, int d, int *S, int i, int *Vu, int nu, int T){

//		cerr<< "Start sub_betaInOut()" << endl;
//		cerr<<" S:"; print_nodes(cerr, S, d);
//		cerr<<" T: " << T << endl;
//		cerr<<" d: " << d << endl;
		//cerr<<"; T:"; print_nodes(cerr, T, Vu, nu); cerr<<endl;



		//HR: add the judge for in and out features

		int compBetaNeeded = 0;

		//cerr << "inFeatureEdges[" << i << "] = " << inFeatureEdges[i] << endl;
		//cerr << "outFeatureEdges[" << i << "] = " << outFeatureEdges[i] << endl;

		//if (each of inFeatureEdges[i] is in T ) && (none of outFeatureEdges[i] is in T )
		if(((inFeatureEdges[i] & T) == inFeatureEdges[i]) && ((outFeatureEdges[i] & T) == 0) ){
			compBetaNeeded = 1;
			//cerr << "inFeatureEdges[" << i << "] = " << inFeatureEdges[i] << endl;
			//cerr << "outFeatureEdges[" << i << "] = " << outFeatureEdges[i] << endl;
			//cerr << "T = " << T << endl;
			//cerr<<" compBetaNeeded = 1;" << endl;

		}

		if(compBetaNeeded == 1){
			double lb;
			//Arguments::option must be 5
			if(Arguments::option < 1){

				if(Arguments::ADtree == 0){
					lb = model->log_prior(i, S, d) + model->log_lcpHR(i, S, d);
				}
				else{
					lb = model->log_prior(i, S, d) + model->log_lcpHR_ADtree(i, S, d);
				}
			}
			else if(Arguments::option >= 1) {
				//HR: for 1 prior for q(G_i) and new Dirichlet hyper-param

				if(Arguments::ADtree == 0){
					lb = 0 + model->log_lcpHR(i, S, d);
				}
				else{
					lb = 0 + model->log_lcpHR_ADtree(i, S, d);
				}
				//cerr << "end: double lb = 0 + model->log_lcpHR(i, S, d)" << endl;
			}

			//double lb = model->log_prior(i, S, d) + model->log_lcp(i, S, d);
			//cerr << "start: double lb = 0 + model->log_lcpHR(i, S, d)" << endl;
			//double lb = 0 + model->log_lcpHR(i, S, d);
			//cerr << "end: double lb = 0 + model->log_lcpHR(i, S, d)" << endl;

			logAdd_New(b[T], lb);
		}
		else{

			double tempLogB = LOG_ZERO;
			logAdd_New(b[T], tempLogB);
		}

		if (d < k){
			for (int j = jprev + 1; j < nu; j ++){

				if (Vu[j] != i){

					S[d] = Vu[j];

					if (j < nh) {

						sub_betaInOut(j, nh, b, d+1, S, i, Vu, nu, T | (1 << j));
					}

					else{

						sub_betaInOut(j, nh, b, d+1, S, i, Vu, nu, T);
					}
				}
			}
		}
	} //end of sub_betaInOut




	void sub_init(int d, int S, double *a, double value, int ones, int nh){
		//HR: restore
		//cerr<<" d: "<<d<<endl;



		a[S] = value;


		if (d < nh){
			sub_init(d+1, S, a, value, ones, nh);

			if (ones < k) sub_init(d+1, S | (1 << d), a, value, ones+1, nh);
		}
	}//void sub_init(int d, int S, double *a, double value, int ones, int nh)



	void init_beta_new(double * b, double value, int nh){
		for(int i = 0; i < (1 << nh); i++){
			b[i] = value;
		}
	}//void init_beta_new(double * b, double value, int nh)




	total_order * sample_one_order(int nh){
		const int n = nh;

		//condSet stores {-1, -1, ..., -1, sigma_{k+1}, sigma_{k+2}, sigma_n}
		int condSet[n];
		for(int i = 0; i < n; i++){
			condSet[i] = -1;
		}

		//condSetBP is the binary representation for (sigma_{k+1}, sigma_{k+2}, sigma_n)
		int condSetBP = 0;


		//U_sigma_k_add_1 is the binary representation for U_sigma_{k+1}
		int U_sigma_k_add_1 = ((1 << n) - 1) - condSetBP;

		//sigma_k_val_vector stores all the valid choices for sigma_k, with length k
		vector<int> sigma_k_val_vector;

		//sigma_k_valBP_vector stores all the valid choices for the binary representation of (sigma_k, sigma_{k+1}, sigma_{k+2}, sigma_n), with length k
		vector<int> sigma_k_valBP_vector;

		//sigma_k_prob_vector stores all the log prob. for p(sigma_k | (sigma_{k+1}, sigma_{k+2}, sigma_n), D), with length k
		vector<double> sigma_k_prob_vector;


		total_order * t_order_ptr1 = new total_order(n);


		for(int sampledOrderLength = 0; sampledOrderLength < n; sampledOrderLength++){
			sigma_k_val_vector.clear();
			sigma_k_valBP_vector.clear();
			sigma_k_prob_vector.clear();

			for(int i = 0; i < n; i++){
				//if node i is in U_sigma_{k+1}
				if(((1 << i) & U_sigma_k_add_1) != 0){
					sigma_k_val_vector.push_back(i);

					//the binary representation of (sigma_k, sigma_{k+1}, sigma_{k+2}, sigma_n)
					int sigma_k_valBP = condSetBP + (1 << i);
					sigma_k_valBP_vector.push_back(sigma_k_valBP);

					int U_sigma_k = ((1 << n) - 1) - sigma_k_valBP;
					//double sigma_k_prob = exp( gf[U_sigma_k] + alpha[i][U_sigma_k] - gf[U_sigma_k_add_1]  );
					double sigma_k_prob = gf[U_sigma_k] + alpha[i][U_sigma_k] - gf[U_sigma_k_add_1] ;
					sigma_k_prob_vector.push_back(sigma_k_prob);

				}
			}
			assert((int) sigma_k_val_vector.size() == n - sampledOrderLength);


			int chosen_index;

			chosen_index = sample_from_dist(sigma_k_prob_vector);

			condSet[ n - 1 - sampledOrderLength] = sigma_k_val_vector[chosen_index];
			condSetBP = sigma_k_valBP_vector[chosen_index];
			U_sigma_k_add_1 = ((1 << n) - 1) - condSetBP;

			(* t_order_ptr1).torder[n - 1 - sampledOrderLength] = condSet[ n - 1 - sampledOrderLength];
			(* t_order_ptr1).torderBP[n - 1 - sampledOrderLength] = U_sigma_k_add_1;
			(* t_order_ptr1).scores += sigma_k_prob_vector[chosen_index];

		}
//		(* t_order_ptr1).print();

		//check the match of scores
		double real_scores = 0;
		for(int i = 0; i < n; i++){
			real_scores += alpha[(* t_order_ptr1).torder[i]][(* t_order_ptr1).torderBP[i]];
		}
		real_scores -= gf[(1 << n) - 1];
//		cout << "real_scores = " << real_scores << endl << endl;
		assert(fabs(real_scores - (* t_order_ptr1).scores) < 0.00001);

		//delete t_order_ptr1;
		return t_order_ptr1;

	}//total_order * sample_one_order(int nh)


	//HR: Generate a random number between 0 and 1
	double HR_unif_rand(){
	    return rand() / ( (double) RAND_MAX + 1 );
	}




	int sample_from_dist(const vector<double> & log_prob_vector){

		double rand_prob = HR_unif_rand();
		int chosen_index = -1;
		double cumu_prob = 0;
		for(unsigned int i = 0; i < log_prob_vector.size(); i++){
			cumu_prob += exp(log_prob_vector[i]);
			if(rand_prob <= cumu_prob){
				chosen_index = i;
				break;
			}
		}
		return chosen_index;
	}//int sample_from_dist(const vector<double> & log_prob_vector)



	void unique_HR(list<total_order> & t_order_list1){
		list<total_order>::iterator it1;
		list<total_order>::iterator it2;
		for (it1 = t_order_list1.begin(); it1 != t_order_list1.end(); it1++){
			it2 = it1;
			it2++;
			//while two orders have the same scores
			while((*it2) == (*it1)){
				if(is_same_order((*it2), (*it1))){
					it2 = t_order_list1.erase(it2);
				}
				else{
					it2++;
				}

			}

		}
	}//void unique_HR(list<total_order> & t_order_list1)





	void unique_HR(list<dag> & dag_list1){
		list<dag>::iterator it1;
		list<dag>::iterator it2;
		for (it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
			it2 = it1;
			it2++;

			while(haveClosescores((*it2), (*it1))){

				if(is_same_dag((*it2), (*it1))){

					it2 = dag_list1.erase(it2);
				}
				else{
					it2++;
				}

			}

		}
	}//void unique_HR(list<dag> & dag_list1)




	bool insert_sorted_dag_list(list<dag> & dag_list1, const dag & dag1){

		bool hasContained = false;

		bool hasInserted = false;

		list<dag>::iterator it1;
		for (it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
			//case 1
			if(*it1 < dag1 ){
				dag_list1.insert(it1, dag1);
				hasInserted = true;
				break;

			}
			//case 2
			else if(*it1 > dag1){
				continue;
			}
			//dag1 == *it1
			else{
				//case 3.1
				if(is_same_dag(*it1, dag1)){
					hasContained = true;

					break;
				}
				//case 3.2: dag1 == *it1 but they are not the same DAG
				else{
					continue;
				}

			}
		}//for (it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++)


		if(hasContained == false && hasInserted == false){
			dag_list1.push_back(dag1);
			hasInserted = true;

		}
		return hasInserted;


	}//bool insert_sorted_dag_list(list<dag> & dag_list1, const dag & dag1)





	//HR: get the sum of the scores of the sampled order stored in t_order_list1
	double get_sum_order_scores(list<total_order> & t_order_list1){
		double sum_scores = LOG_ZERO;
		double temp_order_scores;
		for (list<total_order>::iterator it1 = t_order_list1.begin(); it1 != t_order_list1.end(); it1++){
			temp_order_scores = (*it1).scores;
			logAdd_New(sum_scores, temp_order_scores);
		}
		return sum_scores;

	}


	void print_sigma_k_related_vector(const vector<int> & sigma_k_val_vector, const vector<int> & sigma_k_valBP_vector, const vector<double> & sigma_k_prob_vector){
		cout << endl;
		cout << "sigma_k_val_vector.size() = " << sigma_k_val_vector.size() << endl;
		cout << "  \t" << "sigma_k_val \t" << "sigma_k_valBP \t" << "sigma_k_prob" << "sigma_k_orig_prob" <<endl;
		for(unsigned int i = 0; i < sigma_k_val_vector.size(); i++){
			cout << (i+1) << "\t" << sigma_k_val_vector[i] << "\t"
							  << sigma_k_valBP_vector[i] << "\t"
							  << sigma_k_prob_vector[i] << "\t"
							  << exp(sigma_k_prob_vector[i]) << endl;
		}
		cout << endl;

	}//void print_sigma_k_related_vector()





	vector<family> * get_family_vector(const int i, const int U_i,
			const int n, const int max_indegree, const int numOfFam){

		vector<family> * family_vector_ptr = new vector<family>;

		(*family_vector_ptr).reserve(numOfFam);
		get_family_vector_rec(family_vector_ptr, i, U_i, n, max_indegree, 0, 0, 0);

		return family_vector_ptr;
	}//vector<family> * get_family_vector()




	vector<int> * get_raw_family_vector(const int i, const int U_i,
			const int n, const int max_indegree, const int numOfFam){

		vector<int> * raw_family_vector_ptr = new vector<int>;

		(*raw_family_vector_ptr).reserve(numOfFam);
		get_raw_family_vector_rec(raw_family_vector_ptr, i, U_i, n, max_indegree, 0, 0, 0);

		return raw_family_vector_ptr;
	}//vector<int> * get_raw_family_vector()




	void get_family_vector_rec(vector<family> *family_vector_ptr,
				const int i, const int U_i, const int n, const int max_indegree,
				int depth, int G_i, int G_i_size){
		if(depth == n){
			family fam1;
			fam1.child = i;
			fam1.parents = G_i;

			(*family_vector_ptr).push_back(fam1);
			return;
		}

		else{
			get_family_vector_rec(family_vector_ptr, i, U_i, n, max_indegree, depth+1, G_i, G_i_size);


			if( ( (U_i>>depth) & 1) == 1  &&  G_i_size < max_indegree ){ //this is correct

				get_family_vector_rec(family_vector_ptr, i, U_i, n, max_indegree, depth+1, ( G_i | (1<<depth) ), ( G_i_size + 1 ) );

			}
		}
	}//void get_family_vector_rec()




	void get_raw_family_vector_rec(vector<int> *raw_family_vector_ptr,
				const int i, const int U_i, const int n, const int max_indegree,
				int depth, int G_i, int G_i_size){
		if(depth == n){

			(*raw_family_vector_ptr).push_back(G_i);
			return;
		}

		else{
			get_raw_family_vector_rec(raw_family_vector_ptr, i, U_i, n, max_indegree, depth+1, G_i, G_i_size);

			if( ( (U_i>>depth) & 1) == 1  &&  G_i_size < max_indegree ){ //this is correct
				get_raw_family_vector_rec(raw_family_vector_ptr, i, U_i, n, max_indegree, depth+1, ( G_i | (1<<depth) ), ( G_i_size + 1 ) );

			}
		}
	}//void get_raw_family_vector_rec()




	int sample_one_family_index(vector<family> * U_i_fam_vector_ptr, double alpha_U_i){

		double rand_prob = HR_unif_rand();
		double cumu_prob = 0;

		int family_index;
		//For each (i, Pa_i)
		for(family_index = 0; family_index < (int) (* U_i_fam_vector_ptr).size(); family_index++){
			double curr_fam_scores;

			int curr_fam_child = (* U_i_fam_vector_ptr)[family_index].child;
			int curr_fam_parents = (* U_i_fam_vector_ptr)[family_index].parents;
			curr_fam_scores = beta[curr_fam_child][curr_fam_parents];
			cumu_prob += exp( curr_fam_scores - alpha_U_i );

			if(rand_prob <= cumu_prob){
				return family_index;
			}
		}

		//Impossible to get this point
		cout << "Impossible to get this point in sample_one_family_index()" << endl;
		return family_index;

	}//int sample_one_family_index(vector<family> * U_i_fam_vector_ptr, double alpha_U_i)




	int sample_one_family_index(const vector<double> * prob_cumul_vector_ptr){
		double rand_prob = HR_unif_rand();

		//binary search
		int lower_ind = 0;
		int upper_ind = (*prob_cumul_vector_ptr).size() - 1;
		int mid_ind = (lower_ind + upper_ind) / 2;
		while(!( (rand_prob <= (*prob_cumul_vector_ptr)[mid_ind] && mid_ind == 0)
					|| (rand_prob <= (*prob_cumul_vector_ptr)[mid_ind] && (*prob_cumul_vector_ptr)[mid_ind - 1] < rand_prob)) ){
			if(rand_prob > (*prob_cumul_vector_ptr)[mid_ind]){
				lower_ind = mid_ind + 1;
			}

			else{
				upper_ind = mid_ind - 1;
			}

			assert(lower_ind <= upper_ind);
			mid_ind = (lower_ind + upper_ind) / 2;
		}
		return mid_ind;

	}//int sample_one_family_index(vector<double> * prob_cumul_vector_ptr)




	int sample_one_family_index(const double * prob_cumul_vector_ptr, const int prob_cumul_vector_ptr_length){
		double rand_prob = HR_unif_rand();

		//binary search
		int lower_ind = 0;

		int upper_ind = prob_cumul_vector_ptr_length - 1;
		int mid_ind = (lower_ind + upper_ind) / 2;
		while(!( (rand_prob <= prob_cumul_vector_ptr[mid_ind] && mid_ind == 0)
					|| (rand_prob <= prob_cumul_vector_ptr[mid_ind] && prob_cumul_vector_ptr[mid_ind - 1] < rand_prob)) ){
			if(rand_prob > prob_cumul_vector_ptr[mid_ind]){
				lower_ind = mid_ind + 1;
			}

			else{
				upper_ind = mid_ind - 1;
			}
			//invariant:
			assert(lower_ind <= upper_ind);
			mid_ind = (lower_ind + upper_ind) / 2;
		}
		return mid_ind;

	}//int sample_one_family_index(const double * prob_cumul_vector_ptr, const int prob_cumul_vector_ptr_length)




	vector<double> * get_prob_cumul_vector(vector<family> * U_i_fam_vector_ptr_sigma_i_U_sigma_i, double alpha_sigma_i_U_sigma_i){

		vector<double> * prob_cumul_vector_ptr = new vector<double>;
		(* prob_cumul_vector_ptr).reserve((* U_i_fam_vector_ptr_sigma_i_U_sigma_i).size());
		double cumu_prob = 0;
		for(int family_index = 0; family_index < (int) (* U_i_fam_vector_ptr_sigma_i_U_sigma_i).size(); family_index++){
			double curr_fam_scores;

			int curr_fam_child = (* U_i_fam_vector_ptr_sigma_i_U_sigma_i)[family_index].child;
			int curr_fam_parents = (* U_i_fam_vector_ptr_sigma_i_U_sigma_i)[family_index].parents;
			curr_fam_scores = beta[curr_fam_child][curr_fam_parents];
			cumu_prob += exp( curr_fam_scores - alpha_sigma_i_U_sigma_i );
			(*prob_cumul_vector_ptr).push_back(cumu_prob);
		}
		assert(fabs(cumu_prob - 1) < 1e-5);
		return prob_cumul_vector_ptr;
	}//vector<double> * get_prob_cumul_vector(vector<family> * U_i_fam_vector_ptr_sigma_i_U_sigma_i, double alpha_U_i)




	vector<double> * get_prob_cumul_vector(vector<int> * U_i_fam_vector_ptr_sigma_i_U_sigma_i, double alpha_sigma_i_U_sigma_i, int sigma_i){
		vector<double> * prob_cumul_vector_ptr = new vector<double>;
		(* prob_cumul_vector_ptr).reserve((* U_i_fam_vector_ptr_sigma_i_U_sigma_i).size());
		double cumu_prob = 0;
		for(int family_index = 0; family_index < (int) (* U_i_fam_vector_ptr_sigma_i_U_sigma_i).size(); family_index++){
			double curr_fam_scores;

			int curr_fam_child = sigma_i;

			int curr_fam_parents = (* U_i_fam_vector_ptr_sigma_i_U_sigma_i)[family_index];
			curr_fam_scores = beta[curr_fam_child][curr_fam_parents];
			cumu_prob += exp( curr_fam_scores - alpha_sigma_i_U_sigma_i );
			(*prob_cumul_vector_ptr).push_back(cumu_prob);
		}
		assert(fabs(cumu_prob - 1) < 1e-5);
		return prob_cumul_vector_ptr;
	}//vector<double> * get_prob_cumul_vector(vector<int> * U_i_fam_vector_ptr_sigma_i_U_sigma_i, double alpha_sigma_i_U_sigma_i, int sigma_i)




	double * get_prob_cumul_array(vector<int> * U_i_fam_vector_ptr_sigma_i_U_sigma_i, double alpha_sigma_i_U_sigma_i, int sigma_i){

		double * prob_cumul_array_ptr = new double [(* U_i_fam_vector_ptr_sigma_i_U_sigma_i).size()];
		double cumu_prob = 0;
		for(int family_index = 0; family_index < (int) (* U_i_fam_vector_ptr_sigma_i_U_sigma_i).size(); family_index++){
			double curr_fam_scores;

			int curr_fam_child = sigma_i;

			int curr_fam_parents = (* U_i_fam_vector_ptr_sigma_i_U_sigma_i)[family_index];
			curr_fam_scores = beta[curr_fam_child][curr_fam_parents];
			cumu_prob += exp( curr_fam_scores - alpha_sigma_i_U_sigma_i );

			prob_cumul_array_ptr[family_index] = cumu_prob;
		}
		assert(fabs(cumu_prob - 1) < 1e-5);
		return prob_cumul_array_ptr;
	}//double * get_prob_cumul_array(vector<int> * U_i_fam_vector_ptr_sigma_i_U_sigma_i, double alpha_sigma_i_U_sigma_i, int sigma_i)






	dag sample_one_dag(const total_order & order1, vector<family> *** U_i_fam_vector_ptr){
		int sigma_i;
		int U_sigma_i;
		int family_index;
		//n is num of nodes
		int n = order1.length;
		dag sampled_dag(n);
		for(int ind = 0; ind < order1.length; ind++ ){
			sigma_i = order1.torder[ind];
			U_sigma_i = order1.torderBP[ind];
			family_index = sample_one_family_index(U_i_fam_vector_ptr[sigma_i][U_sigma_i], alpha[sigma_i][U_sigma_i]);
			family sampled_family = (* U_i_fam_vector_ptr[sigma_i][U_sigma_i])[family_index];

			sampled_dag.parents[sampled_family.child] = sampled_family.parents;

			sampled_dag.scores += beta[sampled_family.child][sampled_family.parents];


		}
//		sampled_dag.print();
		return sampled_dag;
	}//dag sample_one_dag()




	dag sample_one_dag(const total_order & order1, vector<family> *** U_i_fam_vector_ptr,
			vector<double> *** U_i_fam_prob_cumul_vector_ptr){
		int sigma_i;
		int U_sigma_i;
		int family_index;
		//n is num of nodes
		int n = order1.length;
		dag sampled_dag(n);
		for(int ind = 0; ind < order1.length; ind++ ){
			sigma_i = order1.torder[ind];
			U_sigma_i = order1.torderBP[ind];

			family_index = sample_one_family_index(U_i_fam_prob_cumul_vector_ptr[sigma_i][U_sigma_i]);
			family sampled_family = (* U_i_fam_vector_ptr[sigma_i][U_sigma_i])[family_index];

			sampled_dag.parents[sampled_family.child] = sampled_family.parents;

			sampled_dag.scores += beta[sampled_family.child][sampled_family.parents];

		}
//		sampled_dag.print();
		return sampled_dag;
	}//dag sample_one_dag()



	dag sample_one_dag(const total_order & order1, vector<int> *** U_i_fam_vector_ptr,
			vector<double> *** U_i_fam_prob_cumul_vector_ptr){
		int sigma_i;
		int U_sigma_i;
		int family_index;
		//n is num of nodes
		int n = order1.length;
		dag sampled_dag(n);
		for(int ind = 0; ind < order1.length; ind++ ){
			sigma_i = order1.torder[ind];
			U_sigma_i = order1.torderBP[ind];

			family_index = sample_one_family_index(U_i_fam_prob_cumul_vector_ptr[sigma_i][U_sigma_i]);

			int sampled_parents = (* U_i_fam_vector_ptr[sigma_i][U_sigma_i])[family_index];

			sampled_dag.parents[sigma_i] = sampled_parents;

			sampled_dag.scores += beta[sigma_i][sampled_parents];

		}
//		sampled_dag.print();
		return sampled_dag;
	}//dag sample_one_dag()




	dag sample_one_dag(const total_order & order1, vector<int> *** U_i_fam_vector_ptr,
			double *** U_i_fam_prob_cumul_vector_ptr){
		int sigma_i;
		int U_sigma_i;
		int family_index;
		//n is num of nodes
		int n = order1.length;
		dag sampled_dag(n);
		for(int ind = 0; ind < order1.length; ind++ ){
			sigma_i = order1.torder[ind];
			U_sigma_i = order1.torderBP[ind];

			family_index = sample_one_family_index(U_i_fam_prob_cumul_vector_ptr[sigma_i][U_sigma_i], (* U_i_fam_vector_ptr[sigma_i][U_sigma_i]).size());

			int sampled_parents = (* U_i_fam_vector_ptr[sigma_i][U_sigma_i])[family_index];

			sampled_dag.parents[sigma_i] = sampled_parents;

			sampled_dag.scores += beta[sigma_i][sampled_parents];

		}
//		sampled_dag.print();
		return sampled_dag;
	}//dag sample_one_dag()



	//HR:
	//	max_indegree is max indegree
	int getNumbOfFam(int max_indegree, int post_ord){
		int numOfFam;
		if(post_ord <= max_indegree){
			//2^post_ord
			numOfFam = (int) pow(2.0, post_ord);
			return numOfFam;
		}
		else{
			double curr_comb = 1;
			double temp_numOfFam = 1;

			for(double denom = 1, numer = post_ord; denom <= max_indegree; denom++, numer--){

				curr_comb = curr_comb * numer / denom;

				temp_numOfFam += curr_comb;
			}
			numOfFam = (int) temp_numOfFam;
			return numOfFam;
		}
	}//int getNumbOfFam(int max_indegree, int post_ord)





		//HR: 	to use the following function to recycle the part of memory when family_setup_count > Max_family_setup_num
		void recycleMemory(
			vector<int> * * * U_i_fam_vector_ptr,
			double * * * U_i_fam_prob_cumul_vector_ptr,
			int & U_i_fam_vector_setup_count,
			int & family_setup_count,
			short * * U_i_fam_vector_used_count,
			const double Max_family_setup_num,
			const double Min_released_family_setup_num,
			const int nh){
		//need to delete some cache if it occupies too much memory (when family_setup_count >= Max_family_setup_num)
		//if(family_setup_count > Max_family_setup_num ){
			//cout << "order_num = " << order_num << endl;
//			cout << "U_i_fam_vector_setup_count = " << U_i_fam_vector_setup_count << endl;
//			cout << "family_setup_count = " << family_setup_count << endl;
//			cout << "family_setup_count > Max_family_setup_num (" << Max_family_setup_num  << ")" << endl;


			int U_i_fam_vector_used_num[7];

			U_i_fam_vector_used_num[1] = 0;
			U_i_fam_vector_used_num[2] = 0;
			U_i_fam_vector_used_num[3] = 0;
			U_i_fam_vector_used_num[4] = 0;
			U_i_fam_vector_used_num[5] = 0;

			//the number of the used U_i_fam_vector is >= 6
			//int U_i_fam_vector_used_ge6_num = 0;
			U_i_fam_vector_used_num[6] = 0;


			int U_i_fam_vector_used_fam_num[7];


			U_i_fam_vector_used_fam_num[1] = 0;
			U_i_fam_vector_used_fam_num[2] = 0;
			U_i_fam_vector_used_fam_num[3] = 0;
			U_i_fam_vector_used_fam_num[4] = 0;
			U_i_fam_vector_used_fam_num[5] = 0;

			//the number of the used the families pointed by the U_i_fam_vector is >= 6
			//int U_i_fam_vector_used_ge6_fam_num = 0;
			U_i_fam_vector_used_fam_num[6] = 0;

			for (int ch = 0; ch < nh; ch++) {
				for (int fam = 0; fam < (1 << nh); fam++) {
					if(U_i_fam_vector_used_count[ch][fam] == 1){
						U_i_fam_vector_used_num[1]++;
						U_i_fam_vector_used_fam_num[1] += (*U_i_fam_vector_ptr[ch][fam]).size();
					}
					else if(U_i_fam_vector_used_count[ch][fam] == 2){
						U_i_fam_vector_used_num[2]++;
						U_i_fam_vector_used_fam_num[2] += (*U_i_fam_vector_ptr[ch][fam]).size();
					}
					else if(U_i_fam_vector_used_count[ch][fam] == 3){
						U_i_fam_vector_used_num[3]++;
						U_i_fam_vector_used_fam_num[3] += (*U_i_fam_vector_ptr[ch][fam]).size();
					}
					else if(U_i_fam_vector_used_count[ch][fam] == 4){
						U_i_fam_vector_used_num[4]++;
						U_i_fam_vector_used_fam_num[4] += (*U_i_fam_vector_ptr[ch][fam]).size();
					}
					else if(U_i_fam_vector_used_count[ch][fam] == 5){
						U_i_fam_vector_used_num[5]++;
						U_i_fam_vector_used_fam_num[5] += (*U_i_fam_vector_ptr[ch][fam]).size();
					}
					else if(U_i_fam_vector_used_count[ch][fam] >= 6){
						U_i_fam_vector_used_num[6]++;
						U_i_fam_vector_used_fam_num[6] += (*U_i_fam_vector_ptr[ch][fam]).size();
					}
				}//for (int fam = 0; fam < (1 << nh); fam++)
			}//for (int ch = 0; ch < nh; ch++)


			int used_num;
			double accum_family_setup_num = 0;
			for(used_num = 1; used_num <= 6; used_num++){
				accum_family_setup_num += U_i_fam_vector_used_fam_num[used_num];
				if(accum_family_setup_num > Min_released_family_setup_num){
					break;
				}
			}
//			cout << "used_num = " << used_num << endl;
//			cout << "accum_family_setup_num = " << accum_family_setup_num << endl;
			if(used_num <= 5){
				for(int ch = 0; ch < nh; ch++){
					for(int fam = 0; fam < (1<<nh); fam++){
						if(U_i_fam_vector_ptr[ch][fam] != NULL && U_i_fam_vector_used_count[ch][fam] <= used_num){
							U_i_fam_vector_setup_count--;
							U_i_fam_vector_used_count[ch][fam] = 0;
							family_setup_count -= (*U_i_fam_vector_ptr[ch][fam]).size();

							(*U_i_fam_vector_ptr[ch][fam]).clear();
							delete U_i_fam_vector_ptr[ch][fam];
							U_i_fam_vector_ptr[ch][fam] = NULL;

							delete [] U_i_fam_prob_cumul_vector_ptr[ch][fam];
							U_i_fam_prob_cumul_vector_ptr[ch][fam] = NULL;
						}

					}
				}

			}//if(used_num <= 5)
			else{
				for(int ch = 0; ch < nh; ch++){
					for(int fam = 0; fam < (1<<nh); fam++){
						if(U_i_fam_vector_ptr[ch][fam] != NULL){
							U_i_fam_vector_setup_count--;
							U_i_fam_vector_used_count[ch][fam] = 0;
							family_setup_count -= (*U_i_fam_vector_ptr[ch][fam]).size();

							(*U_i_fam_vector_ptr[ch][fam]).clear();
							delete U_i_fam_vector_ptr[ch][fam];
							U_i_fam_vector_ptr[ch][fam] = NULL;

							delete [] U_i_fam_prob_cumul_vector_ptr[ch][fam];
							U_i_fam_prob_cumul_vector_ptr[ch][fam] = NULL;
						}

					}
				}

			}
		//}//if(family_setup_count > Max_family_setup_num )

		}//void recycleMemory()





	void read_dag_list_from_files_1(list<dag> & dag_list1,
			const int num_of_runs) {

		for (int run_index = 0; run_index < num_of_runs; run_index++) {
			//test read
			FILE* fin;

			string inStr;

			inStr.assign("dag_list_");
			char listNoStr2[40];
			//int v2 = 0;
			//sprintf(listNoStr2,"%d", v2);
			sprintf(listNoStr2, "%d", run_index);
			inStr.append(listNoStr2);

			fin = fopen(inStr.c_str(), "rb");
			read_DAGs_from_file(dag_list1, fin);

			if(run_index % 4 == 1){
				dag_list1.sort(descend_order_comparison);
				unique_HR(dag_list1);
			}
			//

			fclose(fin);

		}//for(int run_index = 0; run_index < num_of_runs; run_index++)


		dag_list1.sort(descend_order_comparison);
		unique_HR(dag_list1);


	}//void read_dag_list_from_files_1()



	void read_dag_list_from_files_2(list<dag> & dag_list1,
			const int num_of_runs) {


		for (int run_index = 0; run_index < num_of_runs; run_index++) {
			//test read
			FILE* fin;

			string inStr;

			inStr.assign("dag_list_");
			char listNoStr2[40];
			//int v2 = 0;
			//sprintf(listNoStr2,"%d", v2);
			sprintf(listNoStr2, "%d", run_index);
			inStr.append(listNoStr2);

			fin = fopen(inStr.c_str(), "rb");
			read_DAGs_from_file(dag_list1, fin);

			if(run_index % 4 == 1){

				dag_list1.sort(descend_order_comparison);
				unique_HR(dag_list1);
			}
			//

			fclose(fin);

		}//for(int run_index = 0; run_index < num_of_runs; run_index++)


		dag_list1.sort(strict_descend_order_comparison);
		dag_list1.unique(is_same_dag);



	}//void read_dag_list_from_files_2()


	void read_dag_list_from_files_3(list<dag> & dag_list1,
			const int num_of_runs) {


		for (int run_index = 0; run_index < num_of_runs; run_index++) {
			//test read
			FILE* fin;

			string inStr;

			inStr.assign("dag_list_");
			char listNoStr2[40];
			//int v2 = 0;
			//sprintf(listNoStr2,"%d", v2);
			sprintf(listNoStr2, "%d", run_index);
			inStr.append(listNoStr2);

			fin = fopen(inStr.c_str(), "rb");
			read_DAGs_from_file(dag_list1, fin);

			//To reduce dag_list1.size() (the memory requirement)
			if(run_index % 4 == 1){

				dag_list1.sort(descend_order_comparison);
				unique_HR(dag_list1);
			}
			//

			fclose(fin);

		}//for(int run_index = 0; run_index < num_of_runs; run_index++)


		dag_list1.sort(strict_descend_order_comparison);
		unique_HR(dag_list1);

	}//void read_dag_list_from_files_3()



	void read_dag_list_from_files_4(list<dag> & dag_list1,
			const int num_of_runs) {

		list<dag> dag_list_temp;

		for (int run_index = 0; run_index < num_of_runs; run_index++) {
			//test read
			FILE* fin;

			string inStr;

			inStr.assign("dag_list_");
			char listNoStr2[40];
			//int v2 = 0;
			//sprintf(listNoStr2,"%d", v2);
			sprintf(listNoStr2, "%d", run_index);
			inStr.append(listNoStr2);

			fin = fopen(inStr.c_str(), "rb");
			read_DAGs_from_file(dag_list_temp, fin);

			if(run_index % 4 == 1){

				dag_list_temp.sort(descend_order_comparison);
				unique_HR(dag_list_temp);
			}
			//

			fclose(fin);

		}//for(int run_index = 0; run_index < num_of_runs; run_index++)



		DAGHashSet DAGHashSet1;
		DAGHashSet::iterator iter1;
		for (list<dag>::iterator it1 = dag_list_temp.begin(); it1 != dag_list_temp.end(); it1++){
			dag samp_dag = *it1;
			if((iter1 = DAGHashSet1.find(samp_dag)) == DAGHashSet1.end()){
				dag_list1.push_back(samp_dag);
				DAGHashSet1.insert(samp_dag);
			}
		}

	}//void read_dag_list_from_files_4()





	bool has_unprocessed_dags(dag ** curr_comp_dag_array_ptr, const int num_of_runs){
		for (int run_index = 0; run_index < num_of_runs; run_index++){
			if(curr_comp_dag_array_ptr[run_index] != NULL){
				return true;
			}
		}
		return false;
	}//bool has_unprocessed_dags()



	int get_curr_best_dag_index(dag ** curr_comp_dag_array_ptr, const int num_of_runs){
		int best_dag_index = -1;

		int init_run_index;
		for (init_run_index = 0; init_run_index < num_of_runs; init_run_index++){
			if(curr_comp_dag_array_ptr[init_run_index] != NULL){
				best_dag_index = init_run_index;

				break;
			}
		}

		for (int run_index = init_run_index + 1; run_index < num_of_runs; run_index++){

			if(curr_comp_dag_array_ptr[run_index] != NULL && strict_descend_order_comparison(*(curr_comp_dag_array_ptr[run_index]), *(curr_comp_dag_array_ptr[best_dag_index]))){
				best_dag_index = run_index;
			}
		}
		return best_dag_index;
	}//int get_curr_best_dag_index()






	void read_DAGs_from_file_with_threshold_upperCount(
			list<dag> & dag_list_curr_run,
			const int run_index,
			FILE * fin,
			const double DAGScoresThreshold,
			int & availNumOfDAGs,
			DAGHashSet & DAGHashSet1,
			DAGHashSet * DAGHashSetArrWithMarginalScores){

		while(availNumOfDAGs > 0){
			dag dag1(MAX_NUM_OF_VARS);
			read_a_DAG_from_file(dag1, fin);
			availNumOfDAGs--;

			if((run_index > 0
					&& DAGHashSetArrWithMarginalScores[run_index - 1].find(dag1)
							== DAGHashSetArrWithMarginalScores[run_index - 1].end())
				|| run_index == 0){

					if(DAGHashSet1.find(dag1) == DAGHashSet1.end()){
						dag_list_curr_run.push_back(dag1);
						DAGHashSet1.insert(dag1);
						if(fabs(dag1.scores - DAGScoresThreshold) < 0.01){
							DAGHashSetArrWithMarginalScores[run_index].insert(dag1);
						}
					}

			}

			if(dag1.scores <= DAGScoresThreshold){
				break;
			}

		}//while(availNumOfDAGs > 0)

	}//void read_DAGs_from_file_with_threshold_upperCount()


	void read_DAGs_from_file_with_upperCount(
			list<dag> & dag_list_curr_run,
			const int run_index,
			FILE * fin,
			int & availNumOfDAGs,
			DAGHashSet & DAGHashSet1,
			DAGHashSet * DAGHashSetArrWithMarginalScores){

		while(availNumOfDAGs > 0){
			dag dag1(MAX_NUM_OF_VARS);
			read_a_DAG_from_file(dag1, fin);
			availNumOfDAGs--;

			if((run_index > 0
					&& DAGHashSetArrWithMarginalScores[run_index - 1].find(dag1)
							== DAGHashSetArrWithMarginalScores[run_index - 1].end())
				|| run_index == 0){

					if(DAGHashSet1.find(dag1) == DAGHashSet1.end()){
						dag_list_curr_run.push_back(dag1);
						DAGHashSet1.insert(dag1);

					}

			}

		}//while(availNumOfDAGs > 0)

	}//void read_DAGs_from_file_with_upperCount()




	//compute edge posterior by on-the-fly using the DAGs in the hard-disk
	//Using thresholds and Hashtables
	void computeEdgePosteriorInDirectNew(const int nh,
			const int num_of_runs,
			double ** edge_posterior,
			double & sum_posterior_dag_list,
			double * DAGScoresThresholdPtr) {

		const int n = nh;

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				edge_posterior[i][j] = 0;
			}
		}

		double ** log_edge_posterior;
		log_edge_posterior = new double * [n];
		for (int i = 0; i < n; i++){
			log_edge_posterior[i] = new double[n];
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				log_edge_posterior[i][j] = LOG_ZERO;
			}
		}


		sum_posterior_dag_list = LOG_ZERO;


		FILE ** fin_ptr = new FILE *[num_of_runs];
		for (int run_index = 0; run_index < num_of_runs; run_index++) {


			string inStr;

			inStr.assign("dag_list_");
			char listNoStr2[40];
			//int v2 = 0;
			//sprintf(listNoStr2,"%d", v2);
			sprintf(listNoStr2, "%d", run_index);
			inStr.append(listNoStr2);

			fin_ptr[run_index] = fopen(inStr.c_str(), "rb");


		}//for(int run_index = 0; run_index < num_of_runs; run_index++)

		int * dag_list_size_array_ptr = new int [num_of_runs];
		for (int run_index = 0; run_index < num_of_runs; run_index++) {
			fread(&(dag_list_size_array_ptr[run_index]), sizeof(int), 1, fin_ptr[run_index]);
		}

		for (int run_index = 0; run_index < num_of_runs; run_index++) {
//			cerr << "dag_list_size_array_ptr[" << run_index << "] = "
//					<< dag_list_size_array_ptr[run_index] << endl;
		}


		DAGHashSet DAGHashSet1;


		DAGHashSet * DAGHashSetArrWithMarginalScores = new DAGHashSet[num_of_runs - 1];

		//for each segment across all the files
		for (int run_index = 0; run_index <= num_of_runs - 1; run_index++) {
			list<dag> dag_list_curr_run;
			DAGHashSet1.clear();

			if(run_index <= num_of_runs - 2){

				for(int run_index_new = 0; run_index_new <= num_of_runs - 1; run_index_new++){
					read_DAGs_from_file_with_threshold_upperCount(
							dag_list_curr_run,
							run_index,
							fin_ptr[run_index_new],
							DAGScoresThresholdPtr[run_index],
							dag_list_size_array_ptr[run_index_new],
							DAGHashSet1,
							DAGHashSetArrWithMarginalScores);
				}
			}

			else{

				for(int run_index_new = 0; run_index_new <= num_of_runs - 1; run_index_new++){
					read_DAGs_from_file_with_upperCount(
							dag_list_curr_run,
							run_index,
							fin_ptr[run_index_new],
							dag_list_size_array_ptr[run_index_new],
							DAGHashSet1,
							DAGHashSetArrWithMarginalScores);
				}

			}
			if(run_index >= 2){
				DAGHashSetArrWithMarginalScores[run_index - 2].clear();
			}

			double temp_dag_scores;
			for (list<dag>::const_iterator it1 = dag_list_curr_run.begin(); it1
					!= dag_list_curr_run.end(); it1++) {

				temp_dag_scores = (*it1).scores;
				logAdd_New(sum_posterior_dag_list, temp_dag_scores);

				incre_edge_posterior_by_posterior_of_dags(log_edge_posterior, (*it1));
			}

		}//for each segment across all the files
//		printf("sum_posterior_dag_list = %18.6f\n", sum_posterior_dag_list);


		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				edge_posterior[i][j] = exp(log_edge_posterior[i][j] - sum_posterior_dag_list);
			}
		}

		for(int i = 0; i < n; i++){
			delete [] log_edge_posterior[i];
		}
		delete [] log_edge_posterior;


		delete [] dag_list_size_array_ptr;


		for (int run_index = 0; run_index < num_of_runs; run_index++) {
			fclose(fin_ptr[run_index]);
		}
		delete[] fin_ptr;

		DAGHashSet1.clear();

		for(int run_index = 0; run_index <= num_of_runs - 2; run_index++){
			DAGHashSetArrWithMarginalScores[run_index].clear();
		}

		delete [] DAGHashSetArrWithMarginalScores;
	}//void computeEdgePosteriorInDirectNew()






	void get_edge_poster_by_num_of_dags(double ** edge_posterior, int n, const list<dag> & dag_list1){

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				edge_posterior[i][j] = 0;
			}
		}

		assert(n == (dag_list1.front()).length);

		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
			incre_edge_count_by_num_of_dags(edge_posterior, (*it1));
		}

		int num_of_dags = dag_list1.size();
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				edge_posterior[i][j] = edge_posterior[i][j] / ((double) num_of_dags);
			}
		}
	}//void get_edge_poster_by_num_of_dags(double ** edge_posterior, int n, list<dag> & dag_list1)




	void get_edge_poster_by_posterir_of_dags(double ** edge_posterior, int n, const list<dag> & dag_list1){

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				edge_posterior[i][j] = 0;
			}
		}

		double ** log_edge_posterior;
		log_edge_posterior = new double * [n];
		for (int i = 0; i < n; i++){
			log_edge_posterior[i] = new double[n];
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				log_edge_posterior[i][j] = LOG_ZERO;
			}
		}

		assert(n == (dag_list1.front()).length);

		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
			incre_edge_posterior_by_posterior_of_dags(log_edge_posterior, (*it1));
		}

		double sum_posterior_dag_list = LOG_ZERO;
		double temp_dag_scores;
		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){

			temp_dag_scores = (*it1).scores;
			logAdd_New(sum_posterior_dag_list, temp_dag_scores);
		}


		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				edge_posterior[i][j] = exp(log_edge_posterior[i][j] - sum_posterior_dag_list);
			}
		}

		for(int i = 0; i < n; i++){
			delete [] log_edge_posterior[i];
		}
		delete [] log_edge_posterior;



	}//void get_edge_poster_by_posterir_of_dags(double ** edge_posterior, int n, const list<dag> & dag_list1)




	void get_path_poster_by_posterir_of_dags(double ** path_posterior, int n, const list<dag> & dag_list1){

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				path_posterior[i][j] = 0;
			}
		}


		double ** log_path_posterior;
		log_path_posterior = new double * [n];
		for (int i = 0; i < n; i++){
			log_path_posterior[i] = new double[n];
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				log_path_posterior[i][j] = LOG_ZERO;
			}
		}

		assert(n == (dag_list1.front()).length);

		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
			incre_path_posterior_by_posterior_of_dags(log_path_posterior, (*it1));
		}

		double sum_posterior_dag_list = LOG_ZERO;
		double temp_dag_scores;
		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){

			temp_dag_scores = (*it1).scores;
			logAdd_New(sum_posterior_dag_list, temp_dag_scores);
		}


		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				path_posterior[i][j] = exp(log_path_posterior[i][j] - sum_posterior_dag_list);
			}
		}

		for(int i = 0; i < n; i++){
			delete [] log_path_posterior[i];
		}
		delete [] log_path_posterior;


	}//void get_path_poster_by_posterir_of_dags(double ** path_posterior, int n, const list<dag> & dag_list1)





	void get_limited_leng_path_poster_by_posterir_of_dags(double ** limited_leng_path_posterior, int n, const list<dag> & dag_list1, int limited_leng){

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				limited_leng_path_posterior[i][j] = 0;
			}
		}

		//use log version first
		double ** log_limited_leng_path_posterior;
		log_limited_leng_path_posterior = new double * [n];
		for (int i = 0; i < n; i++){
			log_limited_leng_path_posterior[i] = new double[n];
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				log_limited_leng_path_posterior[i][j] = LOG_ZERO;
			}
		}

		assert(n == (dag_list1.front()).length);

		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
			incre_limited_leng_path_posterior_by_posterior_of_dags(log_limited_leng_path_posterior, (*it1), limited_leng);
		}

		double sum_posterior_dag_list = LOG_ZERO;
		double temp_dag_scores;
		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){

			temp_dag_scores = (*it1).scores;
			logAdd_New(sum_posterior_dag_list, temp_dag_scores);
		}


		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				limited_leng_path_posterior[i][j] = exp(log_limited_leng_path_posterior[i][j] - sum_posterior_dag_list);
			}
		}

		for(int i = 0; i < n; i++){
			delete [] log_limited_leng_path_posterior[i];
		}
		delete [] log_limited_leng_path_posterior;


	}//void get_limited_leng_path_poster_by_posterir_of_dags()





	void get_combined_two_path_poster_by_posterir_of_dags(double *** combined_two_path_posterior, int n, const list<dag> & dag_list1){
		//initialize combined_two_path_posterior
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++){
					combined_two_path_posterior[i][j][k] = 0;
				}
			}
		}

		//use log version first
		double *** log_combined_two_path_posterior;
		log_combined_two_path_posterior = new double ** [n];
		for (int i = 0; i < n; i++){
			log_combined_two_path_posterior[i] = new double * [n];
			for(int j = 0; j < n; j++){
				log_combined_two_path_posterior[i][j] = new double [n];
			}
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++){
					log_combined_two_path_posterior[i][j][k] = LOG_ZERO;
				}
			}
		}

		assert(n == (dag_list1.front()).length);

		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
			incre_combined_two_path_posterior_by_posterior_of_dags(log_combined_two_path_posterior, (*it1));
		}

		double sum_posterior_dag_list = LOG_ZERO;
		double temp_dag_scores;
		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){

			temp_dag_scores = (*it1).scores;
			logAdd_New(sum_posterior_dag_list, temp_dag_scores);
		}


		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++){
					combined_two_path_posterior[i][j][k] = exp(log_combined_two_path_posterior[i][j][k] - sum_posterior_dag_list);
				}
			}
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				delete [] log_combined_two_path_posterior[i][j];
			}
			delete [] log_combined_two_path_posterior[i];
		}
		delete [] log_combined_two_path_posterior;


	}//void get_combined_two_path_poster_by_posterir_of_dags()



	//HR: get the estimated posterior of each combination of two paths (i ~~> j and i ~~> k and j ~= k) by weighting the posterior of dags in dag_list1
	//	  combined_two_f4_path_posterior is the pointer to n * n * n matrix,
	//	  combined_two_f4_path_posterior[i][j][k] records the posterior of directed path i ~~> j and i ~~> k and j ~= k
	//    the same as Matlab: C(i, j, k) = 1 iff R(i, j) = 1 and R(i, k) = 1 and j ~= k
	void get_combined_two_f4_path_poster_by_posterir_of_dags(double *** combined_two_f4_path_posterior, int n, const list<dag> & dag_list1){
		//initialize combined_two_f4_path_posterior
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++){
					combined_two_f4_path_posterior[i][j][k] = 0;
				}
			}
		}

		//use log version first
		double *** log_combined_two_f4_path_posterior;
		log_combined_two_f4_path_posterior = new double ** [n];
		for (int i = 0; i < n; i++){
			log_combined_two_f4_path_posterior[i] = new double * [n];
			for(int j = 0; j < n; j++){
				log_combined_two_f4_path_posterior[i][j] = new double [n];
			}
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++){
					log_combined_two_f4_path_posterior[i][j][k] = LOG_ZERO;
				}
			}
		}

		assert(n == (dag_list1.front()).length);

		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
			incre_combined_two_f4_path_posterior_by_posterior_of_dags(log_combined_two_f4_path_posterior, (*it1));
		}

		double sum_posterior_dag_list = LOG_ZERO;
		double temp_dag_scores;
		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){

			temp_dag_scores = (*it1).scores;
			logAdd_New(sum_posterior_dag_list, temp_dag_scores);
		}


		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++){
					combined_two_f4_path_posterior[i][j][k] = exp(log_combined_two_f4_path_posterior[i][j][k] - sum_posterior_dag_list);
				}
			}
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				delete [] log_combined_two_f4_path_posterior[i][j];
			}
			delete [] log_combined_two_f4_path_posterior[i];
		}
		delete [] log_combined_two_f4_path_posterior;


	}//void get_combined_two_f4_path_poster_by_posterir_of_dags()



	//HR: get the estimated posterior of each combination of two paths ( i ~~> j and i not ~~> k and i ~= k) by weighting the posterior of dags in dag_list1
	//	  combined_two_f5_path_posterior is the pointer to n * n * n matrix,
	//	  combined_two_f5_path_posterior[i][j][k] records the posterior of directed path i ~~> j and i not ~~> k and i ~= k
	//    the same as Matlab: C(i, j, k) = 1 iff R(i, j) = 1 and R(i, k) = 0 and i ~= k
	void get_combined_two_f5_path_poster_by_posterir_of_dags(double *** combined_two_f5_path_posterior, int n, const list<dag> & dag_list1){
		//initialize
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++){
					combined_two_f5_path_posterior[i][j][k] = 0;
				}
			}
		}

		//use log version first
		double *** log_combined_two_f5_path_posterior;
		log_combined_two_f5_path_posterior = new double ** [n];
		for (int i = 0; i < n; i++){
			log_combined_two_f5_path_posterior[i] = new double * [n];
			for(int j = 0; j < n; j++){
				log_combined_two_f5_path_posterior[i][j] = new double [n];
			}
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++){
					log_combined_two_f5_path_posterior[i][j][k] = LOG_ZERO;
				}
			}
		}

		assert(n == (dag_list1.front()).length);

		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){
			incre_combined_two_f5_path_posterior_by_posterior_of_dags(log_combined_two_f5_path_posterior, (*it1));
		}

		double sum_posterior_dag_list = LOG_ZERO;
		double temp_dag_scores;
		for (list<dag>::const_iterator it1 = dag_list1.begin(); it1 != dag_list1.end(); it1++){

			temp_dag_scores = (*it1).scores;
			logAdd_New(sum_posterior_dag_list, temp_dag_scores);
		}


		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				for(int k = 0; k < n; k++){
					combined_two_f5_path_posterior[i][j][k] = exp(log_combined_two_f5_path_posterior[i][j][k] - sum_posterior_dag_list);
				}
			}
		}

		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				delete [] log_combined_two_f5_path_posterior[i][j];
			}
			delete [] log_combined_two_f5_path_posterior[i];
		}
		delete [] log_combined_two_f5_path_posterior;


	}//void get_combined_two_f5_path_poster_by_posterir_of_dags()





	//HR: increment the count of each edge based on dag1
	void incre_edge_count_by_num_of_dags(double ** edge_posterior, const dag & dag1){
		//for each node v of dag1
		for(int v = 0; v < dag1.length; v++){
			int parents = dag1.parents[v];
			for(int u = 0; u < dag1.length; u++){
				if( ((parents >> u) & 1) == 1 ) {
					edge_posterior[u][v]++;
				}
			}
		}
	}//void incre_edge_count_by_num_of_dags(double ** edge_posterior, dag & dag1)



	//HR: increment the posterior of each edge based on the posterior of dag1
	void incre_edge_posterior_by_posterior_of_dags(double ** log_edge_posterior, const dag & dag1){
		double dag1_scores;

		//for each node v of dag1
		for(int v = 0; v < dag1.length; v++){
			int parents = dag1.parents[v];
			for(int u = 0; u < dag1.length; u++){
				if( ((parents >> u) & 1) == 1 ) {

					dag1_scores = dag1.scores;
					logAdd_New(log_edge_posterior[u][v], dag1_scores);

				}
			}
		}
	}//void incre_edge_posterior_by_posterior_of_dags(double ** edge_posterior, const dag & dag1)




	//HR: increment the posterior of each path based on the posterior of dag1
	void incre_path_posterior_by_posterior_of_dags(double ** log_path_posterior, const dag & dag1){
		double dag1_scores;

		int n = dag1.length;
		bool ** ance_mat1;
		ance_mat1 = new bool *[n];
		for (int i = 0; i < n; i++) {
			ance_mat1[i] = new bool[n];
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				ance_mat1[i][j] = false;
			}
		}
		get_ancestor_matrix(dag1, ance_mat1);


		//for each node v of dag1
		for(int v = 0; v < dag1.length; v++){
			for(int u = 0; u < dag1.length; u++){
				if(ance_mat1[v][u] == true ) {
					dag1_scores = dag1.scores;
					logAdd_New(log_path_posterior[v][u], dag1_scores);

				}
			}
		}

		for (int i = 0; i < n; i++) {
			delete[] ance_mat1[i];
		}
		delete[] ance_mat1;


	}//void incre_path_posterior_by_posterior_of_dags(double ** log_path_posterior, const dag & dag1)



	//HR: increment the posterior of each limited length path based on the posterior of dag1
	void incre_limited_leng_path_posterior_by_posterior_of_dags(double ** log_limited_leng_path_posterior, const dag & dag1, int limited_leng){
		double dag1_scores;

		int n = dag1.length;


		double ** adj_mat1;
		adj_mat1 = new double *[n];
		for (int i = 0; i < n; i++) {
			adj_mat1[i] = new double[n];
		}
		get_adj_matrix(dag1, adj_mat1, n);

		double ** reachability_limited_leng_mat;
		reachability_limited_leng_mat = new double *[n];
		for (int i = 0; i < n; i++) {
			reachability_limited_leng_mat[i] = new double[n];
		}
		get_reachability_limited_leng(adj_mat1, n, reachability_limited_leng_mat, limited_leng);


		//for each node v of dag1
		for(int v = 0; v < dag1.length; v++){
			//check each possible parent node u
			for(int u = 0; u < dag1.length; u++){
				//if there exists u ~~> v with length <= limited_leng
				if(reachability_limited_leng_mat[u][v] > 0 ) {
					//dag1_scores is used to make to match of the type of
					//parameters of logAdd_New(), since dag1 is const dag
					dag1_scores = dag1.scores;
					logAdd_New(log_limited_leng_path_posterior[u][v], dag1_scores);

				}
			}
		}

		for (int i = 0; i < n; i++) {
			delete[] adj_mat1[i];
		}
		delete[] adj_mat1;

		for (int i = 0; i < n; i++) {
			delete[] reachability_limited_leng_mat[i];
		}
		delete[] reachability_limited_leng_mat;


	}//void incre_limited_leng_path_posterior_by_posterior_of_dags()




	//HR: increment the posterior of each combination of two paths based on the posterior of dag1
	//	  combined_two_path_posterior[i][j][k] records the posterior of directed path i ~~> j ~~> k
	void incre_combined_two_path_posterior_by_posterior_of_dags(double *** log_combined_two_path_posterior, const dag & dag1){
		double dag1_scores;

		int n = dag1.length;
		bool ** ance_mat1;
		ance_mat1 = new bool *[n];
		for (int i = 0; i < n; i++) {
			ance_mat1[i] = new bool[n];
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				ance_mat1[i][j] = false;
			}
		}
		get_ancestor_matrix(dag1, ance_mat1);


		for(int i = 0; i < dag1.length; i++){
			for(int j = 0; j < dag1.length; j++){
				for(int k = 0; k < dag1.length; k++){
					//if i <~~ j exists and j <~~ k exists
					if(ance_mat1[i][j] == true && ance_mat1[j][k] == true) {
						dag1_scores = dag1.scores;
						logAdd_New(log_combined_two_path_posterior[k][j][i], dag1_scores);

					}
				}
			}
		}

		for (int i = 0; i < n; i++) {
			delete[] ance_mat1[i];
		}
		delete[] ance_mat1;


	}//void incre_combined_two_path_posterior_by_posterior_of_dags()


	//HR: increment the posterior of each feature f4 based on the posterior of dag1
	//	  combined_two_f4_path_posterior[i][j][k] records the posterior of directed path i ~~> j and i ~~> k and j ~= k
	//	  the same as Matlab: C(i, j, k) = 1 iff R(i, j) = 1 and R(i, k) = 1 and j ~= k
	void incre_combined_two_f4_path_posterior_by_posterior_of_dags(double *** log_combined_two_f4_path_posterior, const dag & dag1){
		double dag1_scores;

		int n = dag1.length;
		bool ** ance_mat1;
		ance_mat1 = new bool *[n];
		for (int i = 0; i < n; i++) {
			ance_mat1[i] = new bool[n];
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				ance_mat1[i][j] = false;
			}
		}
		get_ancestor_matrix(dag1, ance_mat1);


		for(int i = 0; i < dag1.length; i++){
			for(int j = 0; j < dag1.length; j++){
				for(int k = 0; k < dag1.length; k++){

					//if i ~~> j and i ~~> k and j ~= k
					if(ance_mat1[j][i] == true && ance_mat1[k][i] == true && j != k) {
						dag1_scores = dag1.scores;
						logAdd_New(log_combined_two_f4_path_posterior[i][j][k], dag1_scores);

					}
				}
			}
		}

		for (int i = 0; i < n; i++) {
			delete[] ance_mat1[i];
		}
		delete[] ance_mat1;


	}//void incre_combined_two_f4_path_posterior_by_posterior_of_dags()



	//HR: increment the posterior of each feature f5 based on the posterior of dag1
	//	  combined_two_f5_path_posterior[i][j][k] records the posterior of directed path i ~~> j and i not ~~> k and i ~= k
	//	  the same as Matlab: C(i, j, k) = 1 iff R(i, j) = 1 and R(i, k) = 0 and i ~= k
	void incre_combined_two_f5_path_posterior_by_posterior_of_dags(double *** log_combined_two_f5_path_posterior, const dag & dag1){
		double dag1_scores;

		int n = dag1.length;
		bool ** ance_mat1;
		ance_mat1 = new bool *[n];
		for (int i = 0; i < n; i++) {
			ance_mat1[i] = new bool[n];
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				ance_mat1[i][j] = false;
			}
		}
		get_ancestor_matrix(dag1, ance_mat1);


		for(int i = 0; i < dag1.length; i++){
			for(int j = 0; j < dag1.length; j++){
				for(int k = 0; k < dag1.length; k++){

					//if i ~~> j and i not ~~> k and i ~= k
					if(ance_mat1[j][i] == true && ance_mat1[k][i] != true && i != k) {

						dag1_scores = dag1.scores;
						logAdd_New(log_combined_two_f5_path_posterior[i][j][k], dag1_scores);

					}
				}
			}
		}

		for (int i = 0; i < n; i++) {
			delete[] ance_mat1[i];
		}
		delete[] ance_mat1;


	}//void incre_combined_two_f5_path_posterior_by_posterior_of_dags()






	void add_edge_poster_given_order(double ** edge_posterior, const total_order & order1, int n){
		//for each node sigma_i in the order
		for(int i = 0; i < n; i++){

			for(int j = 0; j < i; j++){

				int U_sigma_i_new = 0;

				int sigma_i = order1.torder[i];
				int U_sigma_i = order1.torderBP[i];
				int sigma_j = order1.torder[j];
				U_sigma_i_new = U_sigma_i - (1 << sigma_j);
				edge_posterior[sigma_j][sigma_i] += 1 - exp(alpha[sigma_i][U_sigma_i_new] - alpha[sigma_i][U_sigma_i]);
			}

		}
	}//void add_edge_poster_given_order()


	// alpha[j][S] = sum_T beta[j][T]
	void compute_alpha(int j, int *Vh, int nh){
		// Use fast truncated upward Möbius transform.
		fumt(beta[j], alpha[j], nh, k);
	}


	void compute_alphaNume(int j, int *Vh, int nh){
		fumt(betaNume[j], alphaNume[j], nh, k);
	}


	void compute_alphaInOut(int j, int *Vh, int nh){
		fumt(betaInOut[j], alphaInOut[j], nh, k);
	}


	void compute_g_forward(int nh){

		sub_gf(0, nh, 0);

	}


	void compute_gNume_forward(int nh){

		sub_gfNume(0, nh, 0);

	}


	void sub_gf(int d, int nh, int S){
		//HR:
		//totalOf2 << "sub_gf(" << d <<", " << nh <<", " << S <<")" << endl;

		if (d < nh){

			sub_gf(d+1, nh, S);
			sub_gf(d+1, nh, S | (1 << d));
		}

		else {
			double sum = MARK;
			int T = S, J = 1;


			for (int j = 0; j < nh; j ++){
				if (T & 1){

					double w = alpha[j][S - J] + gf[S - J];

					logAdd_New(sum, w);
				}

				T >>= 1; J <<= 1;
			}

			if (S == 0)
				gf[0] = 0;
			else
				gf[S] = sum;
		}
	}


	//HR:
		void sub_gfNume(int d, int nh, int S){
		//HR:
		//totalOf2 << "sub_gf(" << d <<", " << nh <<", " << S <<")" << endl;

		if (d < nh){

			sub_gfNume(d+1, nh, S);
			sub_gfNume(d+1, nh, S | (1 << d));
		}

		else {
			double sum = MARK;
			int T = S, J = 1;


			for (int j = 0; j < nh; j ++){

				if (T & 1){

					double w = alphaNume[j][S - J] + gfNume[S - J];

					logAdd_New(sum, w);
				}

				T >>= 1; J <<= 1;
			}
			if (S == 0)
				gfNume[0] = 0;
			else
				gfNume[S] = sum;
		}
	}




	void compute_g_backward(int nh){

		sub_gb(0, nh, 0);

	}
	void sub_gb(int d, int nh, int S){
		if (d < nh){
			sub_gb(d+1, nh, S);
			sub_gb(d+1, nh, S | (1 << d));
		}
		else {
			double sum = MARK;
			int T = S, J = 1, complS = (1 << nh) - 1 - S;
			for (int j = 0; j < nh; j ++){
				if (T & 1){

					double w = alpha[j][complS] + gb[S - J];

					logAdd_New(sum, w);
				}

				T >>= 1; J <<= 1;
			}
			if (S == 0)
				gb[0] = 0;
			else
				gb[S] = sum;
		}
	}


	double eval_edge(int i, int j, double *a, double *b, int nh, int k){

		int T = 1 << i;
		return sub_eval_edge(-1, 1, T, i, j, a, b, nh, k);
	}


	double sub_eval_edge(
		int tprev, int d, int T, int i, int j, double *a, double *b, int nh, int k){

		//cerr<<" T:"; print_set(cerr, T);
		//cerr<<": "<<a[T]<<"; "<<b[T]<<endl;

		if (d == k){
			return a[T] + b[T];
		}


		double sum = a[T] + b[T];


		for (int t = tprev + 1; t < nh; t ++){
			if (t == i || t == j) continue;

			int Tnext = T | (1 << t);
			double w = sub_eval_edge(t, d+1, Tnext, i, j, a, b, nh, k);

			logAdd_New(sum, w);
		}
		return sum;
	}





	double eval_edge_mixIndegree(int i, int j, double * Gamma, double * Beta_j, int nh, int k){
		int T = 1 << i;
		return sub_eval_edge_mixIndegree(-1, 1, T, i, j, Gamma, Beta_j, nh, k);
	}


	double sub_eval_edge_mixIndegree(
		int tprev, int d, int T, int i, int j, double * Gamma, double * Beta_j, int nh, int k){

		//cerr<<" T:"; print_set(cerr, T);
		//cerr<<": "<<a[T]<<"; "<<b[T]<<endl;

		if (d == k){
			double sumRes;
			if(Beta_j[T] == MARK){
				cerr<<"\n**** Beta_j[T] == MARK ";
		   	}


			if(Gamma[T] == LOG_ZERO){
				sumRes = LOG_ZERO;
			}

			else if(Gamma[T] <= 0){
				sumRes = Beta_j[T] + Gamma[T];
			}
			else{
				cerr<<"\n**** Gamma[T] > 0 ";
				sumRes = - ( Beta_j[T] - Gamma[T] );

			}
			return sumRes;

		}


		double sum;
		if(Beta_j[T] == MARK){
			cerr<<"\n**** Beta_j[T] == MARK ";
		}
		if(Gamma[T] == LOG_ZERO){
				sum = LOG_ZERO;
		}

		else if(Gamma[T] <= 0){
			sum = Beta_j[T] + Gamma[T];
		}
		else{
			cerr<<"\n**** Gamma[T] > 0 ";
			sum = - ( Beta_j[T] - Gamma[T] );

		}


		for (int t = tprev + 1; t < nh; t ++){

			if (t == i || t == j) continue;


			int Tnext = T | (1 << t);
			double w = sub_eval_edge_mixIndegree(t, d+1, Tnext, i, j, Gamma, Beta_j, nh, k);
			logAddComp(sum, w);
		}
		return sum;
	}



	double eval_edge_mixIndegreeNoLog(int i, int j, double * Gamma, double * Beta_j, int nh, int k){
		int T = 1 << i;
		return sub_eval_edge_mixIndegreeNoLog(-1, 1, T, i, j, Gamma, Beta_j, nh, k);
	}


	double sub_eval_edge_mixIndegreeNoLog(
		int tprev, int d, int T, int i, int j, double * Gamma, double * Beta_j, int nh, int k){

		//cerr<<" T:"; print_set(cerr, T);
		//cerr<<": "<<a[T]<<"; "<<b[T]<<endl;

		if (d == k){

			double product_res;
			if(Beta_j[T] == MARK){
				cerr<<"***Beta_j[T]  has not been initialized" << endl;
				product_res = 0;
			}
			else{
				if(Gamma[T] == 0){
				 	product_res = exp(Beta_j[T]);
				}
				else if(Gamma[T] > 0){
					product_res = exp(Beta_j[T] + log(Gamma[T]));
				}
				else{
					product_res = - exp(Beta_j[T] + log(-Gamma[T]));
				}

			}

			return product_res;

		}


		double sum;
		double product_res;

			if(Beta_j[T] == MARK){
				cerr<<"***Beta_j[T]  has not been initialized" << endl;
				product_res = 0;
			}
			else{
				if(Gamma[T] == 0){
				 	product_res = exp(Beta_j[T]);
				}
				else if(Gamma[T] > 0){
					product_res = exp(Beta_j[T] + log(Gamma[T]));
				}
				else{
					product_res = - exp(Beta_j[T] + log(-Gamma[T]));
				}

			}

		sum = product_res;
		for (int t = tprev + 1; t < nh; t ++){
			if (t == i || t == j) continue;

			int Tnext = T | (1 << t);
			double w = sub_eval_edge_mixIndegreeNoLog(t, d+1, Tnext, i, j, Gamma, Beta_j, nh, k);

			sum += w;
		}
		return sum;
	}




	Model *model;

	double **alpha, **beta;

	double **alphaNume, **betaNume;

	double *gf, *gb;


	double *gfNume;
	int k;


	double **alphaInOut, **betaInOut;
	int * inFeatureEdges;
	int * outFeatureEdges;


};



#endif
