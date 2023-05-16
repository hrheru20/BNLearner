#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include<stdio.h>
#include<stdlib.h>

using namespace std;

class Arguments {
public:
	static char* datafile;
	static char* layeringfile;
	static char* maxindegree;
	static char* model;
	static char* task;
	static char* maxnumrecords;
	
	static char* inFeasFileName;
	static char* outFeasFileName;

	static int option;
	
	static int ADtree;

	static int num_of_sampled_orders;


	//		rebel_para_choice = 1:
	//				Original rebel
	//				lb = model->log_prior(i, S, d) + model->log_lcp(i, S, d);
	//		rebel_para_choice = 2:
	//				rebelOld
	//				lb = model->log_prior(i, S, d) + model->log_lcpHR(i, S, d);
	//		rebel_para_choice = 3:
	//				rebelNew
	//				lb = 0 + model->log_lcpHR(i, S, d);
	static int rebel_para_choice;
	//


	static void init(int argc, char **args){
		for(int i = 1; i < argc; i ++){
			
			if(args[i][0]=='-' && args[i][1]=='d' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::datafile = args[j];
					j ++;
				} 
			}
			else if(args[i][0]=='-' && args[i][1]=='l' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::layeringfile = args[j]; 
					j ++;
				}
			}
			else if(args[i][0]=='-' && args[i][1]=='m' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::maxindegree = args[j]; 
					j ++;
				}
			}
			else if(args[i][0]=='-' && args[i][1]=='u' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::maxnumrecords = args[j]; 
					j ++;
				}
			}
			else if(args[i][0]=='-' && args[i][1]=='M' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::model = args[j]; 
					j ++;
				}
			}
			else if(args[i][0]=='-' && args[i][1]=='T' && args[i][2]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::task = args[j]; 
					j ++;
				}
			}

			else if(args[i][0]=='-' && args[i][1]=='i' && args[i][2]=='n' && args[i][3]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::inFeasFileName = args[j]; 
					j ++;
				}
			}
			else if(args[i][0]=='-' && args[i][1]=='o' && args[i][2]=='u' && args[i][3]=='t' && args[i][4]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::outFeasFileName = args[j]; 
					j ++;
				}
			}


			else if(args[i][0]=='-' && args[i][1]=='o' && args[i][2]=='p' && args[i][3]=='t' && args[i][4]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::option = atoi(args[j]); 
					j ++;
				}
			}

			else if(args[i][0]=='-' && args[i][1]=='a' && args[i][2]=='d' && args[i][3]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::ADtree = atoi(args[j]); 
					j ++;
				}
			}

			else if(args[i][0]=='-' && args[i][1]=='s' && args[i][2]=='o' && args[i][3]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::num_of_sampled_orders = atoi(args[j]);
					j ++;
				}
			}

			else if(args[i][0]=='-' && args[i][1]=='r' && args[i][2]=='e' && args[i][3]=='\0'){
				int j = i + 1;
				while (j < argc && args[j][0]!='-') {
					Arguments::rebel_para_choice = atoi(args[j]);
					j ++;
				}
			}

		}
	
		print_arguments(stderr);
		
	}
	static void print_arguments(FILE *f){
		fprintf(f, " -d Data file:\n");
		fprintf(f, "   %62s\n", Arguments::datafile);	
//		fprintf(f, "   %62s\n", (Arguments::datafile+6));
		//fprintf(f, " -l Layering file:\n");
		//fprintf(f, "   %62s\n", Arguments::layeringfile);
		fprintf(f, " -m Maximum indegree:\n");
		fprintf(f, "   %62s\n", Arguments::maxindegree);
		fprintf(f, " -u Maximum number of data records read:\n");
		fprintf(f, "   %62s\n", Arguments::maxnumrecords);
		//fprintf(f, " -M Model:\n");
		//fprintf(f, "   %62s\n", Arguments::model);
		//fprintf(f, " -T Task (Infer=I, Generate=G):\n");
		//fprintf(f, "   %62s\n", Arguments::task);
		
		fprintf(f, " -opt Option:\n");
		fprintf(f, "   %62d\n", Arguments::option);
		if(Arguments::option == 5){
			fprintf(f, " -in In_Feature File:\n");
			fprintf(f, "   %62s\n", Arguments::inFeasFileName);
			fprintf(f, " -out Out_Feature File:\n");
			fprintf(f, "   %62s\n", Arguments::outFeasFileName);
		}

		fprintf(f, " -ad ADtree:\n");
		fprintf(f, "   %62d\n", Arguments::ADtree);

		if(Arguments::option >= 8 && Arguments::option <= 12){
			fprintf(f, " -so Number of sampled orders:\n");
			fprintf(f, "   %62d\n", Arguments::num_of_sampled_orders);

		}

		if(Arguments::option ==0){
			fprintf(f, " -re Choice of rebel parameter:\n");
			fprintf(f, "   %62d\n", Arguments::rebel_para_choice);

		}

	}
};

char* Arguments::datafile = "testdata.dat";
char* Arguments::layeringfile = "%";
char* Arguments::maxindegree = "3";
char* Arguments::model = "M";
char* Arguments::task = "I";
char* Arguments::maxnumrecords = "999999";


char* Arguments::inFeasFileName = "inFeature.txt";
char* Arguments::outFeasFileName = "outFeature.txt";
int Arguments::option = 5;


int Arguments::ADtree = 0;


int Arguments::num_of_sampled_orders = 0;


int Arguments::rebel_para_choice = 3;


#endif
