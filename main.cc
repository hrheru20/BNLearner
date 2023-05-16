#include<stdio.h>
#include<stdlib.h>
#include"math.h"

#include"Arguments.h"
#include"Model.h"
#include"Engine.h"

void welcome(){
	fprintf(stderr,
	" ~~~ Welcome ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); 
	fprintf(stderr, 
	" BNLearner:  The Tool for Learning Bayesian Networks\n");
	fprintf(stderr,
	" ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"); 

}
void goodbye(){
	fprintf(stderr,
	"\n ~~~ Goodbye~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
}



int main(int argc, char **argv){
	welcome();

	Arguments::init(argc, argv);
	
	Model model;

	model.init();
	
	Engine engine;

	engine.init(&model);

	engine.compute_edge_probabilities();

	goodbye();
	return 1;
}
