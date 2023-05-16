Readme for BNLearner


I. Contents:
----------------

This package includes the tool BNLearner which is used to compute the exact posterior probabilities of various structural features in Bayesian Networks.

For the detailed algorithm, please refer to the following paper:

"Structure Learning in Bayesian Networks of Moderate Size by Efficient Sampling", Ru He, Jin Tian and Huaiqing Wu. 






II. Acknowledge
-----------------------

For the code reuse as well as the fair performance comparison, this tool modifies and/or reuses some code from the two tools: REBEL tool and POSTER tool.


REBEL tool implements the algorithms proposed in the following two papers:

"Exact Bayesian Structure Discovery in Bayesian Networks", Mikko Koivisto and Kismat Sood. Journal of Machine 
Learning Research, 5:549-573, 2004.

 "Advances in Exact Bayesian Structure Discovery in Bayesian Networks", Mikko Koivisto. In Proceedings of the 
Conference on Uncertainty in Artificial Intelligence (UAI), 2006. 



POSTER tool implements the algorithms proposed in the following paper:

"Computing Posterior Probabilities of Structural Features in Bayesian Networks", Jin Tian and Ru He.  In Proceedings of 
the Conference on Uncertainty in Artificial Intelligence (UAI), 2009.





III. Investigated data cases:
------------------

Note that the BNLearner tool also includes all the data cases investigated the paper. Please see these data cases under the sub-directory "cases/".





IV. Compile 
------------------

g++ -Wall -O2 -o BNLearner main.cc




V. Run
----------------------


<1. Example Runs for the DDS:

./BNLearner -m 5 -d cases/tic_tac_toe.txt  -u 958 -opt 9   -so 20000          > results/tic_tac_toe_dds_m5_ord20000.txt
./BNLearner -m 5 -d cases/letter_samp100.txt  -u 100 -opt 9  -so 20000    > results/letter_samp100_dds_m5_ord20000.txt  
./BNLearner -m 5 -d cases/letter_samp500.txt  -u 500 -opt 9  -so 20000    > results/letter_samp500_dds_m5_ord20000.txt 


Explanation of command options:

 -d Data file:
                                            cases/tic_tac_toe.txt
 -m Maximum indegree:
                                                                5
 -u Maximum number of data records read:
                                                              958
 -opt Option:
                                                                9
 -ad ADtree:
                                                                0
 -so Number of sampled orders:
                                                            20000

(-opt 9 is used to compute the posterior probabilities using the DDS.)


Note:
Inside the sub-directory "results/", 
the file "edge_poster" records the computed posterior probability of each edge feature.




<2. Example Runs for the IW-DDS:

./BNLearner -m 5 -d cases/tic_tac_toe.txt  -u 958 -opt 11  -so 30000         > results/tic_tac_toe_iwdds_m5_ord30000.txt
./BNLearner -m 5 -d cases/letter_samp100.txt -u 100 -opt 11  -so 30000   > results/letter_samp100_iwdds_m5_ord30000.txt



Explanation of command options:

-d Data file:
                                           cases/tic_tac_toe.txt
-m Maximum indegree:
                                                               5
-u Maximum number of data records read:
                                                             958
-opt Option:
                                                              11
-ad ADtree:
                                                               0
-so Number of sampled orders:
                                                           30000

(-opt 11 is used to compute the posterior probabilities using the IW-DDS.)


Note:
Inside the sub-directory "results/", 
the file "edge_poster" records the computed posterior probability of each edge feature;
the file "path_poster" records the computed posterior probability of each f_1 (path feature) mentioned in the paper;
the file "limited_leng_path_poster" records the computed posterior probability of each f_2 (limited-length path feature) mentioned in the paper;
the file "combined_two_path_poster" records the computed posterior probability of each f_3 mentioned in the paper;
the file "combined_two_f4_path_poster" records the computed posterior probability of each f_4 mentioned in the paper;
the file "combined_two_f5_path_poster" records the computed posterior probability of each f_5 mentioned in the paper.




<3. Example Runs for the IW-DDS with a huge N_o (i.e., N_o >= 2,000,000)

./BNLearner -m 5 -d cases/InsuranceSub19_samp200.txt  -u 200 -opt 12   -so 10000000 > results/InsuranceSub19_samp200_iwdds_m5_ord10000000.txt

(-opt 12 is used to compute the posterior probabilities using the IW-DDS with a huge N_o, that is, N_o >= 2,000,000.)




<4. Example Runs for the DOS:

./BNLearner -m 5 -d cases/tic_tac_toe.txt  -u 958 -opt 8   -so 20000          > results/tic_tac_toe_dos_m5_ord20000.txt
./BNLearner -m 5 -d cases/letter_samp500.txt  -u 500 -opt 8  -so 20000    > results/letter_samp500_dos_m5_ord20000.txt 

(-opt 8 is used to compute the posterior probabilities using the DOS.)





Note that BNLearner tool also includes the relevant functions in POSTER tool and REBEL tool .


<5. Example Runs for the exact method in the paper "Computing Posterior Probabilities of Structural Features in Bayesian Networks".

./BNLearner -m 5 -d cases/tic_tac_toe.txt  -u 958 -opt 3  > results/tic_tac_toe_mix_m5.txt 
./BNLearner -m 5 -d cases/letter_samp500.txt  -u 500 -opt 3  > results/letter_samp500_mix_m5.txt  

(-opt 3 is used to compute the exact posterior probability of each edge without the bias from the order-modular prior.)




<6. Example Runs for the exact method in REBEL tool:

(-opt 0 is used to compute the exact posterior probability of each edge with the bias from the order-modular prior.)


For K2 score with prior  1 / nchoosek(n-1, |Pa_i|)  for \rho_i(Pa_i)
./BNLearner -m 5 -d cases/tic_tac_toe.txt  -u 958  -opt 0   	-re 1            > results/tic_tac_toe_K2_m5.txt 
./BNLearner -m 5 -d cases/letter_samp500.txt  -u 500 -opt 0        	-re 1 	> results/letter_samp500_K2_m5.txt  


For BDeu score with prior  1 / nchoosek(n-1, |Pa_i|)  for \rho_i(Pa_i)
./BNLearner -m 5 -d cases/tic_tac_toe.txt  -u 958  -opt 0   	-re 2           > results/tic_tac_toe_rebelOld_m5.txt 
./BNLearner -m 5 -d cases/letter_samp500.txt  -u 500 -opt 0        	-re 2 	> results/letter_samp500_rebelOld_m5.txt  


For BDeu score with prior  1  for \rho_i(Pa_i)
./BNLearner -m 5 -d cases/tic_tac_toe.txt  -u 958  -opt 0   	-re 3           > results/tic_tac_toe_rebelNew_m5.txt 
./BNLearner -m 5 -d cases/letter_samp500.txt  -u 500 -opt 0        	-re 3 	> results/letter_samp500_rebelNew_m5.txt  



VI. Contacts
----------------------

Please send bug reports to Ru He by his email hrheru AT gmail DOT com.

