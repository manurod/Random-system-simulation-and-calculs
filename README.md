# Simulation and calcul of Random Systems

This python librarie contain several fonction allowing us to simulate a continus time markov chaine (CMTC) an to compute the 1/) terme C : E[X]=x*+C/N+O(1/N^2).

## List of fonctions available 

- simu_cmtc(taux,liste_transition,X_0,N,final_time,my_file.out,fix) 
the last two arguments are optionals, the argument "my_file.out" will creat a file and will store result in this file by default this argument is "" then the fonction will return an array. The last argument allow us to fix the seed of randoms fonctions.

- file_to_array(file,n) return the array which has been store with simu_cmtc (if we choose to store it)

- fixed_point(taux,liste_transition,n) compute the fix point with the help of the fonction odeint.

- theorique(taux,liste_transitions,n,dim) compute the coefficient C describe at the begining with theoricals formulas
