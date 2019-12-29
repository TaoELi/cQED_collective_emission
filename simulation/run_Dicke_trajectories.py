# First copy data from original path to new path
ORIGIN_PATH=src
NEW_PATH=SE_Dicke

# Under new path, we modify some basic parameters
cp do_many_trajs.sh $NEW_PATH

cd $NEW_PATH

for N in 7 9 11 15 19 25 35 45

do
    # First, copy data from old to new
    CURRENT_SUB_PATH=N_$N\_statics
    mkdir $CURRENT_SUB_PATH
    cp N_$N/* $CURRENT_SUB_PATH
    cp do_many_trajs.sh $CURRENT_SUB_PATH
    cd $CURRENT_SUB_PATH
    pwd 
    # change parameters.hpp file
    sed -i".bak" "s/NTRAJ_MULTITRAJ_Eh = 12800/NTRAJ_MULTITRAJ_Eh = 1/g" parameters.hpp
    sed -i".bak" "s/PI \* 3.0/PI \* 0.3/g" parameters.hpp
    sed -i".bak" "s/5000/500/g" parameters.hpp
    sed -i".bak" "s/#define __MPI_ENVIRONMENT__/\/\/#define __MPI_ENVIRONMENT__/g" parameters.hpp
    #compile
    mpic++ mmst_solver.cpp -o run_mmst -std=c++11 -I ~/local/armadillo-9.600.5/include -L ~/local/armadillo-9.600.5/lib64/ -larmadillo -O3
    wait
    
    for i in {1..1000}
    do
        ./run_mmst
        wait
        cp traj_MultiEh.txt traj_MultiEh_$i.txt
    done
    
    wait
    cd ..
done

cd ..
