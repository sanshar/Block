/* 
 * File:   ripdm_mpi.h
 * Author: roemelt
 *
 * Created on August 16, 2013, 4:34 PM
 */



#ifndef NEVPT2_MPI_H
#define	NEVPT2_MPI_H

//==============================================================================
// A set if functions that facilitates the parallelization of the BlockNEVPT2 
// module
//==============================================================================
int mpi_world_size();
int mpi_rank();
void mpi_barrier();
void mpi_divide_loop(int &start_loc, int &stop_loc, int start_global, int stop_global);
void mpi_debug_break();
void mpi_message(char *msg);
void PALDivideLoop(int &start_loc, int &stop_loc, int start_global, int stop_global);


//==============================================================================
//A class that allows us to synchronize a loop that has a different length for 
//different processes
//All processes undergo the maximum length loop. This means that some might have
//to undergo some dummy iterations
//==============================================================================
class LoopSynchronizer {
private:
    int MaxDimGlobal;
    int LocalDim;
    bool HaveGlobalDim;
public:
    //the constructor and destructor
    LoopSynchronizer(){
        MaxDimGlobal = -1;
        LocalDim = -1;
        HaveGlobalDim = false;
    }
    LoopSynchronizer(int dl){
        LocalDim = dl;
        HaveGlobalDim = false;
    }
    ~LoopSynchronizer(){};
    
    //Get the global (maximum) dimension. This function implies the calculation
    //of the global dimension
    int GetGlobalDim(){
        if (!HaveGlobalDim) DetermineGlobalDim();
        return MaxDimGlobal;
    }
    
    //determine the (maximum) global dimension
    void DetermineGlobalDim();
    
    //get the local dimension
    int GetLocalDim() const {return LocalDim;}
    
    //set the local dimension
    void SetLocalDim(int ld){LocalDim = ld;}
    
    //tell if we have determined the global dim yet
    bool GetHaveGlobalDim() const {return HaveGlobalDim;}
    
    //tell if the current iteration is a dummy iteration
    bool DummyIter(int iter){return(iter>=LocalDim);}
    
    //reset the Synchronizer
    void Reset(){
        MaxDimGlobal = -1;
        LocalDim = -1;
        HaveGlobalDim = false;
    }
    
};

#endif	/* RIPDM_MPI_H */

