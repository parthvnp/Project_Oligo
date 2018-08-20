#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TurboFold_object.h"
#include "Alignment.h"
#include "../src/rna_library.h"
#include "../Rsample/Rsample.h"
#include <string>
#include "../src/phmm/phmm_aln.h"
#include "../src/phmm/structure/structure_object.h"
#include "../src/phmm/utils/xmath/matrix/matrix.h"
#include "../src/score.h"
#ifdef COMPILE_SMP
    #include "../src/phmm/utils/mutex/ansi_mutex.h"
    #include "../src/phmm/utils/ansi_thread/ansi_thread.h"
#endif
//#include "TurboFold_thread.h"

#ifdef _WINDOWS
    #include <float.h>
    #define isnan _isnan
#endif
//#define debug
//#define SUM


TurboFold::TurboFold(const char _fasta_fp[])
{
    is_RSample_mode = false;
    err_code = 0;
    progress=NULL;
    multiple_sequences=NULL;
    multiple_alignment=NULL;

    char*fasta_fp = new char[strlen(_fasta_fp)+1];
    strcpy(fasta_fp, _fasta_fp);
    vector<t_structure*> &new_sequences = *t_structure::read_multi_seq(fasta_fp);
    const int count = new_sequences.size();
    if(count==0) {
        setError(RNA_LIB_RNA_CONSTRUCTOR_ERROR, "Need at least 1 sequence to predict structure for.");
        return;
    }
    sequences.resize(count);
    folders.resize(count);
    saves.resize(count);
    
    if (readThermo()!=0) return; //error message already set.

    for(unsigned int i_str = 0; i_str < new_sequences.size(); i_str++)
    {
        sequences[i_str] = new_sequences[i_str]; //new t_structure(new_sequences[i_str]);
        RNA* new_folder = new RNA(&(new_sequences[i_str]->nucs[1]), SEQUENCE_STRING, thermo);
        if(new_folder->GetErrorCode() != 0)
        {
            setError(RNA_LIB_RNA_CONSTRUCTOR_ERROR, new_folder->GetFullErrorMessage());
            return;
        }
        folders[i_str]=new_folder;
        //delete(new_sequences[i_str]); // delete the t_structure
    } 

    allocate_extrinsic_information();

    initialize_multiple_sequences();

    turbofold_threads_mutex = NULL;

    if(initialize_alignment_information() != 0)
    {
        err_code = ALIGNMENT_INFO_COMPUTATION_ERROR;
        return;
    }

     err_code = ALIGNMENT_INFO_COMPUTATION_ERROR;
}

TurboFold::TurboFold(vector<string>* _sequences, vector<string>* _saves, string _aln_save)
{
    is_RSample_mode = false;
    err_code = 0;
    turbofold_threads_mutex = NULL;
    progress=NULL;
    multiple_sequences=NULL;
    multiple_alignment=NULL;
    
    vector<string> &seqFiles = *_sequences;
    sequence_names.resize(seqFiles.size());
    for(unsigned int i=0;i<seqFiles.size();i++)
        sequence_names[i] = seqFiles[i];

    const int count = seqFiles.size();
    sequences.resize(count);
    folders.resize(count);
    refolding_job_queue.resize(count);
    saves.resize(count);

    // Alignment output name.
    aln_save_name = _aln_save;
     // data tables shared by all RNA objects
    if (readThermo()!=0) return; //error message already set.
    for(unsigned int i_seq = 0; i_seq < seqFiles.size(); i_seq++)
    {
        string &seqFile = seqFiles[i_seq];
        sequences[i_seq]=new t_structure(seqFile.c_str());
        folders[i_seq]=new RNA(seqFile.c_str(), FILE_SEQ, thermo); 
        if (_saves!=NULL) saves[i_seq]=copy_cstr((*_saves)[i_seq].c_str());
    }

    this->allocate_extrinsic_information();

    this->initialize_multiple_sequences();

    if(this->initialize_alignment_information() != 0)
    {
        this->err_code = ALIGNMENT_INFO_COMPUTATION_ERROR;
        return;
    }
}

int TurboFold::setupRsample(vector<string> *shapeFiles, RsampleData *rsdata, int num_samples, int rsample_seed, double Cparam, double Offset) {
    is_RSample_mode = true;
    this->rsdata = rsdata;
    this->Cparam = Cparam;
    this->Offset = Offset;
    this->rsample_numsamples = num_samples;
    this->rsample_seed = rsample_seed;
    if (err_code==0) // only read restraints if there is no error.
        err_code=read_shape_restraints(*shapeFiles);
    else
        setError(err_code, sfmt("Skipped reading restraints due to prior error: %d", err_code)); // add an error message
    return err_code;
}

void TurboFold::allocate_extrinsic_information()
{
    
    const int count = sequences.size();
    // Allocate extrinsic info.
    basepairing_extrinsic_info_list.resize(count);
    upstream_list.resize(count);
    downstream_list.resize(count);
    unpairing_list.resize(count);
    match_score_list.resize(count);
    
    for(unsigned int i_seq = 0; i_seq < count; i_seq++)
    {
        const int size = sequences[i_seq]->numofbases + 1; // allocate 1 more because numofbases doesn't include the element at position 0.
        basepairing_extrinsic_info_list[i_seq]=new t_matrix(size, size, true);
        upstream_list[i_seq].resize(size, 0);
        downstream_list[i_seq].resize(size, 0);
        unpairing_list[i_seq].resize(size, 0);
        for(unsigned int i_seq2 = i_seq + 1; i_seq2 < sequences.size(); i_seq2++)
            match_score_list[i_seq].push_back(new t_matrix(size, sequences[i_seq2]->numofbases+1, false) );
    }
}

int TurboFold::read_shape_restraints(const vector<string> &shape_files)
{
    const int count = sequences.size();
    shape_reactivities.resize(count);
    
    for(int i=0; i<count; i++) {
        if (shape_files[i].empty())
            shape_reactivities[i]=NULL;
        else {
            shape_reactivities[i]=new vector<double>(sequences[i]->numofbases); // this returns a 0-indexed vector (i.e. base 1 is at position 0)
            int err = ReadRestraints(*shape_reactivities[i], shape_files[i].c_str());
            if (err != 0) {
                cerr << "ReadRestraints Error: " << err << endl;
                cerr << "File: " << shape_files[i] << endl;
                cerr << "Message: " << (string(RNA::GetErrorMessage(err))+"File: \""+shape_files[i]+"\".") << endl;
                return setError(RNALIB_READSHAPE_ERROR,string(RNA::GetErrorMessage(err))+"File: \""+shape_files[i]+"\".");
            }
        }
    }
    return 0; //success
}

TurboFold::~TurboFold()
{
    if (thermo!=NULL) delete thermo;

    const int n_seq = sequences.size();
    for(int i_seq = 0; i_seq < n_seq; i_seq++) {
        if (saves[i_seq]!=NULL) delete[] saves[i_seq];
        delete folders[i_seq];
        delete(basepairing_extrinsic_info_list[i_seq]);
        for(unsigned int i_seq2 = i_seq+1; i_seq2 < n_seq; i_seq2++)
            delete match_score_list[i_seq][i_seq2-i_seq-1];
        free(similarities[i_seq]);
    } 
    free(similarities);

    for(unsigned int i_seq1 = 0; i_seq1 < n_seq; i_seq1++)
    {
        for(unsigned int i_seq2 = i_seq1+1; i_seq2 < n_seq; i_seq2++)
        {
            if(i_seq1 != i_seq2)
            {
                t_phmm_aln* phmm_aln = new t_phmm_aln(sequences[i_seq1], sequences[i_seq2]);

                // Free alignment mapping probabilities.
                for(int i = 1; i <= sequences[i_seq1]->numofbases; i++)
                {
                    free(this->aln_mapping_probs[i_seq1][i_seq2][i]);
                    free(this->aln_probs[i_seq1][i_seq2][i]);
                } // i loop
                free(this->aln_mapping_probs[i_seq1][i_seq2]);
                free(this->aln_probs[i_seq1][i_seq2]);

                // Free alignment envelope result.
                phmm_aln->free_aln_env_result(this->aln_env_results[i_seq1][i_seq2]);
                phmm_aln->free_aln_env_result(this->aln_env_results[i_seq2][i_seq1]);
                delete(phmm_aln);
            }
            else
            {
            }
        } // i_seq2 loop.

        free(aln_env_results[i_seq1]);
        aln_mapping_probs[i_seq1] += i_seq1;
        aln_probs[i_seq1] += i_seq1;
        free(this->aln_mapping_probs[i_seq1]);
        free(this->aln_probs[i_seq1]);
    } // i_seq1 loop.
    

    // sequences are used in memory free'ing. Must make sure it is free'ed at the very last step.
    for(int i_seq = 0; i_seq < n_seq; i_seq++)
    {
        if (sequences[i_seq]!=NULL) delete sequences[i_seq];
    }

    if (is_RSample_mode) {
        for(int i=0; i<n_seq; i++) 
            if (shape_reactivities[i]!=NULL)
                delete shape_reactivities[i]; //delete the vector<double>
    }

    free(aln_env_results);
    free(aln_mapping_probs);
    free(aln_probs);
    #if COMPILE_SMP
    if (turbofold_threads_mutex!=NULL) delete turbofold_threads_mutex;
    #endif
    if (multiple_sequences!=NULL) delete multiple_sequences;
    if (multiple_alignment!=NULL) delete multiple_alignment;
}


//Set the maximum pairing distance for all sequences.
int TurboFold::SetMaxPairingDistance(int distance) {
    int error;

    for(unsigned int i_seq = 0; i_seq < sequences.size(); i_seq++) {

        error = folders[i_seq]->ForceMaximumPairingDistance(distance);
        if (error!=0) {
            err_code = RNALIB_SETDISTANCE_ERROR;
            return RNALIB_SETDISTANCE_ERROR;
        }
    }
    return 0;
}

#ifdef COMPILE_SMP
inline THREAD_HANDLE threadID() { return t_ansi_thread::current_thread_handle(); }
#else
inline unsigned long int threadID() { return 0; }
#endif

/////////////////////////////////////////////////////
// fold - The entrypoint to the TurboFold operation.
/////////////////////////////////////////////////////
int TurboFold::fold(double gamma, int num_iterations, int num_parallel_threads, string AlnFormat, int numColumns)
{
#ifdef COMPILE_SMP
    // Allocate the mutex.
    this->turbofold_threads_mutex = new t_ansi_mutex();
    refolding_threads.resize(min(sequences.size(),num_parallel_threads)-1);
#else
    this->turbofold_threads_mutex = NULL;
#endif
    
    // Generate the alignment information.
    this->n_iterations = num_iterations;

    // Initialize the loops.
    for(int i_iter = 0; i_iter <= n_iterations; i_iter++)
    {
#if defined debug
        cerr<<"i_iter "<<i_iter<<"\n";
#endif
        //Add a coarse update of progress:
        if (progress!=NULL) {
            progress->update((int)((100.0*((double) i_iter))/((double) n_iterations+1)));
        }

        // Set the extrinsic information for each sequence.
        if(i_iter == 0)
        {
            // Initialize the extrinsic information to 1.0 for all the sequences.
            for(int i_seq = 0; i_seq < sequences.size(); i_seq++)
            {
                for(int i = 1; i <= sequences[i_seq]->numofbases; i++)
                {
                    for(int j = i+1; j <= sequences[i_seq]->numofbases; j++)
                    {
                        folders[i_seq]->SetExtrinsic(i, j, 1.0f);
                    } // j loop
                } // i loop         
            } // i_Seq loop 
        }
        else
        {
            // Compute the extrinsic information using the base pairing probabilities just computed for each sequence.
            if(this->generate_alignment_extrinsic_information(i_iter) != 0)
            {
                return(this->err_code);
            }
            if(this->run_phmm_alignment(true) != 0)
            {
                return(this->err_code);
            }
            if(this->generate_folding_extrinsic_information(i_iter, gamma, is_RSample_mode))
            {
                return(this->err_code);
            }
                                                      
            // Set ExtrinsicInformation(const int i, const int j, double ext_info), for each base pair in each sequence.    
            for(int i_seq = 0; i_seq < sequences.size(); i_seq++)
            {
                for(int i = 1; i <= sequences[i_seq]->numofbases; i++)
                {
#if defined debug
                    cerr << "i_iter "<<i_iter << " i_seq " << i_seq << " i " << i << " singlestranded_extrinsic_info_final " << this->singlestranded_extrinsic_info_list->at(i_seq)->at(i) << "\n";
#endif
                    for(int j = i+1; j <= sequences[i_seq]->numofbases; j++)
                    {
                        folders[i_seq]->SetExtrinsic(i, j, basepairing_extrinsic_info_list[i_seq]->x(i, j));
                    } // j loop
                } // i loop             
            } // i_Seq loop              

        } 
        
        //DEBUG_THREADS: cout << "Main Thread Handle: " << threadID() << endl;
        resetJobQueue(i_iter); // initialize the queue use for sequence refolding.
        startThreads(); // start parallel threads (does nothing in serial mode)
        refoldSequences();  //begin a loop to refold all sequences. in SMP mode, multiple threads all run the same loop.
        endThreads(); // wait for all parallel threads to complete (does nothing in serial mode)

        if (this->err_code!=0)
           return this->err_code;

        // Output final alignment
        if(i_iter == n_iterations)
        {
            if(this->generate_alignment_extrinsic_information(i_iter) != 0)
            {
                return(this->err_code);
            }
            if(this->run_phmm_alignment(true) != 0)
            {
                return(this->err_code);
            }
            if(this->run_multiple_alignment() != 0)
            {
                return(this->err_code);
            } 
        }
    } // i_iter loop.

    //Add a coarse update of progress:
    if (progress!=NULL) {
        progress->update(100);
    }
    
    // Write Alignment file
    if (!aln_save_name.empty()) {
        ofstream outfile(aln_save_name.c_str());
        if (!outfile.good())
            return 2; //failed to open file. TODO: see if there is a better error constant to return.

        bool isClustal = AlnFormat != "Fasta";
        multiple_alignment->WriteALN (outfile, numColumns, isClustal);
        if (!outfile.good())
            return 2; //failed to open file. TODO: see if there is a better error constant to return.        
        outfile.close();
    }
    return 0;
}

// This is called by ansi_thread at its start.
void* TurboFold::thread_start(void*arg) {
    TurboFold *t = (TurboFold *)arg;
    t->refoldSequences();
    return NULL;
}

void TurboFold::refoldSequences() {
    int jobID;
    while (-1 != (jobID = nextRefoldingJob()))  // exits the thread if job == NULL
        refoldSequence(jobID);
    //DEBUG_THREADS: cout << "Thread Completed: " << threadID() << endl;
}

// threads call this to get the next available folding job.
// a return value of -1 indicates no more jobs.
int TurboFold::nextRefoldingJob() {
    int found = -1;
    if (this->err_code != 0) return -1; // DO NOT return any more jobs if there was an error. exit immediately.
    #ifdef COMPILE_SMP
        // prevent other threads from entering this function.
        this->turbofold_threads_mutex->lock_mutex();
    #endif // COMPILE_SMP
        // Check if there is a sequence to start folding.
        for(int i = 0; i < refolding_job_queue.size(); i++) {
            if (refolding_job_queue[i].owner == 0) {
                found = i;
                refolding_job_queue[i].owner = 1; // TODO: make owner the thread ID etc.
                break;
            }
        }

    #ifdef COMPILE_SMP
        // allow other threads to enter now.
        this->turbofold_threads_mutex->release_mutex();
    #endif // COMPILE_SMP

    return found;
}

void TurboFold::startThreads() {
    #ifdef COMPILE_SMP
    for(int i = 0; i <refolding_threads.size(); i++) {
        refolding_threads[i] = new t_ansi_thread(&thread_start, this);
        refolding_threads[i]->run_thread();
        //DEBUG_THREADS: cout << "Started thread: " << refolding_threads[i]->get_handle() << endl;
    }
    #endif
}

void TurboFold::endThreads() {
    #ifdef COMPILE_SMP
    for(int i = 0; i <refolding_threads.size(); i++) {
        refolding_threads[i]->wait_thread();
        delete refolding_threads[i];
        refolding_threads[i]=NULL;
    }
    #endif
}

void TurboFold::resetJobQueue(int iter) {
    for(int i = 0; i < refolding_job_queue.size(); i++) {
        refolding_job &job = refolding_job_queue[i];
        job.owner = 0;
        job.iter = iter;
        job.seq = i;
    }
}

void TurboFold::refoldSequence(int jobID)
{
    //DEBUG_THREADS: cout << "refoldSequence " << jobID << " by thread " << threadID() << endl;
    refolding_job &job = refolding_job_queue[jobID];
    const int i_seq = job.seq, i_iter=job.iter;
    RNA &rna = *folders[i_seq];
    const int seed = this->rsample_seed + i_seq + i_iter * sequences.size();

    // Recalculate the sequence using either PartitionFunction or Rsample.
    // On the last iteration, save to PFS if the user has specified one.

    const char* saveFile = (i_iter == n_iterations) ? saves[i_seq] : NULL;
    int ret = 0;
    if(is_RSample_mode && shape_reactivities[i_seq]!=NULL) {
        //DEBUG_THREADS: cout << "Rsample i_seq=" << i_seq << " seed=" << seed << " save: " << (saveFile==NULL?"NULL":saveFile) << endl;
        //DEBUG_SHAPE: WriteRestraints(*(shape_reactivities[i_seq]), sfmt("RampleSHAPE_seq%i_iter%i.txt", i_seq, i_iter));
        ret = rna.Rsample(*(shape_reactivities[i_seq]), *rsdata, seed, saveFile, Cparam, Offset, rsample_numsamples);
    } else {
        //DEBUG_THREADS: cout << "PartitionFunction i_seq=" << i_seq << " save: " << (saveFile==NULL?"NULL":saveFile) << endl;
        //DEBUG: rna.GetStructure()->WriteSHAPE(sfmt("shapedat_seq%i_iter%i.txt", i_seq, i_iter).c_str());
        ret = rna.PartitionFunction(saveFile);
    }
    if(ret != 0)
        this->err_code = RNALIB_PARTITIONFUNCTION_ERROR;
}

int TurboFold::generate_alignment_extrinsic_information(int i_iter)
{
        for(int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
        {
            fill(upstream_list[i_seq1].begin(), upstream_list[i_seq1].end(), 0); 
            fill(downstream_list[i_seq1].begin(), downstream_list[i_seq1].end(), 0);
            fill(unpairing_list[i_seq1].begin(), unpairing_list[i_seq1].end(), 0);

            for(int i_seq2 = i_seq1 + 1; i_seq2 < sequences.size(); i_seq2++)
            {
                for(int i = 0; i <= sequences[i_seq1]->numofbases ; i++)
                {
                    for(int j = 0; j <= sequences[i_seq2]->numofbases ; j++)
                    {
                        match_score_list[i_seq1][i_seq2-i_seq1-1]->x(i,j) = 1;
                    }
                }
            }
        }

        // Main computation loop.
        // upstream_; downstream_; unpairing_list; for each nt for each sequence
        for(int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
        {
            RNA* RNA1 = folders[i_seq1];
            for(int i = 1; i <= sequences[i_seq1]->numofbases; i++)
            {
                // upstream
                for(int j = 0; j < i; j++)
                {
                    upstream_list[i_seq1][i] += RNA1->GetPairProbability(j,i);
                }
                // downstream
                for(int j = i + 1; j <= sequences[i_seq1]->numofbases; j++)
                {
                    downstream_list[i_seq1][i] += RNA1->GetPairProbability(i,j);
                }
                // unpairing
                unpairing_list[i_seq1][i] += 1 - upstream_list[i_seq1][i] - downstream_list[i_seq1][i];
            }
        }

        // match_score_list 
        for(int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
        {
            for(int i_seq2 = i_seq1+1; i_seq2 < sequences.size(); i_seq2++)
            {
                for(int i = 1; i <= sequences[i_seq1]->numofbases; i++)
                {
                    for(int k = 1; k <= sequences[i_seq2]->numofbases; k++)
                    {
                        double temp_upstream_info = sqrt(upstream_list[i_seq1][i] * upstream_list[i_seq2][k]);
                        double temp_downstream_info = sqrt(downstream_list[i_seq1][i] * downstream_list[i_seq2][k]);
                        double temp_unpairing_info = sqrt(unpairing_list[i_seq1][i] * unpairing_list[i_seq2][k]);

                        match_score_list[i_seq1][i_seq2-i_seq1-1]->x(i,k) *= (temp_upstream_info + temp_downstream_info) * 1.0 + temp_unpairing_info * 0.8 + 0.5;
                    }
                }
            }
        }

        return(0);

}    


int TurboFold::generate_folding_extrinsic_information(int i_iter, const double gamma, bool is_RSample_mode)
{

    // Reset the extrinsic information.
    for(unsigned int i_seq = 0; i_seq < sequences.size(); i_seq++)
    {
           
        for(int i = 1; i <= sequences[i_seq]->numofbases; i++)
        {
                        
            for(int j = i+1; j <= sequences[i_seq]->numofbases; j++)
            {
                basepairing_extrinsic_info_list[i_seq]->x(i,j) = 0.0f;
            } // j loop
        } // i loop 
    } // i_seq loop


    // Main computation loop.
    for(unsigned int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
    {
        RNA* RNA1 = folders[i_seq1];

        for(unsigned int i_seq2 = i_seq1+1; i_seq2 < sequences.size(); i_seq2++)
        {       
            RNA* RNA2 = folders[i_seq2];
     
            for(int i = 1; i <= sequences[i_seq1]->numofbases; i++)
            {
                for(int j = i+1; j <= sequences[i_seq1]->numofbases; j++)
                {
                    int low_k = max(1, this->aln_env_results[i_seq1][i_seq2]->low_limits[i]);
                    int high_k = this->aln_env_results[i_seq1][i_seq2]->high_limits[i];
                    // min_k starts from 0, because low_k starts from 1, we need this to make the adjustment to get the coresponding aln_mapping_probs[i_seq1][i_seq2][i][k-min_k+1]
                    int min_k = this->aln_env_results[i_seq1][i_seq2]->low_limits[i];

                    for(int k = low_k; 
                            k <= high_k;
                            k++)
                    {
                        int low_l = max(k+1, this->aln_env_results[i_seq1][i_seq2]->low_limits[j]);
                        int high_l = this->aln_env_results[i_seq1][i_seq2]->high_limits[j];
                        // min_l starts from 0 
                        int min_l = this->aln_env_results[i_seq1][i_seq2]->low_limits[j];

                        for(int l = low_l; l <= high_l; l++)
                        {
                            double seq1_seq2_mapping_probability = this->aln_mapping_probs[i_seq1][i_seq2][i][k-min_k+1] * this->aln_mapping_probs[i_seq1][i_seq2][j][l-min_l+1];
                            double seq1_seq2_seq_similarity_weight = (1.0f - this->similarities[i_seq1][i_seq2]);
                            

                            // Note: shape_reactivities[N] is the shape reactivity data for sequence N
                            if(is_RSample_mode&&shape_reactivities[i_seq1]!=NULL)
                            {
                                // To make the extrinsic_info from the sequence with SHAPE to other sequences not impact by seq_similarity_weight, 
                                // seq1_seq2_mapping_probability for first line and not second.
                                // where i_seq1 is the sequence index and i is the position in shape_sequence_indexes (which is now just the same as i_seq)
                                basepairing_extrinsic_info_list[i_seq2]->x(k,l) += seq1_seq2_mapping_probability * RNA1->GetPairProbability(i,j) * (sequences.size() - 1);
                                basepairing_extrinsic_info_list[i_seq1]->x(i,j) += seq1_seq2_seq_similarity_weight * seq1_seq2_mapping_probability * RNA2->GetPairProbability(k,l);                                
                            } else  {
                                basepairing_extrinsic_info_list[i_seq2]->x(k,l) += seq1_seq2_seq_similarity_weight * seq1_seq2_mapping_probability * RNA1->GetPairProbability(i,j);
                                basepairing_extrinsic_info_list[i_seq1]->x(i,j) += seq1_seq2_seq_similarity_weight * seq1_seq2_mapping_probability * RNA2->GetPairProbability(k,l);
                            }

                            basepairing_extrinsic_info_list[i_seq2]->x(k,l) += seq1_seq2_seq_similarity_weight * seq1_seq2_mapping_probability * RNA1->GetPairProbability(i,j);
                            basepairing_extrinsic_info_list[i_seq1]->x(i,j) += seq1_seq2_seq_similarity_weight * seq1_seq2_mapping_probability * RNA2->GetPairProbability(k,l);
                            
                            if(RNA1->GetErrorCode() != 0)
                                return setError(RNALIB_GETPAIRPROBABILITY_ERROR, sfmt("Problem getting pairing probability for (%d, %d) in sequence %d\n", i, j, i_seq1));

                            if(RNA2->GetErrorCode() != 0)
                                return setError(RNALIB_GETPAIRPROBABILITY_ERROR, sfmt("Problem getting pairing probability for (%d, %d) in sequence %d\n", k, l, i_seq2));
                        } // l loop
                    } // k loop
                } // j loop
            } // i loop.
        } // i_seq2 loop.
    } // i_seq1 loop.


        // Post process. For normalization and powerizing.
        for(unsigned int i_seq = 0; i_seq < sequences.size(); i_seq++)
        {
            // Now normalize everything. The case where maximum is 0 is handled inside normalize_by_max function.
            basepairing_extrinsic_info_list[i_seq]->normalize_by_max();
            basepairing_extrinsic_info_list[i_seq]->powerize_each_element(gamma);
        } // i_seq loop.

    return(0);

}

int TurboFold::allocate_phmm()
{
    this->aln_env_results = (t_aln_env_result***)malloc(sizeof(t_aln_env_result**) * (sequences.size() + 1));
    this->aln_mapping_probs = (double****)malloc(sizeof(double***) * (sequences.size() + 1));
    this->aln_probs = (double****)malloc(sizeof(double***) * (sequences.size() + 1));

    this->similarities = (double**)malloc(sizeof(double*) * (sequences.size() + 2));
    
    for(unsigned int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
    {

        this->similarities[i_seq1] = (double*)malloc(sizeof(double) * (sequences.size() + 2));

        aln_env_results[i_seq1] = (t_aln_env_result**)malloc((sequences.size() + 2) * sizeof(t_aln_env_result*));                         
        
        this->aln_mapping_probs[i_seq1] = (double***)malloc((sequences.size() + 2) * sizeof(double**));
        this->aln_probs[i_seq1] = (double***)malloc((sequences.size() + 2) * sizeof(double**));

        // Move the pointers.
        this->aln_mapping_probs[i_seq1] -= i_seq1;
        this->aln_probs[i_seq1] -= i_seq1;

        for(unsigned int i_seq2 = i_seq1+1; i_seq2 < sequences.size(); i_seq2++)
        {
            this->similarities[i_seq1][i_seq2] = 0.0f;
                            
            if(i_seq1 != i_seq2)
            {
                this->aln_mapping_probs[i_seq1][i_seq2] = (double**)malloc(sizeof(double*) * (sequences[i_seq1]->numofbases + 2));  
                this->aln_probs[i_seq1][i_seq2] = (double**)malloc(sizeof(double*) * (sequences[i_seq1]->numofbases + 2));  

                for(int i = 1; i <= sequences[i_seq1]->numofbases; i++)
                {
                    this->aln_mapping_probs[i_seq1][i_seq2][i] = NULL;
                    this->aln_probs[i_seq1][i_seq2][i] = NULL;                    
                }
            }
            else
            {
                aln_env_results[i_seq1][i_seq2] = NULL;
            }
        }
    }
    return(0);
}
    


int TurboFold::run_phmm_alignment(bool using_prior)
{

    for(unsigned int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
    {

        for(unsigned int i_seq2 = i_seq1+1; i_seq2 < sequences.size(); i_seq2++)
        {
            if(i_seq1 != i_seq2)
            {
                t_phmm_aln* phmm_aln = new t_phmm_aln(sequences[i_seq1], sequences[i_seq2]);

                if(using_prior)
                {
                    phmm_aln->set_match_priors(match_score_list[i_seq1][i_seq2-i_seq1-1] );
                } 

                t_pp_result* cur_pp_results = phmm_aln->compute_posterior_probs();

                t_aln_env_result* cur_aln_env_result = phmm_aln->compute_alignment_envelope(PROB_ALN_ENV, cur_pp_results, cur_pp_results->fam_threshold, 7);

                // Copy the ML similarity.
                this->similarities[i_seq1][i_seq2] = cur_pp_results->ml_similarity;

                if(using_prior)
                {
                    free(aln_env_results[i_seq1][i_seq2]->low_limits);
                    free(aln_env_results[i_seq1][i_seq2]->high_limits);                           
                    free(aln_env_results[i_seq1][i_seq2]);
                    free(aln_env_results[i_seq2][i_seq1]->low_limits);
                    free(aln_env_results[i_seq2][i_seq1]->high_limits);                          
                    free(aln_env_results[i_seq2][i_seq1]);
                }

                aln_env_results[i_seq1][i_seq2] = cur_aln_env_result;
                aln_env_results[i_seq2][i_seq1] = (t_aln_env_result*)malloc(sizeof(t_aln_env_result));
                aln_env_results[i_seq2][i_seq1]->low_limits = (int*)malloc(sizeof(int) * (sequences[i_seq2]->numofbases + 2));
                aln_env_results[i_seq2][i_seq1]->high_limits = (int*)malloc(sizeof(int) * (sequences[i_seq2]->numofbases + 2));
                // Initialize loop limits.
                for(int i = 0; i <= sequences[i_seq2]->numofbases + 1; i++)
                {
                    aln_env_results[i_seq2][i_seq1]->low_limits[i] = sequences[i_seq1]->numofbases;
                    aln_env_results[i_seq2][i_seq1]->high_limits[i] = 1;
                }

                for(int k = 1; k <= sequences[i_seq2]->numofbases; k++)               
                {
                    for(int i = 1; i <= sequences[i_seq1]->numofbases; i++)
                    {
                        if(aln_env_results[i_seq2][i_seq1]->low_limits[k] > i && cur_aln_env_result->high_limits[i] >= k)
                        {
                            aln_env_results[i_seq2][i_seq1]->low_limits[k] = i;
                            break;
                        }
                    }
                    
                    for(int i = sequences[i_seq1]->numofbases; i >= 1; i--)
                    {
                        if(aln_env_results[i_seq2][i_seq1]->high_limits[k] < i && cur_aln_env_result->low_limits[i] <= k)
                        {
                            aln_env_results[i_seq2][i_seq1]->high_limits[k] = i;
                            break;
                        }
                    }
                }

                aln_env_results[i_seq2][i_seq1]->low_limits[0] = aln_env_results[i_seq2][i_seq1]->low_limits[1];
                aln_env_results[i_seq2][i_seq1]->high_limits[0] = aln_env_results[i_seq2][i_seq1]->high_limits[1];

                for(int i = 0; i <= sequences[i_seq2]->numofbases; i++)
                {
                    if(aln_env_results[i_seq2][i_seq1]->low_limits[i] == 1)
                    {
                        aln_env_results[i_seq2][i_seq1]->low_limits[i] = 0;
                    }
                }

                for(int i = 1; i <= sequences[i_seq1]->numofbases; i++)
                {

                    int min_k = cur_aln_env_result->low_limits[i];
                    int max_k = cur_aln_env_result->high_limits[i];
                        
                    if(this->aln_mapping_probs[i_seq1][i_seq2][i] != NULL)
                    {
                        free(this->aln_mapping_probs[i_seq1][i_seq2][i]);
                    }
                    
                    if(this->aln_probs[i_seq1][i_seq2][i] != NULL)
                    {
                        free(this->aln_probs[i_seq1][i_seq2][i]);
                    }                    
                            
                    this->aln_mapping_probs[i_seq1][i_seq2][i] = (double*)malloc(sizeof(double) * (max_k - min_k + 2));
                    this->aln_probs[i_seq1][i_seq2][i] = (double*)malloc(sizeof(double) * (max_k - min_k + 2));
                    // 
                       
                    for(int k = min_k; k <= max_k; k++)
                    {
                        // Mapping probability is coincidence probability of the nucleotides.
                        double aln_prob = exp(cur_pp_results->aln_probs[i][k]);
    
                        double ins1_prob = exp(cur_pp_results->ins1_probs[i][k]);
                        double ins2_prob = exp(cur_pp_results->ins2_probs[i][k]);
                        //double mapping_probability;
                        double mapping_probability = aln_prob + ins1_prob + ins2_prob;

                        this->aln_mapping_probs[i_seq1][i_seq2][i][k - min_k + 1] = mapping_probability;
                        this->aln_probs[i_seq1][i_seq2][i][k - min_k + 1] = aln_prob;

#if defined debug
                        cerr<<"i_seq1 "<<i_seq1<<" i_seq2 "<<i_seq2<<" i k "<<i<< " "<<k<< " aln_prob "<<aln_prob<<"\n";
#endif
    
                    } // k loop
                } // i loop

                // Free current pp result.
                phmm_aln->free_pp_result(cur_pp_results);
                delete(phmm_aln);
            }
        } // i_seq2 loop.
    } // i_seq1 loop.

#if defined debug
        for(int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
        { 
            for(int i_seq2 = i_seq1 + 1; i_seq2 < sequences.size(); i_seq2++)
            {

                for(int i = 1; i <= sequences[i_seq1]->numofbases; i++) {
                    cerr <<"i_seq1 "<<i_seq1<<" i_seq2 "<<i_seq2<<" 1 i indel_prob "<<i<<" "<<this->indel_prob_list->at(i_seq1)->at(i_seq2-i_seq1-1)->at(i)<<"\n";
                    
                }

                for(int k = 1; k <= sequences[i_seq2]->numofbases; k++) {
                    cerr <<"i_seq1 "<<i_seq1<<" i_seq2 "<<i_seq2<<" 2 k indel_prob "<<k<<" "<<this->indel_prob_list->at(i_seq1)->at(i_seq2-i_seq1-1)->at(k+sequences[i_seq1]->numofbases)<<"\n";
                    
                }
            }
        }
#endif
    return(0);
}

int TurboFold::run_multiple_alignment()
{
    SafeVector<SafeVector<SparseMatrix *> > sparseMatrices( this -> GetNumberSequences(), SafeVector<SparseMatrix *>(this->GetNumberSequences(), NULL));
    SafeVector<SafeVector<float> > distances (this -> GetNumberSequences(), SafeVector<float> (this -> GetNumberSequences(), 0));
    ProbabilisticModel model;
    for(unsigned int i_seq1 = 0; i_seq1 < sequences.size(); i_seq1++)
    {

        for(unsigned int i_seq2 = i_seq1+1; i_seq2 < sequences.size(); i_seq2++)
        {

            if(i_seq1 != i_seq2)
            {
                SafeVector<float> *posteriorPtr = new SafeVector<float>((sequences[i_seq1]->numofbases+1) * (sequences[i_seq2]->numofbases+1),0); assert (posteriorPtr);
                SafeVector<float> &posterior = *posteriorPtr;
            
                int pointer=0+(sequences[i_seq2]->numofbases+1);
                SafeVector<float>::iterator ptr = posterior.begin()+(sequences[i_seq2]->numofbases+1);

                for(int i = 1; i <= sequences[i_seq1]->numofbases; i++)
                {
                    int min_k = max(1,aln_env_results[i_seq1][i_seq2]->low_limits[i]);
                    int max_k = aln_env_results[i_seq1][i_seq2]->high_limits[i];
                    pointer+=min_k;
                    ptr+=min_k;
                    for(int k = min_k; k <= max_k; k++)
                    {
                        *ptr=static_cast<float>(this->aln_probs[i_seq1][i_seq2][i][k - min_k + 1]);
                        ptr++;
                        pointer++;
                    }
                    ptr+=sequences[i_seq2]->numofbases-max_k;
                    pointer+=sequences[i_seq2]->numofbases-max_k;
                } // i loop

                if(i_seq1!=sequences.size()-1){ 
                    sparseMatrices[i_seq1][i_seq2]= new SparseMatrix(sequences[i_seq1]->numofbases,sequences[i_seq2]->numofbases, *posteriorPtr);
                    sparseMatrices[i_seq2][i_seq1]=NULL;
                }

                pair<SafeVector<char> *, float> pair_alignment = model.ComputeAlignment(sequences[i_seq1]->numofbases,sequences[i_seq2]->numofbases,
                                                                                          *posteriorPtr);
                 
                // compute "expected accuracy" distance for evolutionary tree computation
                float distance = pair_alignment.second / min (sequences[i_seq1]->numofbases,sequences[i_seq2]->numofbases);
                distances[i_seq1][i_seq2] = distances[i_seq2][i_seq1] = distance;
                delete pair_alignment.first;
                delete posteriorPtr;

            }
        } // i_seq2 loop.
    } // i_seq1 loop.


    for (int r = 0; r<numConsistencyReps; r++ ) {
        SafeVector<SafeVector<SparseMatrix *> > newSparseMatrices = MultiConsistencyTransform(this->multiple_sequences,sparseMatrices);

        // now replace the old posterior matrices
        for (int i = 0; i < this -> GetNumberSequences(); i++){
            for (int j = 0; j < this -> GetNumberSequences(); j++){
                delete sparseMatrices[i][j];
                sparseMatrices[i][j] = newSparseMatrices[i][j];
            }
        }
    }
    
    this->multiple_sequences->SaveOrdering();

    this->multiple_alignment=NULL;
    TreeNode *tree = TreeNode::ComputeTree(distances);
  
    // make the final alignment
    this->multiple_alignment = ComputeFinalAlignment(tree, this->multiple_sequences, sparseMatrices, model);
  
    delete tree;

    // delete sparse matrices
    for (int a = 0; a < this->GetNumberSequences()-1; a++){
        for (int b = a+1; b < this->GetNumberSequences(); b++){
            delete sparseMatrices[a][b];
            delete sparseMatrices[b][a];
        }
    }
    return(0);
}

int TurboFold::initialize_alignment_information()
{
    this->allocate_phmm();
    this->run_phmm_alignment();
    return(0);
}

void TurboFold::initialize_multiple_sequences()
{
    this->multiple_sequences=new MultiSequence();
    for(int i = 0 ; i < this->GetNumberSequences(); i++){
    
    SafeVector<char>* data= new SafeVector<char>;
    data->push_back('@');
    for (int cnt=1; cnt<=sequences[i]->numofbases;cnt++){
        data->push_back(toupper(sequences[i]->nucs[cnt]));
    }
    Sequence* temp_seq=new Sequence(data,string(sequences[i]->ctlabel),sequences[i]->numofbases,i,i);
    this->multiple_sequences->AddSequence(temp_seq);
    }

}

int TurboFold::ProbKnot(const int i_seq, const int n_iterations, const int min_helix_length)
{
    if(i_seq > this->GetNumberSequences())
    {
        this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;

    }
    else
    {

    int ret = this->folders[i_seq-1]->ProbKnot(n_iterations, min_helix_length);

    if(ret != 0)
    {
        this->err_code = RNALIB_PROBKNOT_ERROR;
    }
        else
        {
            this->err_code = 0;
        }
    }

    return(this->err_code);
}

int TurboFold::PredictProbablePairs(const int i_seq, const float probability)
{
    if(i_seq > this->GetNumberSequences())
    {
        this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;

    }
    else
    {
        int ret = this->folders[i_seq-1]->PredictProbablePairs(probability);

        if(ret != 0)
        {
            this->err_code = RNALIB_THRESHOLDING_ERROR;
        }
        else
        {
        this->err_code = 0;
        }
    }
    return(this->err_code);

}

int TurboFold::MaximizeExpectedAccuracy(const int i_seq, const double maxPercent, const int maxStructures, const int window, const double gamma)
{
    if(i_seq > this->GetNumberSequences())
    {
    this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;
    
    }
    else
    {
        int ret = this->folders[i_seq-1]->MaximizeExpectedAccuracy(maxPercent, maxStructures, window, gamma);

        if(ret != 0)
        {
            this->err_code = RNALIB_MEA_ERROR;
        }
        else
        {
            this->err_code = 0;
        }
    }
        return(this->err_code);

}

int TurboFold::GetPair(const int i_seq, const int i, const int structurenumber)
{
    if(i_seq > this->GetNumberSequences())
    {
            this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;
    }
    else
    {
        int ret = this->folders[i_seq-1]->GetPair(i, structurenumber);

        if(ret != 0)
        {
            this->err_code = RNALIB_GETPAIR_ERROR;
        }
        else
        {
            this->err_code = 0;
        }
    }

    return(this->err_code);
}

int TurboFold::WriteCt(const int i_seq, const char fp[])
{
    if(i_seq > this->GetNumberSequences())
    {
        this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;
    }
    else
    {
        int ret = this->folders[i_seq-1]->WriteCt(fp);

        if(ret != 0)
        {
            this->err_code = RNALIB_WRITECT_ERROR;
        }
        else
        {
            this->err_code = 0;
        }
    }       

    return(this->err_code);

}

// int TurboFold::OutputEmpty()
// {
//
//     int existed = 0;
//     for(unsigned int i_seq = 0;i_seq < sequences.size();++i_seq){
//         size_t position = this->sequence_names[i_seq].rfind(".");
//         if (position!=string::npos) this->sequence_names[i_seq].erase(pos);
//         this->sequence_names[i_seq].append(".ct");
//         RNA* temp_rna = new RNA(this->sequence_names[i_seq].c_str(),1);
//         structure* existed_structure = temp_rna->GetStructure();
//         for(int i = 1; i<=existed_structure->GetSequenceLength();i++) {
//             if(existed_structure->GetPair(i,1)>i)existed++;
//         }
//         delete temp_rna;
//     }
//     return existed;
// }

double TurboFold::GetPairProbability(const int i_seq, const int i, const int j)
{
    if(i_seq > this->GetNumberSequences())
    {
        this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;
        return(0.0f);
    }
    else
    {
        double pair_prob = this->folders[i_seq-1]->GetPairProbability(i,j);
        int ret = this->folders[i_seq-1]->GetErrorCode();

        if(ret != 0)
        {
            this->err_code = RNALIB_GETPAIRPROBABILITY_ERROR;
            return(0.0f);
        }
        else
        {
            this->err_code = 0;
            return(pair_prob);
        }
    }
}

int TurboFold::ReadSHAPE(const int i_seq, const char fp[], const double par1, const double par2)
{
    if(i_seq > this->GetNumberSequences())
    {
        this->err_code = ISEQ_ARGUMENT_OVERFLOW_ERROR;
    }
    else
    {
        //Note that 0.0 and 0.0 are hard-wired because these are single-stranded terms that do not apply using current practices.
        int ret = this->folders[i_seq-1]->ReadSHAPE(fp, par1, par2, 0.0, 0.0);
        if(ret != 0)
            return setError(RNALIB_READSHAPE_ERROR, this->folders[i_seq-1]->GetFullErrorMessage());
        else
            this->err_code = 0;
    }
    return(this->err_code);
}

int TurboFold::GetNumberSequences()
{
    return(sequences.size());
}

int TurboFold::readThermo() {
   thermo = new Thermodynamics();
   int ret = thermo->ReadThermodynamic();
   if (ret != 0)
     return setError(CONSTRUCTOR_ERROR, RNA::GetErrorMessage(ret));
   return 0;
}

int TurboFold::SetTemperature(const double temp)
{
    int ret = thermo->SetTemperature(temp);
    if(ret!=0)
        return setError(RNALIB_SETTEMPERATURE_ERROR);
    return 0;
}


int TurboFold::GetErrorCode()
{
    return(this->err_code);
}

const char* TurboFold::GetErrorMessage(const int err_code)
{
    if (err_code>=0&&err_code<N_TF_ERRORS)
        return err_strings[err_code];
    return "Unknown Error Code";
}

string TurboFold::GetErrorString(const int err_code)
{
    return GetErrorMessage(err_code);
}

//Provide a TProgressDialog for following calculation progress.
//A TProgressDialog class has a public function void update(int percent) that indicates the progress of a long calculation.
void TurboFold::SetProgress(ProgressHandler& Progress) {
    progress = &Progress;
}

//Provide a means to stop using a TProgressDialog.
//StopProgress tells the RNA class to no longer follow progress.  This should be called if the TProgressDialog is deleted, so that this class does not make reference to it.
void TurboFold::StopProgress() {
    progress=NULL;
}

// Set the internal error code and optionally set a detailed error message.
// Note that if the error code has ALREADY been set, calling this function will NOT modify it.
//    however the new error message (if any) will be appended to any existing error message.
// Setting overwrite to true causes the new error code and error message to completely overwrite any existing values.
// \return The current error code. (Usually this will be the one passed in, unless the error code was already set before this call.
int TurboFold::setError(const int code, string errorDetails, const bool overwrite) {
    if (err_code==0||overwrite) err_code = code;
    if (!errorDetails.empty()) {
        if (lastErrorDetails.empty() || overwrite)
            lastErrorDetails = errorDetails;
        else
            lastErrorDetails = lastErrorDetails + "\n" + errorDetails;
    }
    return err_code;
}