#include <iostream>
#include <vector>
#include <list>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h> // Random Number Generation
#include <gsl/gsl_randist.h> // Random Number Distributions

using namespace std;

struct params {
    double F, R, mu, imu, initF;
    int N, seed, runtime, selftime;
};

typedef list < vector <int> > parent;

void copy_csome ( unsigned ind, unsigned random_ind, parent *p1, parent * p2);

void remove_fixed(  parent *mom, parent *dad, list<double> *freq, params *pop );

void mutate ( gsl_rng *r, int mutants, vector<int> *mutant_list, struct params *pop, parent *mom, parent *dad );

void measure_p ( struct params *pop, list< double > *freq, parent *mom, parent *dad );

void print_out ( parent *p1, parent *p2, params *pop );

int main(int argc, char *argv[]) {

    int opt=1;
    params pop;

    // get variables
    while ( opt < argc ) {
        switch ( argv[opt][1] ) {
        case 'R': opt++; pop.R = atof(argv[opt++]); break; //relative mutation rate
        case 'F': opt++; pop.F = atof(argv[opt++]); break; //observed data file name
        case 'M': opt++; pop.mu = atof(argv[opt++]); break; //mutation rate normal
        case 'I': opt++; pop.imu = atof(argv[opt++]); break; //mutation rate indel
        case 'T': opt++; pop.runtime = atoi(argv[opt++]); break; //generations to run
        case 'O': opt++; pop.initF = atof( argv[opt++]); break; //original selfing rate 
        case 'Q': opt++; pop.selftime = atoi( argv[opt++]); break; //when selfing evolves
        case 'S': opt++; pop.seed = atoi( argv[opt++]); break; //random seed
        case 'N': opt++; pop.N = atoi(argv[opt++]); break; //diploid popsize
        }
    } 

    if ( argc < 18 ) {
        cerr << "\nIncorrect command line specifications.\nOptions are:\n" 
        << "R: relative mutation rate for indel hets\n"
        << "F: selfing rate\n"
        << "O: original selfing rate\n"
        << "M: mutations rate per locus at SNPs\n"
        << "I: mutations rate per site at the indel\n"
        << "T: generations to run\n"
        << "Q: generations at which selfing evolves\n"
        << "S: random seed\n"
        << "N: diploid popsize\n\n";
        exit(1);
    }

    // intialize gsl random stuff
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, pop.seed);

    // list of int vectors -- each vector is a site with state for each individual
    // vectors for current (M & P) and next (M_ and P_) gens. 
    parent P( 1, vector<int>(pop.N, 0) ); // list of one element
    parent M( P ); // copy P to M

    // initial selfing rate to use during sims (changes at time pop.selftime)
    double simF = pop.initF;

    // run long to get to pseudo-equilibrium
    for ( int gen = 0; gen < pop.runtime; gen++ ) {

        // p: frequency of indel allele
        list <double> p(M.size(), 0.0);

        // keep track of who mutates and how many
        vector <int> mutations (pop.N, 0);
        int muts = 0;

        // measure frequency and remove all sites that are fixed
        measure_p ( &pop, &p, &M, &P );
        if ( M.size()>1 )
            remove_fixed( &M, &P, &p, &pop );

        // create new lists for copying purposes
        parent M_ (M.size(), vector<int>(pop.N,0)); 
        parent P_ (P.size(), vector<int>(pop.N,0));

        // go through each individual in pop
        for ( int i=0; i<pop.N; i++ ) {

            // count if indel is in homo vs. hetero and populate mutation vector
            if ( (*M.begin())[i] == (*P.begin())[i] ) {
                //2*N*mu mutations per generation
                mutations[i] = gsl_ran_poisson(r,pop.mu)+gsl_ran_poisson(r,pop.mu); 
                muts += mutations[i];
            } else {
                mutations[i] = gsl_ran_poisson(r,pop.mu*pop.R)+gsl_ran_poisson(r,pop.mu*pop.R); 
                muts += mutations[i]; 
            } 

            //copy a chromosome from random parent
            int dude = int ( gsl_ran_flat(  r, 0, pop.N  ) );
            if ( gsl_rng_uniform(r)<=0.5 )
                copy_csome( i, dude, &M_, &M );
            else
                copy_csome( i, dude, &M_, &P );

            // if not inbred than both alleles from different parents

            if ( gen>=pop.selftime )
                simF=pop.F;
            if ( gsl_rng_uniform(r)>simF )
                dude = int ( gsl_ran_flat(  r, 0, pop.N  ) );
            if ( gsl_rng_uniform(r)<=0.5 )
                copy_csome( i, dude, &P_, &M );
            else
                copy_csome ( i, dude, &P_, &P );    
        }

        // add new sites based on number of mutations
        mutate( r, muts, &mutations, &pop, &M_, &P_ );

        //copy lists back to originals;
        M=M_; P=P_;     

        //print_out( &M_, &P_, &pop);
    }

    //freq frequency of indel allele
    list <double> freq(M.size(),0.0);

    //final measurement of frequency
    measure_p ( &pop, &freq, &M, &P );
    if ( M.size()>1 )
        remove_fixed( &M, &P, &freq, &pop );

    //calculate pi at nonindel sites and output
    double pi=0;
    list<double>::iterator it1 = freq.begin(); it1++;
    for ( list<double>::iterator it = it1; it != freq.end(); ) {
        //below is output for SFS
        //cout << (*it)*2*pop.N << " ";
        pi += 2.0*(*it)*(1-*it);
        it++;
    }

    double isindel = (*freq.begin()) - int( (*freq.begin()) );
    if ( isindel > 0 )
        isindel = 1;

    cout << pop.F << "\t" << pop.R << "\t" << pi << "\t" << M.size() << "\t"  << isindel << "\n";
}

void mutate ( gsl_rng *r, int mutants, vector<int> *mutant_list, struct params *pop, parent *mom, parent *dad ) {
    //mutate first at the indel position
    int imutants = gsl_ran_poisson(r,pop->imu*pop->N*2); // 2*N*imu mutants coming in to the indel position
    list< vector < int > >::iterator i_m = (*mom).begin(); 
    list< vector < int > >::iterator i_p = (*dad).begin();

    //add imutants mutations
    for ( int i=0; i<imutants; i++ ) {
        int guy = int ( gsl_ran_flat(  r, 0, pop->N  ) );
        int new_mutant;

        if ( *min_element((*i_m).begin(), (*i_m).end()) > 0 && *min_element((*i_p).begin(), (*i_p).end()) > 1 ) {
            if ( *min_element((*i_m).begin(), (*i_m).end()) < *min_element((*i_p).begin(), (*i_p).end()) )
                new_mutant = *min_element((*i_m).begin(), (*i_m).end()) - 1;
            else new_mutant = *min_element((*i_p).begin(), (*i_p).end()) - 1;
        } else {
            if ( *max_element((*i_m).begin(), (*i_m).end()) > *max_element((*i_p).begin(), (*i_p).end()) )
                new_mutant = *max_element((*i_m).begin(), (*i_m).end()) + 1;
            else new_mutant = *max_element((*i_p).begin(), (*i_p).end()) + 1;
        }

        if ( gsl_rng_uniform(r) < 0.5 )
            (*i_m)[guy]=new_mutant;
        else
            (*i_p)[guy]=new_mutant;
    }

    //makes parental lists bigger by mutants # of mutants for SNPS
    list< vector < int > >::iterator itm = (*mom).begin(); 
    list< vector < int > >::iterator itp = (*dad).begin();
    itm++; itp++;

    (*mom).insert(itm,mutants,vector<int> (pop->N,0) );
    (*dad).insert(itp,mutants,vector<int> (pop->N,0) );

    //puts all mutations in their place
    list< vector < int > >::iterator im = (*mom).begin(); 
    list< vector < int > >::iterator ip = (*dad).begin();
    im++; ip++;

    for ( int i=0; i<pop->N; i++ ) {
        for ( int j=0; j<(*mutant_list)[i]; j++ ) {
            if ( gsl_rng_uniform(r)<=0.5 )
                (*im)[i]=1;
            else
                (*ip)[i]=1;
            ip++; im++;
        }
    } 
} 

void copy_csome ( unsigned ind, unsigned random_ind, parent *p1, parent * p2) {
    // makes p1 same as p2
    list< vector < int > >::iterator p2_loc = (*p2).begin();
    for ( list< vector < int > >::iterator p1_loc = (*p1).begin(); p1_loc != (*p1).end(); ) {
        (*p1_loc)[ind]=(*p2_loc)[random_ind];   
        p1_loc++; p2_loc++;
    } 
}

void remove_fixed(  parent *mom, parent *dad, list<double> *freq, params *pop ) {
    // find empty vectors and remove from list
    list< vector < int > >::iterator p_ind, m_ind;
    list<double>::iterator ind;
    p_ind=(*dad).begin(); 
    m_ind=(*mom).begin(); 
    ind=(*freq).begin();    

    //return indel state to 0's
    if ( *ind == 1.0 ) {
        for ( int j=0; j< pop->N; j++ ) {
            (*p_ind)[j]=0; (*m_ind)[j]=0;       
        }
    }
    ind++; p_ind++; m_ind++;

    for ( list<double>::iterator it = ind; it != (*freq).end(); ) {
        if ( *it == 0.0 || *it == 1.0 ) {
            it = (*freq).erase(it); 
            m_ind = (*mom).erase(m_ind); 
            p_ind = (*dad).erase(p_ind); 
        } else {
            it++; m_ind++; p_ind++; 
        }
    } 
}

void measure_p( struct params *pop, list< double > *freq, parent *mom, parent *dad ) {
    // measure frequency at each site
    list< vector < int > >::iterator ip = (*dad).begin();
    list< double >::iterator ifr = (*freq).begin();

    //over sites
    for ( list< vector < int > >::iterator im = (*mom).begin(); im != (*mom).end(); ) {
        int sum_p=0;
        for ( int i=0; i<pop->N; i++ ) {
            sum_p=sum_p+(*im)[i]+(*ip)[i];
        }
        (*ifr)=sum_p/(2.0*pop->N);
        ip++; im++; ifr++;
    }
}

void print_out ( parent *p1, parent *p2, params *pop ) {
    // prints out mom and dad for each locus
    cout << "\n\n";
    for ( int i=0; i<pop->N; i++ ) {
        for ( list< vector < int > >::iterator ip = (*p1).begin(); ip != (*p1).end(); ) {
            cout  << " "<< (*ip)[i];
            ip++;
        }
        cout << endl;
        for ( list< vector < int > >::iterator im = (*p2).begin(); im != (*p2).end(); ) {
            cout  << " " << (*im)[i];
            im++;
        }
        cout << endl;
    } 
}

