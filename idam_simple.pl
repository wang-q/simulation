#!/usr/bin/perl
use strict;
use warnings;

package IDAM;
use Moose;
use Carp;

use Math::GSL::RNG qw(:all);
use Math::GSL::Randist qw(:all);
use PDL;
use PDL::Ufunc;

use YAML qw(Dump Load DumpFile LoadFile);

# relative mutation rate
has 'R' => ( is => 'ro', isa => 'Num', default => 1.5 );

# M, mutation rate normal
has 'mu' => ( is => 'ro', isa => 'Num', default => 10 );

# I, mutation rate indel
has 'imu' => ( is => 'ro', isa => 'Num', default => 15 );

# T, generations to run
has 'runtime' => ( is => 'ro', isa => 'Int' );

# diploid popsize
has 'N' => ( is => 'ro', isa => 'Int', default => 10 );

# S, random seed
has 'seed' => ( is => 'ro', isa => 'Int', default => int( rand(10000) ) );

# Random Number Generators
has 'rng' => ( is => 'ro', isa => 'Object', );

sub BUILD {
    my $self = shift;

    # intialize gsl random stuff
    # ¡°Mersenne Twister¡± generator
    my $seed = $self->seed;
    my $rng  = gsl_rng_alloc($gsl_rng_mt19937);
    gsl_rng_set( $rng, $seed );
    $self->{rng} = $rng;

    # running time to get to pseudo-equilibrium
    # default is 8N
    unless ( $self->runtime ) {
        $self->{runtime} = 8 * $self->N;
    }

    return;
}

sub DEMOLISH {
    my $self = shift;
    gsl_rng_free( $self->{rng} );
    return;
}

sub run {
    my $self = shift;

    # list of int vectors
    #   -- each vector is a site with state for each individual
    # vectors for current ($mom & $dad) and next (M_ and P_) gens.
    my $pop_size     = $self->N;
    my $default_site = zeroes($pop_size);

    my $rng = $self->rng;
    my $mu  = $self->mu;
    my $R   = $self->R;

    my $mom = [ pdl $default_site ];
    my $dad = [ pdl $default_site ];

    for my $gen ( 0 .. $self->runtime - 1 ) {
        print "On generation [$gen]\n";
        
        # p: frequency of indel allele, double
        my $site_number = scalar @$mom;
        my $p = [ (0.0) x $site_number ];
        print "At the beginning, there are [$site_number] sites\n";

        # measure frequency and remove all sites that are fixed
        $p = $self->measure_p( $mom, $dad );
        #print "Mutation freq is @$p\n";
        print "Indel freq is $p->[0]\n";
        
        if ( $site_number > 1 ) {
            $self->remove_fixed( $mom, $dad, $p );
        }
        $site_number = scalar @$mom;
        print "After remove fixed sites, there are [$site_number] sites\n";

        # keep track of who mutates and how many, int
        my $mutations = [ (0) x $pop_size ];
        my $muts = 0;

        # create new lists for copying purposes
        my $mom_ = [];
        my $dad_ = [];
        for my $i ( 0 .. $site_number - 1 ) {
            $mom_->[$i] = pdl $default_site;
            $dad_->[$i] = pdl $default_site;
        }

        # go through each individual in pop
        for my $i ( 0 .. $pop_size - 1 ) {

            # count if indel is in homo vs hetero and populate mutation vector
            if ( $mom->[0]->index($i) == $dad->[0]->index($i) ) {

                # 2 * N * mu mutations per generation
                $mutations->[$i] = gsl_ran_poisson( $rng, $mu )
                    + gsl_ran_poisson( $rng, $mu );
                $muts += $mutations->[$i];
            }
            else {
                $mutations->[$i] = gsl_ran_poisson( $rng, $mu * $R )
                    + gsl_ran_poisson( $rng, $mu * $R );
                $muts += $mutations->[$i];
            }

            # copy a chromosome from random parent
            my $dude = int gsl_ran_flat( $rng, 0, $pop_size );
            if ( gsl_rng_uniform($rng) <= 0.5 ) {
                $self->copy_csome( $i, $dude, $mom_, $mom );
            }
            else {
                $self->copy_csome( $i, $dude, $mom_, $dad );
            }

            $dude = int gsl_ran_flat( $rng, 0, $pop_size );
            if ( gsl_rng_uniform($rng) <= 0.5 ) {
                $self->copy_csome( $i, $dude, $dad_, $mom );
            }
            else {
                $self->copy_csome( $i, $dude, $dad_, $dad );
            }
        }

        # add new sites based on number of mutations
        $self->mutate( $mom_, $dad_, $mutations, $muts );

        # copy lists back to originals;
        $mom = $mom_;
        $dad = $dad_;
    }

    # freq: frequency of indel allele
    # final measurement of frequency
    my $freq = $self->measure_p( $mom, $dad );
    if ( @$mom > 1 ) {
        $self->remove_fixed( $mom, $dad, $freq );
    }

    # calculate pi at nonindel sites and output
    my $pi = 0;
    for ( 1 .. scalar @$freq - 1 ) {
        $pi += 2.0 * $freq->[$_] * ( 1 - $freq->[$_] );
    }
    
    my $site_number = scalar @$mom;

    my $isindel = $freq->[0] - int( $freq->[0] );
    $isindel = 1 if $isindel > 0;

    print "R $self->R\t pi $pi \t site number $site_number \t isindel $isindel \n";

    return;
}

sub measure_p {
    my $self = shift;
    my $mom  = shift;
    my $dad  = shift;

    my $freq = [];

    my $pop_size    = $self->N;
    my $site_number = scalar @$mom;
    if ( $site_number != scalar @$dad ) {
        confess "unequal elements in mom & dad\n";
    }

    for my $site ( 0 .. $site_number - 1 ) {
        my $sum_p = 0;
        for my $i ( 0 .. $pop_size - 1 ) {
            $sum_p += $mom->[$site]->index($i) + $dad->[$site]->index($i);
        }
        $freq->[$site] = $sum_p / ( 2 * $pop_size );
    }

    return $freq;
}

sub remove_fixed {
    my $self = shift;
    my $mom  = shift;
    my $dad  = shift;
    my $freq = shift;

    my $site_number = scalar @$mom;
    if ( $site_number != scalar @$dad ) {
        confess "unequal elements in mom & dad\n";
    }
    if ( $site_number != scalar @$freq ) {
        confess "unequal elements in parents & freqs\n";
    }

    # return indel state to 0's
    if ( $freq->[0] == 1.0 ) {
        for my $i ( 0 .. $self->N - 1 ) {
            $mom->[0]->index($i) = 0;
            $dad->[0]->index($i) = 0;
        }
    }

    # find empty vectors and remove from list
    # the first is indel, so we skip it
    my @remove;
    for my $site ( 1 .. $site_number - 1 ) {
        if ( $freq->[$site] == 0.0 or $freq->[$site] == 1.0 ) {
            push @remove, $site;
        }
    }

    # remove elements backward
    for my $site ( reverse @remove ) {
        splice @$mom,  $site, 1;
        splice @$dad,  $site, 1;
        splice @$freq, $site, 1;
    }

    return;
}

# makes p1 same as p2
sub copy_csome {
    my $self       = shift;
    my $idx        = shift;
    my $random_idx = shift;
    my $p1         = shift;
    my $p2         = shift;

    my $site_number = scalar @$p2;
    if ( $site_number != scalar @$p1 ) {
        confess "unequal elements when copying parents\n";
    }

    for my $site ( 0 .. $site_number - 1 ) {
        $p1->[$site]->index($idx) .= $p2->[$site]->index($random_idx);
    }

    return;
}

# mutate (  int mutants, vector<int> *mutant_list,  parent *mom, parent *dad )
sub mutate {
    my $self        = shift;
    my $mom         = shift;
    my $dad         = shift;
    my $mutant_list = shift;
    my $mutants     = shift;

    my $rng          = $self->rng;
    my $mu           = $self->mu;
    my $pop_size     = $self->N;
    my $default_site = zeroes($pop_size);

    # mutate first at the indel position
    # 2 * N * imu mutants coming in to the indel position
    my $imutants = gsl_ran_poisson( $rng, $self->imu * $pop_size * 2 );

    # add imutants mutations
    for my $i ( 0 .. $imutants - 1 ) {
        my $guy = int( gsl_ran_flat( $rng, 0, $pop_size ) );
        my $new_mutant;

        if ( $mom->[0]->min > 0 and $dad->[0]->min > 1 ) {
            if ( $mom->[0]->min < $dad->[0]->min ) {
                $new_mutant = $mom->[0]->min - 1;
            }
            else {
                $new_mutant = $dad->[0]->min - 1;
            }
        }
        else {
            if ( $mom->[0]->max > $dad->[0]->max ) {
                $new_mutant = $mom->[0]->max + 1;
            }
            else {
                $new_mutant = $dad->[0]->max + 1;
            }
        }

        if ( gsl_rng_uniform($rng) < 0.5 ) {
            $mom->[0]->index($guy) .= $new_mutant;
        }
        else {
            $dad->[0]->index($guy) .= $new_mutant;
        }
    }

    # makes parental lists bigger by mutants No. of mutants for SNPs
    for ( 1 .. $mutants ) {
        splice @$mom, 1, 0, pdl $default_site;
        splice @$dad, 1, 0, pdl $default_site;
    }

    # puts all mutations in their place
    my $idx = 1;
    for my $i ( 0 .. $pop_size - 1 ) {
        for my $j ( 0 .. $mutant_list->[$i] - 1 ) {
            if ( gsl_rng_uniform($rng) <= 0.5 ) {
                $mom->[$idx]->index($i) .= 1;
            }
            else {
                $dad->[$idx]->index($i) .= 1;
            }
            $idx++;
        }
    }

    return;
}

package main;
use YAML qw(Dump Load DumpFile LoadFile);

my $idam = IDAM->new;
$idam->run;

#print Dump $idam;
#print $idam->rng->get, "\n" for 1 .. 10;
