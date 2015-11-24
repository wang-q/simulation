#!/usr/bin/perl
use strict;
use warnings;

package IDAM;
use Moose;
use Carp;

use Math::GSL::RNG qw(:all);
use Math::GSL::Randist qw(:all);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

use YAML qw(Dump Load DumpFile LoadFile);

# relative mutation rate
has 'R' => ( is => 'ro', isa => 'Num', default => 1.5 );

# selfing rate
has 'F' => ( is => 'ro', isa => 'Num', default => 0 );

# O, original selfing rate
has 'init_F' => ( is => 'ro', isa => 'Num', default => 0 );

# M, mutation rate normal
has 'mu' => ( is => 'ro', isa => 'Num', default => 0.0001 );

# I, mutation rate indel
has 'imu' => ( is => 'ro', isa => 'Num', default => 0.00015 );

# T, generations to run
has 'runtime' => ( is => 'ro', isa => 'Int' );

# Q, when selfing evolves
has 'selftime' => ( is => 'ro', isa => 'Int', default => 10 );

# diploid popsize
has 'N' => ( is => 'ro', isa => 'Int', default => 100 );

# S, random seed
has 'seed' => ( is => 'ro', isa => 'Int', default => int( rand(10000) ) );

# Random Number Generators
has 'rng' => ( is => 'ro', isa => 'Object', );

sub BUILD {
    my $self = shift;

    # intialize gsl random stuff
    # ¡°Mersenne Twister¡± generator
    my $seed = $self->seed;
    my $rng = Math::GSL::RNG->new( $gsl_rng_mt19937, $seed );
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
    $self->rng->free;
    return;
}

sub run {
    my $self = shift;

    my $pop_size    = $self->N;
    my @default_pop = (0) x $pop_size;

    my $rng = $self->rng;
    my $mu  = $self->mu;
    my $R   = $self->R;

    my $mom = [ [@default_pop] ];
    my $dad = [ [@default_pop] ];

    # initial selfing rate to use during sims (changes at selftime)
    my $sim_F = $self->init_F;

    for my $gen ( 0 .. $pop_size - 1 ) {
        my $parent_size = scalar $mom;

        # p: frequency of indel allele
        my $p = [ (0.0) x $parent_size ];

        # keep track of who mutates and how many
        my $mutations = [ (0) x $pop_size ];
        my $muts = 0;

        # measure frequency and remove all sites that are fixed
        $p = $self->measure_p( $mom, $dad );
        if ( $parent_size > 1 ) {
            $self->remove_fixed( $mom, $dad, $p );
        }

        # create new lists for copying purposes
        my $mom_ = [];
        my $dad_ = [];
        for my $i ( 0 .. $parent_size - 1 ) {
            $mom_->[$i] = [@default_pop];
            $dad_->[$i] = [@default_pop];
        }

        # go through each individual in pop
        for my $i ( 0 .. $pop_size - 1 ) {

            # count if indel is in homo vs hetero and populate mutation vector
            if ( $mom->[0][$i] == $dad->[0][$i] ) {

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

            # if not inbred than both alleles from different parents
            if ( $gen >= $self->selftime ) {
                $sim_F = $self->F;
            }
            if ( gsl_rng_uniform($rng) > $sim_F ) {
                $dude = int( gsl_ran_flat( $rng, 0, $pop_size ) );
            }
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

    my $isindel = $freq->[0] - int( $freq->[0] );
    $isindel = 1 if $isindel > 0;

    print $self->F, "\t", $self->R, "\t", $pi, "\t", scalar @$mom, "\t",
        $isindel, "\n";

    return;
}

sub measure_p {
    my $self = shift;
    my $mom  = shift;
    my $dad  = shift;

    my $freq = [];

    my $pop_size = $self->N;
    my $count    = scalar @$mom;
    if ( $count != scalar @$dad ) {
        confess "unequal elements in mom & dad\n";
    }

    for my $i ( 0 .. $count - 1 ) {
        my $sum_p = 0;
        for my $n ( 0 .. $pop_size - 1 ) {
            $sum_p += $mom->[$n] + $dad->[$n];
        }
        $freq->[$i] = $sum_p / ( 2 * $pop_size );
    }

    return $freq;
}

sub remove_fixed {
    my $self = shift;
    my $mom  = shift;
    my $dad  = shift;
    my $freq = shift;

    my $count = scalar @$mom;
    if ( $count != scalar @$dad ) {
        confess "unequal elements in mom & dad\n";
    }
    if ( $count != scalar @$freq ) {
        confess "unequal elements in parents & freqs\n";
    }

    # return indel state to 0's
    if ( $freq->[0] == 1.0 ) {
        for my $i ( 0 .. $self->N - 1 ) {
            $mom->[0][$i] = 0;
            $dad->[0][$i] = 0;
        }
    }

    # find empty vectors and remove from list
    my @remove;
    for my $i ( 0 .. $count - 1 ) {
        if ( $freq->[$i] == 0.0 or $freq->[$i] == 1.0 ) {
            push @remove, $i;
        }
    }

    # remove elements backward
    for my $i ( reverse @remove ) {
        splice @$mom, $i, 1;
        splice @$mom, $i, 1;
        splice @$mom, $i, 1;
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

    my $count = scalar @$p2;
    if ( $count != scalar @$p1 ) {
        confess "unequal elements when copying parents\n";
    }

    for my $i ( 0 .. $count - 1 ) {
        $p1->[$i][$idx] = $p2->[$i][$random_idx];
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

    my $rng         = $self->rng;
    my $mu          = $self->mu;
    my $pop_size    = $self->N;
    my @default_pop = (0) x $pop_size;

    # mutate first at the indel position
    # 2 * N * imu mutants coming in to the indel position
    my $imutants = gsl_ran_poisson( $rng, $self->imu * $pop_size * 2 );

    # add imutants mutations
    for my $i ( 0 .. $imutants - 1 ) {
        my $guy = int( gsl_ran_flat( $rng, 0, $pop_size ) );
        my $new_mutant;

        if ( min( @{ $mom->[0] } ) > 0 and min( @{ $dad->[0] } ) > 1 ) {
            if ( min( @{ $mom->[0] } ) < min( @{ $dad->[0] } ) ) {
                $new_mutant = min( @{ $mom->[0] } ) - 1;
            }
            else {
                $new_mutant = min( @{ $dad->[0] } ) - 1;
            }
        }
        else {
            if ( max( @{ $mom->[0] } ) > max( @{ $dad->[0] } ) ) {
                $new_mutant = max( @{ $mom->[0] } ) + 1;
            }
            else {
                $new_mutant = max( @{ $dad->[0] } ) + 1;
            }
        }

        if ( gsl_rng_uniform($rng) < 0.5 ) {
            $mom->[0][$guy] = $new_mutant;
        }
        else {
            $dad->[0][$guy] = $new_mutant;
        }
    }

    # makes parental lists bigger by mutants No. of mutants for SNPs
    for ( 1 .. $mutants ) {
        splice @$mom, 1, 0, [@default_pop];
        splice @$dad, 1, 0, [@default_pop];
    }

    # puts all mutations in their place
    my $idx = 1;
    for my $i ( 0 .. $pop_size - 1 ) {
        for my $j ( 0 .. $mutant_list->[$i] - 1 ) {
            if ( gsl_rng_uniform($rng) <= 0.5 ) {
                $mom->[$idx][$i] = 1;
            }
            else {
                $dad->[$idx][$i] = 1;
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
