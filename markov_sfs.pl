#!/usr/bin/perl
use strict;
use warnings;

package MarkovSFS;
use Moose;
with 'MooseX::Getopt';

use Math::GSL::RNG qw(:all);
use Math::GSL::Randist qw(:all);
use PDL;
use PDL::IO::FITS;

use YAML qw(Dump Load DumpFile LoadFile);

has 'pop_size' => (
    is            => 'ro',
    isa           => 'Int',
    default       => 100,
    traits        => ['Getopt'],
    cmd_aliases   => [qw{ n N }],
    documentation => "diploid popsize",
);
has 'runtime' => (
    is            => 'ro',
    isa           => 'Int',
    traits        => ['Getopt'],
    cmd_aliases   => [qw{ t T }],
    documentation => "generations to run",
);
has 'max_allele' => (
    is            => 'ro',
    isa           => 'Int',
    traits        => ['Getopt'],
    cmd_aliases   => [qw{ m M }],
    documentation => "max allele number",
);
has 'mu' => (
    is            => 'ro',
    isa           => 'Num',
    default       => 0.05,
    traits        => ['Getopt'],
    cmd_aliases   => [qw{ u U }],
    documentation => "mutation rate normal",
);
has 'epsilon' => (
    is            => 'ro',
    isa           => 'Num',
    default       => 0.001,
    traits        => ['Getopt'],
    cmd_aliases   => [qw{ e E }],
    documentation => "gene conversion rate",
);
has 'output' => (
    is            => 'rw',
    isa           => 'Str',
    traits        => ['Getopt'],
    cmd_aliases   => [qw{ o O }],
    documentation => "output dir name",
);

# suffix for
has 'suffix' => ( is => 'rw', isa => 'Str', default => "001" );

# random seed
has 'seed' => (
    is      => 'ro',
    isa     => 'Int',
    default => int( rand(10000) ),
    traits  => ['Getopt']
);

# Random Number Generators
has '_rng' => ( is => 'ro', isa => 'Object', );

# internal storage
has '_mom'   => ( is => 'ro', isa => 'Object', );
has '_dad'   => ( is => 'ro', isa => 'Object', );
has '_freq'  => ( is => 'ro', isa => 'ArrayRef[Int]', default => sub { [] } );
has '_used'  => ( is => 'ro', isa => 'ArrayRef[Int]', default => sub { [] } );
has '_empty' => ( is => 'ro', isa => 'ArrayRef[Int]', default => sub { [] } );

# fixed and all allele count
has '_all_allele'   => ( is => 'ro', isa => 'Int', default => 0 );
has '_drift_allele' => ( is => 'ro', isa => 'Int', default => 0 );
has '_fixed_allele' => ( is => 'ro', isa => 'Int', default => 0 );
has '_lost_allele'  => ( is => 'ro', isa => 'Int', default => 0 );

# current generation
has '_gen' => ( is => 'ro', isa => 'Int', default => 0 );

# allele dynamic
has '_dynamic' =>
    ( is => 'ro', isa => 'HashRef[Ref]', default => sub { {} } );

sub BUILD {
    my $self = shift;

    # intialize gsl random stuff
    # ¡°Mersenne Twister¡± generator
    my $seed = $self->seed;
    my $rng  = gsl_rng_alloc($gsl_rng_mt19937);
    gsl_rng_set( $rng, $seed );
    $self->{_rng} = $rng;

    # running time to get to pseudo-equilibrium
    # default is 8N
    unless ( $self->runtime ) {
        $self->{runtime} = 8 * $self->pop_size;
    }

    # max drifting alleles
    unless ( $self->max_allele ) {
        $self->{max_allele} = 2 * $self->mu * $self->pop_size * 20;
    }

    # default output filename
    unless ( $self->output ) {
        $self->{output}
            = "Freq[N"
            . $self->pop_size . "][T"
            . $self->runtime . "][M"
            . $self->max_allele . "][Mu"
            . $self->mu . "][Eps"
            . $self->epsilon . "]";
    }

    # use matrix
    # rows represent alleles, columns represent individual chromosomes
    #           0   1   2   3   4    5   6   7   8   9
    #        mom0   1   2   3   4 dad0   1   2   3   4
    # allele0   1   1   0   0   0    1   1   0   0   0
    # allele1   0   0   0   0   1    1   0   0   0   0
    # allele2   1   0   0   0   0    0   0   0   0   0
    # allele3   1   0   0   0   0    0   0   1   0   1
    $self->{_mom} = byte( zeroes( $self->pop_size, $self->max_allele ) );
    $self->{_dad} = byte( zeroes( $self->pop_size, $self->max_allele ) );

    $self->{_freq} = [ (0) x $self->max_allele ];
    $self->{_empty} = [ 0 .. $self->max_allele - 1 ];

    return;
}

sub DEMOLISH {
    my $self = shift;
    gsl_rng_free( $self->{_rng} );
    return;
}

sub run {
    my $self = shift;

    for my $gen ( 0 .. $self->runtime - 1 ) {
        $self->{_gen} = $gen;

        # measure frequency and remove all alleles that are fixed
        $self->measure_freq;
        $self->reset_fixed;

        $self->report;

        #----------------------------#
        # gene convert
        #----------------------------#
        $self->gene_convert;

        #----------------------------#
        # mutate
        #----------------------------#
        $self->mutate;

        #----------------------------#
        # reproduce
        #----------------------------#
        $self->reproduce;
    }

    # last calc and report
    $self->measure_freq;
    $self->reset_fixed;

    $self->report("END");
    $self->record;

    return;
}

sub measure_freq {
    my $self = shift;
    my $freq = ( $self->_mom->daverage + $self->_dad->daverage ) / 2;
    $self->{_freq} = [ $freq->list ];
    return;
}

sub reset_fixed {
    my $self = shift;

    my $freq = $self->_freq;

    # find fixed allele in used
    # In back order, so the indexes don't change
    my $used = $self->_used;
    for my $i ( reverse( 0 .. scalar @$used - 1 ) ) {
        my $allele = $used->[$i];
        if ( $freq->[$allele] == 1 ) {
            $freq->[$allele] = 0;
            splice @$used, $i, 1;
            unshift @{ $self->_empty }, $allele;
            $self->{_mom}->slice(",$allele") .= 0;
            $self->{_dad}->slice(",$allele") .= 0;
            $self->{_fixed_allele}++;
        }
        elsif ( $freq->[$allele] == 0 ) {
            splice @$used, $i, 1;
            unshift @{ $self->_empty }, $allele;
            $self->{_lost_allele}++;
        }
    }

    $self->{_used}         = $used;
    $self->{_freq}         = $freq;
    $self->{_drift_allele} = scalar @$used;

    return;
}

sub mu_poisson {
    my $self = shift;
    my $muts = $self->mu * $self->pop_size * 2;
    return gsl_ran_poisson( $self->_rng, $muts );
}

sub epsilon_poisson {
    my $self  = shift;
    my $count = shift;
    my $cons  = $self->epsilon * $count;
    return gsl_ran_poisson( $self->_rng, abs $cons );
}

sub random_pick {
    my $self = shift;
    my $max  = shift;
    return int gsl_ran_flat( $self->_rng, 0, $max );
}

sub random_pick_n {
    my $self = shift;
    my $max  = shift;
    my $n    = shift;
    return map { int gsl_ran_flat( $self->_rng, 0, $max ) } ( 1 .. $n );
}

sub new_allele_idx {
    my $self = shift;

    my $idx = shift @{ $self->_empty };
    if ( defined $idx ) {
        push @{ $self->_used }, $idx;
        return $idx;
    }
    else {
        return;
    }
}

sub mutate {
    my $self = shift;

    my $mu       = $self->mu;
    my $pop_size = $self->pop_size;

    # 2 * N * mu mutants coming in to all chromosomes
    my $muts  = $self->mu_poisson;
    my $empty = scalar @{ $self->_empty };

    return if $empty <= 0;

    if ( $empty <= $muts ) {
        $muts = $empty;
    }
    $self->{_all_allele} += $muts;
    $self->{_dynamic}->{ $self->_gen }->{newMuts} = $muts;
    print "# Add $muts mutations\n";

    my $max = 2 * $self->pop_size - 1;
    for ( 1 .. $muts ) {
        my $dude_idx   = $self->random_pick($max);
        my $allele_idx = $self->new_allele_idx;
        if ( $dude_idx < $pop_size ) {
            $self->{_mom}->slice("$dude_idx,$allele_idx") .= 1;
        }
        else {
            $dude_idx = $dude_idx - $pop_size;
            $self->{_dad}->slice("$dude_idx,$allele_idx") .= 1;
        }
    }

    return;
}

sub gene_convert {
    my $self = shift;

    my $epsilon  = $self->epsilon;
    my $pop_size = $self->pop_size;

    my $hetero = which( $self->_mom != $self->_dad );
    my $count  = $hetero->nelem;

    my $cons = $self->epsilon_poisson($count);

    $self->{_dynamic}->{ $self->_gen }->{newCons}
        = $epsilon >= 0 ? $cons : -$cons;
    print "# Convert $cons ",
        $epsilon >= 0 ? "wild types to mutants\n" : "mutants to wild types\n";

    my %seen;
    for ( 1 .. $cons ) {
        my $idx = $self->random_pick($count);
        redo if exists $seen{$idx};
        $seen{$idx}++;
        my $real_idx = $hetero->index($idx);

        my ( $wild_type, $mutant ) = ( "_mom", "_dad" );
        if ( $self->{_mom}->flat->index($real_idx) == 1 ) {
            ( $wild_type, $mutant ) = ( "_dad", "_mom" );
        }

        if ( $epsilon >= 0 ) {
            $self->{$wild_type}->flat->index($real_idx) .= 1;
        }
        else {
            $self->{$mutant}->flat->index($real_idx) .= 0;
        }
    }

    return;
}

sub reproduce {
    my $self = shift;

    my $pop_size = $self->pop_size;

    my $matrix = $self->_mom->append( $self->_dad );

    # go through each individual in pop
    my $max = 2 * $self->pop_size - 1;

    # pick chromosomes from random parents to mom and dad
    my @moms = $self->random_pick_n( $max, $pop_size );
    my @dads = $self->random_pick_n( $max, $pop_size );

    $self->{_mom} = $matrix->dice_axis( 0, \@moms )->sever;
    $self->{_dad} = $matrix->dice_axis( 0, \@dads )->sever;

    return;
}

sub report {
    my $self = shift;
    my $gen = shift || $self->_gen;

    printf "Gen:[%6s]\tAll:%8d\tLost:%8d\tFixed:%6d\tDrift:%6d\n", $gen,
        $self->_all_allele, $self->_lost_allele, $self->_fixed_allele,
        $self->_drift_allele;

    $self->{_dynamic}->{$gen} = {
        All   => $self->_all_allele,
        Lost  => $self->_lost_allele,
        Fixed => $self->_fixed_allele,
        Drift => $self->_drift_allele,
    };

    return;
}

sub record {
    my $self = shift;

    my $outdir;
    while (1) {
        $outdir = $self->output . $self->suffix;
        if ( -d $outdir ) {
            $self->{suffix}++;
        }
        else {
            mkdir $outdir;
            last;
        }
    }

    my @drifting;
    for my $allele ( @{ $self->_used } ) {
        push @drifting, $self->_freq->[$allele];
    }
    DumpFile( "$outdir/freq.yml", \@drifting );
    DumpFile(
        "$outdir/opt.yml",
        {   N   => $self->pop_size,
            T   => $self->runtime,
            M   => $self->max_allele,
            Mu  => $self->mu,
            Eps => $self->epsilon,
        }
    );
    DumpFile( "$outdir/dynamic.yml", $self->_dynamic );

    $self->_mom->wfits("$outdir/mom.fits");
    $self->_dad->wfits("$outdir/dad.fits");

    my $matrix = $self->_mom->append( $self->_dad );
    my $order  = pdl $self->_freq;
    $order = $order->qsorti;
    $matrix = $matrix->dice_axis( 1, $order );
    $matrix->wfits("$outdir/matrix.fits");

    eval { require Image::Magick; };
    if ( !$@ ) {
        for (qw{matrix.fits mom.fits dad.fits}) {
            my $model = Image::Magick->new;
            $model->ReadImage("$outdir/$_");
            $model->Write("$outdir/$_.png");
        }
    }

    return;
}

package main;

my $markov_sfs = MarkovSFS->new_with_options;
$markov_sfs->run;
