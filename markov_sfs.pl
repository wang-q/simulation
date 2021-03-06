#!/usr/bin/perl
use strict;
use warnings;
use autodie;

package MarkovSFS;
use Moo;
use MooX::Options;

use Math::Random::MT::Auto qw(:!auto);
use PDL;
use PDL::IO::FITS;
use GD;
use Path::Tiny;
use Statistics::R;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

option 'pop_size' => (
    is      => 'ro',
    format  => 'i',
    short   => 'n',
    default => 100,
    doc     => "diploid popsize, default is 100",
);
option 'runtime' => (
    is     => 'ro',
    format => 'i',
    short  => 't',
    doc    => "generations to run, default is 8N",
);
option 'max_allele' => (
    is     => 'ro',
    format => 'i',
    short  => 'm',
    doc    => "max allele number, default is 2N",
);
option 'mu' => (
    is      => 'ro',
    format  => 'f',
    short   => 'u',
    default => 0.01,
    doc     => "mutation rate normal, default is 0.01",
);
option 'epsilon' => (
    is      => 'ro',
    format  => 'f',
    short   => 'e',
    default => 0.001,
    doc     => "gene conversion rate, default is 0.001",
);
option 'output' => (
    is     => 'rw',
    format => 's',
    short  => 'o',
    doc    => "output dir name",
);

# suffix for
has 'suffix' => ( is => 'rw', default => "01" );    # Str

# Random Number Generators
has 'rng' => ( is => 'ro', );                       # Object

# internal storage
has 'mom'   => ( is => 'ro', );                           # Object
has 'dad'   => ( is => 'ro', );                           # Object
has 'freq'  => ( is => 'ro', default => sub { [] }, );    # ArrayRef[Int]
has 'used'  => ( is => 'ro', default => sub { [] }, );    # ArrayRef[Int]
has 'empty' => ( is => 'ro', default => sub { [] }, );    # ArrayRef[Int]

# current generation
has 'gen' => ( is => 'ro', default => 0, );               # Int

# allele dynamic
has 'dynamic' => ( is => 'ro', default => sub { {} }, );    # HashRef[Ref]

# fixed and all allele count
has 'allele_of' => (
    is      => 'ro',
    default => sub {
        {   all   => 0,
            drift => 0,
            fixed => 0,
            lost  => 0,
        };
    },
);                                                          # HashRef

sub BUILD {
    my $self = shift;

    # Mersenne Twister generator
    my $rng = Math::Random::MT::Auto->new;
    $rng->srand;
    $self->{rng} = $rng;

    # running time to get to pseudo-equilibrium
    unless ( $self->runtime ) {
        $self->{runtime} = 8 * $self->pop_size;
    }

    # max drifting alleles
    unless ( $self->max_allele ) {
        $self->{max_allele} = $self->pop_size * 2;
    }

    # output directory
    unless ( $self->output ) {
        $self->{output} = sprintf "Freq[N%s][T%s][M%s][Mu%s][Eps%s]", $self->pop_size,
            $self->runtime, $self->max_allele, $self->mu, $self->epsilon;
    }

    # use matrix
    # rows represent alleles, columns represent individual chromosomes
    #           0   1   2   3   4    5   6   7   8   9
    #        mom0   1   2   3   4 dad0   1   2   3   4
    # allele0   1   1   0   0   0    1   1   0   0   0
    # allele1   0   0   0   0   1    1   0   0   0   0
    # allele2   1   0   0   0   0    0   0   0   0   0
    # allele3   1   0   0   0   0    0   0   1   0   1
    $self->{mom} = byte( zeroes( $self->pop_size, $self->max_allele ) );
    $self->{dad} = byte( zeroes( $self->pop_size, $self->max_allele ) );

    $self->{freq} = [ (0) x $self->max_allele ];
    $self->{empty} = [ 0 .. $self->max_allele - 1 ];

    return;
}

sub DEMOLISH {
    my $self = shift;
    undef $self->{rng};
    return;
}

sub run {
    my $self = shift;

    for my $gen ( 0 .. $self->runtime - 1 ) {
        $self->{gen} = $gen;

        # measure frequency and remove all alleles that are fixed
        $self->measure_freq;
        $self->reset_fixed;
        $self->report;

        # gene convert
        $self->gene_convert;

        # mutate
        $self->mutate;

        # reproduce
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
    my $freq = ( $self->mom->daverage + $self->dad->daverage ) / 2;
    $self->{freq} = [ $freq->list ];
    return;
}

sub reset_fixed {
    my $self = shift;

    my $freq = $self->freq;

    # find fixed allele in used
    # In back order, so the indexes don't change
    my $used = $self->used;
    for my $i ( reverse( 0 .. scalar @$used - 1 ) ) {
        my $allele = $used->[$i];
        if ( $freq->[$allele] == 1 ) {
            $freq->[$allele] = 0;
            splice @$used, $i, 1;
            unshift @{ $self->empty }, $allele;
            $self->{mom}->slice(",$allele") .= 0;
            $self->{dad}->slice(",$allele") .= 0;
            $self->{allele_of}->{fixed}++;
        }
        elsif ( $freq->[$allele] == 0 ) {
            splice @$used, $i, 1;
            unshift @{ $self->empty }, $allele;
            $self->{allele_of}->{lost}++;
        }
    }

    $self->{used}               = $used;
    $self->{freq}               = $freq;
    $self->{allele_of}->{drift} = scalar @$used;

    return;
}

sub mu_poisson {
    my $self = shift;
    my $muts = $self->mu * $self->pop_size * 2;
    return $muts == 0 ? 0 : $self->rng->poisson($muts);
}

sub epsilon_poisson {
    my $self  = shift;
    my $count = shift;
    my $cons  = $self->epsilon * $count;
    return $cons == 0 ? 0 : $self->rng->poisson( abs $cons );
}

sub random_pick {
    my $self = shift;
    my $max  = shift;
    return int $self->rng->rand($max);
}

sub random_pick_n {
    my $self = shift;
    my $max  = shift;
    my $n    = shift;
    return map { int $self->rng->rand($max) } ( 1 .. $n );
}

sub new_allele_idx {
    my $self = shift;

    my $idx = shift @{ $self->empty };
    if ( defined $idx ) {
        push @{ $self->used }, $idx;
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
    my $empty = scalar @{ $self->empty };

    return if $empty <= 0;

    if ( $empty <= $muts ) {
        $muts = $empty;
    }
    $self->{allele_of}->{all} += $muts;
    $self->{dynamic}->{ $self->gen }->{newMuts} = $muts;
    print "# Add $muts mutations\n";

    my $max = 2 * $self->pop_size - 1;
    for ( 1 .. $muts ) {
        my $dude_idx   = $self->random_pick($max);
        my $allele_idx = $self->new_allele_idx;
        if ( $dude_idx < $pop_size ) {
            $self->{mom}->slice("$dude_idx,$allele_idx") .= 1;
        }
        else {
            $dude_idx = $dude_idx - $pop_size;
            $self->{dad}->slice("$dude_idx,$allele_idx") .= 1;
        }
    }

    return;
}

sub gene_convert {
    my $self = shift;

    my $epsilon  = $self->epsilon;
    my $pop_size = $self->pop_size;

    my $hetero = which( $self->mom != $self->dad );
    my $count  = $hetero->nelem;

    my $cons = $self->epsilon_poisson($count);

    $self->{dynamic}->{ $self->gen }->{newCons}
        = $epsilon >= 0 ? $cons : -$cons;
    print "# Convert $cons ", $epsilon >= 0 ? "wild types to mutants\n" : "mutants to wild types\n";

    my %seen;
    for ( 1 .. $cons ) {
        my $idx = $self->random_pick($count);
        redo if exists $seen{$idx};
        $seen{$idx}++;
        my $real_idx = $hetero->index($idx);

        my ( $wild_type, $mutant ) = ( "mom", "dad" );
        if ( $self->{mom}->flat->index($real_idx) == 1 ) {
            ( $wild_type, $mutant ) = ( "dad", "mom" );
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

    my $matrix = $self->mom->append( $self->dad );

    # go through each individual in pop
    my $max = 2 * $self->pop_size - 1;

    # pick chromosomes from random parents to mom and dad
    my @moms = $self->random_pick_n( $max, $pop_size );
    my @dads = $self->random_pick_n( $max, $pop_size );

    $self->{mom} = $matrix->dice_axis( 0, \@moms )->sever;
    $self->{dad} = $matrix->dice_axis( 0, \@dads )->sever;

    return;
}

sub report {
    my $self = shift;
    my $gen = shift || $self->gen;

    printf "Gen:[%6s]\tAll:%8d\tLost:%8d\tFixed:%6d\tDrift:%6d\n", $gen,
        $self->{allele_of}->{all}, $self->{allele_of}->{lost},
        $self->{allele_of}->{fixed},
        $self->{allele_of}->{drift};

    $self->{dynamic}->{$gen} = {
        All   => $self->{allele_of}->{all},
        Lost  => $self->{allele_of}->{lost},
        Fixed => $self->{allele_of}->{fixed},
        Drift => $self->{allele_of}->{drift},
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
    for my $allele ( @{ $self->used } ) {
        push @drifting, $self->freq->[$allele];
    }
    my $freq_file = path( $outdir, "freq.csv" );

    print "Write freq file.\n";
    $freq_file->spew( map {"$_\n"} ("freq", @drifting) );

    print "Write freq histogram.\n";
    $self->freq_hist( $freq_file->stringify );

    print "Write opt.yml.\n";
    DumpFile(
        "$outdir/opt.yml",
        {   N   => $self->pop_size,
            T   => $self->runtime,
            M   => $self->max_allele,
            Mu  => $self->mu,
            Eps => $self->epsilon,
        }
    );

    print "Write dynamic.yml.\n";
    DumpFile( "$outdir/dynamic.yml", $self->dynamic );

    print "Write fits files.\n";
    $self->mom->wfits("$outdir/mom.fits");
    $self->dad->wfits("$outdir/dad.fits");

    my $matrix = $self->mom->append( $self->dad );
    my $order  = pdl $self->freq;
    $order = $order->qsorti;
    $matrix = $matrix->dice_axis( 1, $order );
    $matrix->wfits("$outdir/matrix.fits");

    print "Write matrix.png.\n";
    $self->pdl_via_gd( $matrix, "$outdir/matrix.png" );

    return;
}

sub freq_hist {
    my $self = shift;
    my $file = shift;

    my $R = Statistics::R->new;
    $R->set( 'file_csv',   $file );
    $R->set( 'file_chart',   "$file.pdf" );

    my $r_code = <<'EOF';
        library(ggplot2)
        library(scales)
        plotdata <- read.csv(file_csv, header = TRUE)
        pdf( file_chart, width = 3, height = 3)
        ggplot(plotdata, aes(x=freq)) +
            geom_histogram(binwidth = 0.1) +
            theme_bw(base_size = 10)
        dev.off()
EOF

    $R->run($r_code);
    print $R->result;
    $R->stop;
}

sub pdl_via_gd {
    my $self = shift;
    my $pdl  = shift;
    my $file = shift;

    print $pdl;

    my ( $dimx, $dimy ) = $pdl->dims;

    my $image = GD::Image->new( $dimx, $dimy );
    my $white = $image->colorAllocate( 255, 255, 255 );
    my $black = $image->colorAllocate( 0,   0,   0 );

    for my $x ( 0 .. $dimx - 1 ) {
        for my $y ( 0 .. $dimy - 1 ) {
            if ( $pdl->slice("$x,$y") == 1 ) {
                $image->setPixel( $x, $y, $white );
            }
            elsif ( $pdl->slice("$x,$y") == 0 ) {
                $image->setPixel( $x, $y, $black );
            }
            else {
                print $pdl->slice("$x,$y");
            }
        }
    }
    open my $fh, '>', $file;
    binmode $fh;
    print {$fh} $image->png;
    close $fh;
}

package main;

my $markov_sfs = MarkovSFS->new_with_options;
$markov_sfs->run;
