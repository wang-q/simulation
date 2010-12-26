#!/usr/bin/perl

package World;
use Moose;
with 'MooseX::Getopt';

has 'size' => ( isa => 'Int', is => 'rw', default => 100, );
has 'rate' => ( isa => 'Num', is => 'rw', default => 0.1, );
has 'target' => (
    isa       => 'Str',
    is        => 'ro',
    default   => 'METHINKS IT IS LIKE A WEASEL',
    metaclass => 'NoGetopt',
);
has current_generation => (
    isa       => 'ArrayRef[Weasel]',
    is        => 'rw',
    builder   => 'first_generation',
    metaclass => 'NoGetopt',
);

sub best {
    my $self = shift;
    $self->current_generation->[0];
}

sub first_generation {
    my $self = shift;

    my @population;
    for ( 1 .. $self->size ) {
        my $individual = Weasel->new(
            target => $self->target,
            rate   => $self->rate,
        );
        push @population, $individual;
    }

    return [ sort { $b->fitness <=> $a->fitness } @population ];
}

sub new_generation {
    my $self = shift;

    my @population;
    for ( 1 .. $self->size ) {
        my $child = $self->best->spawn;
        push @population, $child;
    }
    @population = sort { $b->fitness <=> $a->fitness } @population;

    $self->current_generation( [@population] );
}

sub run {
    my $self = shift;
    $self->new_generation until $self->perfect_offspring;
}

sub perfect_offspring {
    my $self = shift;
    return $self->best->perfect;
}

after new_generation => sub {
    my $self = shift;
    print $self->best->to_string;
};

after run => sub {
    my $self = shift;
    print "\n[Report]\n";
    print "Population size:  ${\$self->size}\n";
    print "Mutation rate:    ${\$self->rate}\n";
    print "Final generation: ${\$self->best->generation}\n";
};

1;

package Mutations;
use Moose::Role;
requires qw(mutate);

use String::Compare;

has fitness => ( isa => 'Num', is => 'rw', lazy_build => 1 );

sub _build_fitness {
    my $self = shift;
    String::Compare::char_by_char( $self->string, $self->target );
}

sub inherit_string {
    my $self = shift;
    return join '',
        map { $self->mutate($_) } 0 .. ( length $self->parent->string ) - 1;
}

sub random_str {
    my $self = shift;
    return join '',
        map { $self->random_char } 0 .. ( length $self->target ) - 1;
}

sub random_char {
    my $self = shift;
    return ( 'A' .. 'Z', ' ' )[ rand(27) ];
}

1;

package NonLockingMutations;
use Moose::Role;
with 'Mutations';

sub mutate {
    my $self   = shift;
    my $idx    = shift;
    my $target = substr( $self->parent->string, $idx, 1 );
    return $target unless rand() < $self->rate;
    return $self->random_char;
}

1;

package LockingMutations;
use Moose::Role;
with 'Mutations';

sub mutate {
    my $self   = shift;
    my $idx    = shift;
    my $target = substr( $self->parent->string, $idx, 1 );
    return $target if $target eq substr( $self->target, $idx, 1 );
    return $target unless rand() < $self->rate;
    return $self->random_char;
}

1;

package Weasel;
use Moose;
with 'NonLockingMutations';

has 'target'   => ( isa => 'Str',    is => 'ro', );
has 'rate'     => ( isa => 'Num',    is => 'rw', );
has parent     => ( isa => 'Weasel', is => 'ro', );
has generation => ( isa => 'Int',    is => 'rw', builder => 'my_generation' );

sub my_generation {
    my $self = shift;
    return 0 unless $self->parent;
    $self->parent->generation + 1;
}

has string => ( isa => 'Str', is => 'ro', lazy_build => 1 );

sub _build_string {
    my $self = shift;
    return $self->inherit_string if $self->parent;
    return $self->random_str;
}

sub _build_generation {
    my $self = shift;
    return 0 unless $self->parent;
    $self->parent->generation + 1;
}

sub perfect {
    my $self = shift;
    return $self->fitness > 0.9999;
}

sub spawn {
    my $self  = shift;
    my $child = Weasel->new(
        parent => $self,
        target => $self->target,
        rate   => $self->rate,
    );
    return $child;
}

sub to_string {
    my $self = shift;
    return
          " ${\sprintf('%04d', $self->generation)}: "
        . "[${ \$self->string }] "
        . "(${\sprintf('%1.4f', $self->fitness)})\n";
}

1;

package main;
my $world = World->new_with_options;
$world->run;

__END__
>perl newweasels.pl --size 50 --rate 0.025
