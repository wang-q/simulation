#!/usr/bin/perl -w
##############################
#	Overview
##############################
# This program implements a simple yet easily generalizable genetic algorithm, which is modified
#	by altering certain parameters and providing a fitness function. 
#
#		By Seth Just <seth.just@gmail.com> Dec. 10, 2008
#		Licensed under Creative Commons Non-Commercial, Share-Alike, Attribution license
#		You may distribute this file for non-commercial purposes
#			as long as this message remains in place
#			
# BASIC OUTLINE:
#	1. A randomized population is created. An array holds references to the specified number of
#		"individuals". Each individual is an array of alleles (scalar values) that are chosen by
#		evaluating a specified string. These are currently positive real numbers.
#	2. A loop runs for a specified number of iterations. Each iteration mutates the population
#		according to specified parameters, and then breeds the entire population as described below.
#
# WRITING A FITNESS FUNCTION:
#	The fitness function takes a reference to an individual, which is an array of scalars. It is
#		expected to return a scalar value that is used to sort the individuals by performance.
#		To write a fitness function, either iterate over the individual (Ex. 1)  or by addressing
#		each allele individually (Ex. 2).
#	Example 1:
#	 	for (@$ind) {					# This is a standart problem: the sphere model
# 			$result += ($_ - 1)**2;		#	Minimize f(x) = sum((x(i)-1)^2)
# 		}
#	Example 2:
#		$result = sqrt($ind->[0]**2+$ind->[1]**2);
#
# FURTHER MODIFICATION:	
#	It is possible through modification of the initialization, fitness, and possibly 
#		mutation/breeding routines, to apply this algorithm to any problem. For example, it is
#		possible to treat individuals as strings simply by using chr(int(255*$ind[$i])). In cases
#		like this it is particularly important to play with $mut_weight and $mut_offset in order to
#		allow mutations to be picked up on by the breeding algorithm.

##############################
#	Pragmas
##############################

use strict;
use warnings;

##############################
#	Modules
##############################

# Add any needed includes here

##############################
#	Global Variables
##############################

my $population 	= 	100;			# Number of individuals
my $alleles		=	5;				# Number of alleles for each individual
my $init = sub { rand(10) };		# Each allele is initialized by executing this

my $iterations 	= 	500;			# Number of iterations to run

my $mut_weight	=	.01;			# Random mutation parameters
my $mut_offset	=	.01;

my $reset_prob	=	0.0001;			# Percentage of allele crosses that result in re-initialization
									# 	instead of crosing. As a rule, keep this very, very low

my $minimize	=	1;				# Set true to minimize fitness, otherwise it maximizes

my $verbose 	= 	0;				# Level of verbosity:
									# 	1 => Print mean and best fitnesses
									# 	2 => Print average characteristics, mean, and best
									# 	3 => Print individuals, average characteristics, mean, and best

##############################
#	Fitness Function
##############################

sub fitness{
	my $ind = shift;
	my $result = 1;

	for (@$ind) {					# This is a standard problem: the sphere model
		$result += ($_ - 1)**2;		#	Minimize f(x) = sum((x(i)-1)^2)
	}

	return $result;
}

##############################
#	Subroutines
##############################

sub gen_indvs{	# Initializes the individuals list with randomly initialized individuals
	my $indvs = shift;
	for (1..$population) {
		my @temp = ();
		for (1..$alleles) {
			push @temp, $init->();
		}
		push(@$indvs, \@temp);
	}
	return 0;
}

sub mutate{		# Mutates each allele of an individual by a random, weighted amount,
				#	controlled by certain parameters. It re-initializes if the allele
				#	is zero, and takes abs()
	my $ind = shift;
	
	for (@$ind) {
		$_ = $init->() unless ($_);
		$_ += (rand(2) - 1)*($mut_weight*($_)+$mut_offset);
		$_ = -$_ unless ($_>0);
	}
	
	return 0;
}

sub breed{		# Breeds each individual with the top 20 percent of the population randomly
				#	by averaging each allele of the individual with the corresponding allele
				#	of a random individual in the top quintile. Note that this means that
				#	each individual breeds with, at most (and ideally), as many individuals 
				#	as it has alleles, which breaks down the parent / child model slightly.
	my $indvs = shift;
	my @list;
	($minimize) || (@list = sort {fitness($b) <=> fitness($a)} @$indvs);	# Sort asc. vs. desc.
	($minimize) && (@list = sort {fitness($a) <=> fitness($b)} @$indvs);
	
	@list = @list[1..(int(scalar(@list)/5)+1)];
	
	for (@$indvs) {				# Iterate through individuals
		my $ind_ref = $_;
		my $i = 0;
		for (@$ind_ref) {		# Iterate through alleles
			$_ = ((@list[int(rand(length(@list)-1))])->[$i] + $_)/2;	# Average given allele with 
																		#	a random, fit,
																		#	corresponding allele
			$_ = (rand(1)<$reset_prob?$init->():$_);			# Re-init the allele some
																		#	percent of the time
			$i++;
		}
	}
	
	return 0;
}

sub average{	# Caluclates and prints the fitness of the average individual.
				#	NOTE: This is not the average fitness, it is the fitness of an average
				#	individual. In some cases this is very misleading.
	my $indvs = shift;
	my @average;
	
	for (@$indvs) {
		my $index = 0;
		for (@$_) {
			$average[$index] += $_;
			$index++;
		}
	}
	
	for (@average) {
		$_ /= scalar @$indvs;
		($verbose>=2) && print $_ . ' 'x(25-length($_));
	}
	
	my $fitness = fitness(\@average);
	
	print $fitness . ' 'x(25-length($fitness));
	return 0;
}

sub best_fit{	# Finds the best individual and prints its fitness
	my $indvs = shift;
	my @list;
	
	($minimize) || (@list = sort {fitness($b) <=> fitness($a)} @$indvs);
	($minimize) && (@list = sort {fitness($a) <=> fitness($b)} @$indvs);
	
	print fitness($list[0]) . "\n";
	return 0;
}

sub best{		# Finds the best individual and prints its characteristics
	my $indvs = shift;
	my @list;
	
	($minimize) || (@list = sort {fitness($b) <=> fitness($a)} @$indvs);
	($minimize) && (@list = sort {fitness($a) <=> fitness($b)} @$indvs);
	
	my $best = $list[0];
	
	for (@$best){
		print "$_\n";
	}
	return 0;
}

##############################
#	Main
##############################

main: {
	my @indvs;
	gen_indvs(\@indvs);
		
	for (0..$iterations) {			# Main loop set to number of iterations
		($verbose>=1) && print $_ . "\t";
		
		for (@indvs) {				# Iterate through the individuals		
			mutate($_);				#	and mutate them
			
			($verbose>=3) && do { foreach (@$_){
				print "$_," . ' 'x(25-length($_));
			}};
		}
		
		($verbose>=1) && average(\@indvs);
		($verbose>=1) && best_fit(\@indvs);
		
		breed(\@indvs);				# Breed the generation
	}								# END main loop
	
	($verbose) || do{
		average(\@indvs);
		best_fit(\@indvs);
	};
	
	best(\@indvs);
}