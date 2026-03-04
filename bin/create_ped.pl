#! /usr/bin/env perl 

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;


my %opt = ();
GetOptions( \%opt, 'father=s', 'mother=s', 'id=s', 'group=s', 'sex=s');
my %const;
my @lines;
my $father;
my $mother;
push @lines,$opt{id};
unless ($opt{father} eq '0') {
    $father = $opt{father};
    push @lines,$opt{father};
}
unless ($opt{mother} eq '0') { 
    $mother = $opt{mother};
    push @lines,$opt{mother};
}
my $sex;
if ($opt{sex} eq 'F') { $sex = 2;}
elsif ($opt{sex} eq 'M') {$sex = 1;}
else {$sex = 1;}

open (BASE, '>',$opt{group}."_base.ped" ) or die $!;

if (@lines > 1) {
    if ($father) {
        my $group = $opt{group}."_fa.ped";
        open (FA, '>', $group) or die $!;
        # trio
        if ($mother) {
            $group = $opt{group}."_ma.ped";
            open (MA, '>', $group) or die $!;
            print BASE $opt{group}."\t".$opt{id}."\t$opt{father}\t$opt{mother}\t".$sex."\t2\n";
            print BASE $opt{group}."\t".$opt{father}."\t0\t0\t1\t1\n";
            print BASE $opt{group}."\t".$opt{mother}."\t0\t0\t2\t1\n";
            
            print FA $opt{group}."\t".$opt{id}."\t$opt{father}\t$opt{mother}\t".$sex."\t2\n";
            print FA $opt{group}."\t".$opt{father}."\t0\t0\t1\t2\n";
            print FA $opt{group}."\t".$opt{mother}."\t0\t0\t2\t1\n";

            print MA $opt{group}."\t".$opt{id}."\t$opt{father}\t$opt{mother}\t".$sex."\t2\n";
            print MA $opt{group}."\t".$opt{father}."\t0\t0\t1\t1\n";
            print MA $opt{group}."\t".$opt{mother}."\t0\t0\t2\t2\n";
            close FA;
            close MA;
            close BASE;

        }
        # duo with father
        else {
            print BASE $opt{group}."\t".$opt{id}."\t$opt{father}\t0\t".$sex."\t2\n";
            print BASE $opt{group}."\t".$opt{father}."\t0\t0\t1\t1\n";

            print FA $opt{group}."\t".$opt{id}."\t$opt{father}\t$opt{mother}\t".$sex."\t2\n";
            print FA $opt{group}."\t".$opt{father}."\t0\t0\t1\t2\n";
            close FA;
            close BASE;
        }
    }
    # duo with mother
    else {
        my $group = $opt{group}."_ma.ped";
        open (MA, '>', $group) or die $!;
        print BASE $opt{group}."\t".$opt{id}."\t0\t$opt{mother}\t".$sex."\t2\n";
        print BASE $opt{group}."\t".$opt{mother}."\t0\t0\t2\t1\n";
        
        print MA $opt{group}."\t".$opt{id}."\t0\t$opt{mother}\t".$sex."\t2\n";
        print MA $opt{group}."\t".$opt{mother}."\t0\t0\t2\t2\n";
        close MA;
        close BASE;
    }

}

## Single sample, no parents, affected child
else {
    print BASE $opt{group}."\t".$opt{id}."\t0\t0\t".$sex."\t2\n";
    close BASE;
}









