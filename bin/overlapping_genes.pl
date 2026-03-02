#!/usr/bin/env perl
use strict;

use warnings;
my $in_bed = $ARGV[0];
my $genes_bed = $ARGV[1];

my $rnd = int rand 1000000000;
my $tmp_infile = "input.$rnd.bed";
system("grep -v '^\@' $in_bed| grep -v ^CONTIG| grep -v ^REF > $tmp_infile");
my @overlap = `bedtools intersect -a $tmp_infile -b $genes_bed -loj`;
unlink $tmp_infile;

foreach my $line (@overlap) {
    chomp $line;
    my @f = split /\t/, $line;
    print "$f[0]\t$f[1]\t$f[2]\t$f[3]\t$f[-1]\n";
}

