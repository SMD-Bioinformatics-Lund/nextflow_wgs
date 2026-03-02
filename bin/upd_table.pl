#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
open(IN, $ARGV[0]);
my %sum;
my %tot;
while(<IN>) {
    chomp;
    my($chr, $start, $end, $what ) = split /\t/;
    $sum{$chr}->{$what}++;
    $tot{$chr}++;
}

my @chroms = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X);

print "Chromosome\tTotal SNPs\tNon-informative\tMismatch father\tMismatch mother\tAnti-UPD\n";
for my $chr (@chroms) {

    my $paternal = $sum{$chr}->{UPD_PATERNAL_ORIGIN};
    my $maternal = $sum{$chr}->{UPD_MATERNAL_ORIGIN};
    my $anti     = $sum{$chr}->{ANTI_UPD};
    my $tot      = $tot{$chr};

    my $non_inf = $tot-($paternal+$maternal+$anti);
	
    printf "%s\t%d\t%d (%.1f%%)\t%d (%.1f%%)\t%d (%.1f%%)\t%d (%.1f%%)\n", $chr, $tot, $non_inf, 100*$non_inf/$tot, $maternal, 100*$maternal/$tot, $paternal, 100*$paternal/$tot, $anti, 100*$anti/$tot;

}
