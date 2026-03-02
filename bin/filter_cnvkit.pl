#!/usr/bin/env perl
use strict;

use warnings;
my $tsv = $ARGV[0];
my $cov = $ARGV[1];


open (TSV, $tsv) or die $!;

my $row = 0;
while (<TSV>) {
    my @row = split/\t/;
    $row++;
    # print first row, don't do anything more on this row
    if ($row == 1) {print; next;}
    # ignore call outside of genes
    if ($row[3] eq '-') {next;}
    # ignore normal cn calls
    if ($row[6] == 2) {next;}
    if ($row[2] - $row[1] > 1000000) { next; }
    # ignore hetdels with depths below 10% of average
    elsif ($row[6] == 1) {
        if ($row[9]/$cov <= 0.10) {
            next;
        }
        
    }
    print join("\t",@row);
    
}

