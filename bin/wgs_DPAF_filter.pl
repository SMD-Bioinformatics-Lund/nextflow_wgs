#!/usr/bin/env perl
use strict;

use warnings;
if( $ARGV[0] =~ /\.gz$/ ) {
    open( VCF, "zcat $ARGV[0]|" );
}
else {
    open( VCF, "$ARGV[0]" );
}

while( <VCF> ) {
    if( /^#/ ) {
	print;
    }
    else {
	chomp;

	my @a = split /\t/;
	my $keep = 1;
	for(9..$#a) {
	    my @a = split /:/, $a[$_];
	    
	    my @allele_counts = split /,/, $a[1];
	    my $dp = $a[2];

	    if( $dp eq "." or $dp == 0 ) {
		$keep = 0;
		last;
	    }
	    elsif( $a[0] ne "0/0" and $allele_counts[1] < 3 ) {
		$keep = 0;
		last;
	    }
	    elsif( $a[0] ne "0/0" and $allele_counts[1]/$dp < 0.2 ) {
		$keep = 0;
		last;
	    }
	}
	print $_."\n" if $keep;
    }
    
}

close VCF;
