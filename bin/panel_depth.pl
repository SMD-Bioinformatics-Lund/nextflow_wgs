#!/usr/bin/env perl
use strict;

use warnings;
my $cutoff = 500;
die "USAGE: panel_depth.pl BAM BED\n" unless @ARGV == 2;

my( $bam, $bed ) = ( $ARGV[0], $ARGV[1] );

die "No file $bam" unless -s $bam;
die "No file $bed" unless -s $bed;


open( DEPTH, "sambamba depth base $bam -L $bed |" );

my( $start_pos, $start_chr, $last_low_pos, $last_low_chr, $low_cov_sum );
while( <DEPTH> ) {
    my @a = split /\t/;
    my( $chr, $pos, $depth ) = ( $a[0], $a[1], $a[2] );

    if( $depth < $cutoff ) {

	# Prev low position was right before
	if( $last_low_chr and $last_low_pos and $last_low_chr eq $chr and $last_low_pos == $pos-1 ) {
	    # Skip along low depth region
	    $low_cov_sum += $depth;
	}

	# Prev low postion was somewhere else
	else {
	    if( $start_pos and $start_chr ) {
		print $start_chr."\t".$start_pos."\t".$last_low_pos."\t".($low_cov_sum/($last_low_pos-$start_pos+1))."\n";
	    }
	    $start_chr = $chr;
	    $start_pos = $pos;
	    $low_cov_sum = $depth;
	}
	$last_low_chr = $chr;
	$last_low_pos = $pos;
    }
    else {
	if( $start_chr and $start_pos ) {
	    print $start_chr."\t".$start_pos."\t".$last_low_pos."\t".($low_cov_sum/($last_low_pos-$start_pos+1))."\n";
	    undef $start_chr;
	    undef $start_pos;
	}
    }
}
