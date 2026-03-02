#!/usr/bin/env perl
use strict;


use warnings;
if( is_gzipped($ARGV[0]) ) {
    open( VCF, "zcat $ARGV[0] |" );
}
else {
    open( VCF, $ARGV[0] );
}

while( <VCF> ) {
    if( /^#/ ) {
	print;
	if( /ID=CLNSIG_MOD,/ ) {
            print "##INFO=<ID=SPLICE_INDEL,Number=1,Type=Integer,Description=\"Marks if variant is an INDEL in a splice acceptor or donor region.\">\n";
	}
    }
    else {
	my @a = split /\t/;

	# Print all VCF fields prior to INFO unchanged
	print join "\t", @a[0..6] ;

	if( length($a[3]) == length($a[4]) ) { # Not INDEL
	    print "\t".$a[7];
	}
	else {
	
	    my @INFO = split /;/, $a[7];

	    # Check if splice acceptor/donor variant
	    my $splice_var = 0;
	    foreach( @INFO ) {
		my( $key, $val ) = split /=/;
		if( $key eq "CSQ" ) {
		    if( $val =~ /splice_acceptor_variant|splice_donor_variant/ ) {
			$splice_var = 1;
		    }
		}
	    }

	    # Add INFO flag if splice variant
	    if( $splice_var ) {
		push @INFO, "SPLICE_INDEL=1";
	    }

	    # Print INFO field
	    print "\t";
	    print join ';', @INFO;
	}

	# Print remaining VCF fields unchanged.
	print "\t";
	print join "\t", @a[8..$#a];
    }
}



sub is_gzipped {
    my $fn = shift;
    
    my $file_str = `file $fn`;
    return 1 if $file_str =~ /gzip compressed/;
    return 0;
}
