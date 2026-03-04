#!/usr/bin/env perl 
use strict;
use warnings;
use Data::Dumper;

my $fn = $ARGV[0];

my( $header1, $dels ) = merge("DEL", $fn);
my( $header2, $dups ) = merge("DUP", $fn);

my $rnd = int rand 10000000;
open(TMP, ">tmp.$rnd");
print TMP join("", @$dels);
print TMP join("", @$dups);
close TMP;

print join("", @$header1);
system("sort -k1,1V -k2,2n tmp.$rnd");
unlink("tmp.$rnd");


sub merge {
    my( $type, $fn ) = @_;
    
    my $vcf;
	open($vcf, "zcat $fn |" ) or die "Could not open file $fn";
    
    my @info_fields_to_comma_merge = ("gatkCN");
    
    my %agg_info;
    my @agg;
    my @vars;
    my @header;
    while(<$vcf>) {
	if( /^#/ ) {
	    
	    # Modify header entries to allow multiple values for modified fields
	    if( /<ID=(.*?),/ ) {
		if( grep /^$1$/, @info_fields_to_comma_merge ) {
		    s/Number=1,/Number=\.,/;
		}
	    }
	    push @header, $_;
	    next;
	}
	
	
	my @data = split /\t/;
	my $info = $data[7];
	my @info = split /;/, $info;
	my %info;
	foreach ( @info ) {
	    my( $key, $val ) = split /=/;
	    $info{$key} = ( defined($val) ? $val : "defined" );
	}

	next unless $info{SVTYPE} eq $type;
	
	if( @agg ) {
	    my $dist = abs($data[1] - $agg_info{END});
	    
	    my $agg_rd_first = ( split /,/, $agg_info{gatkCN} )[0];
	    
	    if( $data[0] eq $agg[0] and 
		$dist/abs($info{SVLEN}) < 0.1 and 
		$dist/abs($agg_info{SVLEN}) < 0.1 and
		$dist < 50000 and
		abs($agg_rd_first - $info{gatkCN}) < 0.2 and 
		$info{SVTYPE} eq $agg_info{SVTYPE} ) {
		
		$agg[2] .= ",".$data[2];
		$agg_info{END} = $info{END};
		$agg_info{SVLEN} += $info{SVLEN};
		for (@info_fields_to_comma_merge) { 
		    $agg_info{$_} .= ",".$info{$_};
		}
	    }
	    else {
		# PRINT AGG
		my $var = join("\t", @agg[0..6]);
		$var .= "\t";
		my @info_fields;
		for ("END", "SVLEN", "SVTYPE","gatkCN") {
		    push @info_fields, "$_=$agg_info{$_}";
		}
		$var .= join(";", @info_fields);
		$var .= "\t";
		$var .= join("\t", @agg[8..$#agg]);
		push(@vars, $var);
		@agg = @data;
		%agg_info = %info;
	    }
	}
	else {
	    %agg_info = %info;
	    @agg = @data;
	}
    }
    
	## if no del or no dup ##
	if (!@agg) {
		@vars = ();
		return( \@header, \@vars );
	}
	
    my $var = join("\t", @agg[0..6]);
    $var .= "\t";
    my @info_fields;
    for ("END", "SVLEN", "SVTYPE", "gatkCN") {
	push @info_fields, "$_=$agg_info{$_}";
    }
    $var .= join(";", @info_fields);
    $var .= "\t";
    $var .= join("\t", @agg[8..$#agg]);
    push @vars, $var;


    return( \@header, \@vars );
}


