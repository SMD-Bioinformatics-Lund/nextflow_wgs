#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;


my %opt = ();
GetOptions( \%opt, 'vcf=s', 'mother=s', 'father=s', 'out=s' );
my $vcf_file = $opt{vcf};
my $mother = $opt{mother};
my $father = $opt{father};
my $out = $opt{out};
open (OUT,'>',$out);

my @extra_individuals;
if ($father ne "null") {
	push @extra_individuals,$father;
}
if ($mother ne "null") {
	push @extra_individuals,$mother;
}

my $vcf = CMD::vcf2->new( 'file'=>$vcf_file );
my @header = split/\n/,$vcf->{header_str};

## Add the relatives to VCF-header
foreach my $line (@header) {
	if ($line =~ /^#CHROM/) {
		$line = $line."\t".join("\t",@extra_individuals);
	}
	print OUT $line,"\n";
}
while ( my $a = $vcf->next_var() ) {
	## Save all keys of proband GT
	my @keys;
	foreach my $key (keys %{ $a->{GT} ->[0] }) {
		push @keys,$key;
	}
	## Create empty hashes for mother and father
	## Set _sample_id and GT.
	for (my $i = 1; $i <= scalar(@extra_individuals); $i++ ) {
		foreach my $s (@keys) {
			#print $s,"\n";
			if ($s eq 'GT') {
				$a->{GT} ->[$i] -> {$s} = './.';
			}
			elsif ($s eq '_sample_id') {
				$a->{GT} ->[$i] -> {$s} = $extra_individuals[$i-1];
			}
			else {
				$a->{GT} ->[$i] -> {$s} = '.';
			}
			
		}
		
	}
	#print Dumper($a->{GT});

	my $vcf_str = vcfstr($a,[]);
	print OUT $vcf_str;
}
close OUT;

sub vcfstr {
	my( $v, $sample_order ) = @_;
	
	my @all_info;
	my $tot_str = $v->{CHROM}."\t".$v->{POS}."\t".$v->{ID}."\t".$v->{REF}."\t".$v->{ALT}."\t".$v->{QUAL}."\t".$v->{FILTER}."\t";

	# Generate and print INFO field
	for my $info_key (@{$v->{INFO_order}}) {
		if($info_key eq "CSQ") {
			push @all_info, $info_key."=".$v->{_CSQstr};
		}
		else {
			push @all_info, $info_key."=".$v->{INFO}->{$info_key};
		}
	}
	$tot_str = $tot_str.join(";", @all_info)."\t";

	# Print FORMAT field
	$tot_str = $tot_str.join(":", @{$v->{FORMAT}})."\t";


	my %order;
	my $i=0;
	if( @$sample_order > 0 ) {
		$order{$_} = $i++ foreach @{$sample_order};
	}
	else {
		$order{$_->{_sample_id}} = $i++ foreach @{$v->{GT}};
	}
	# Print GT fields for all samples
	for my $gt ( sort {$order{$a->{_sample_id}} <=> $order{$b->{_sample_id}}} @{$v->{GT}}) {
		my @all_gt;
		for my $key ( @{$v->{FORMAT}} ) {
			push @all_gt, ( defined $gt->{$key} ? $gt->{$key} : "");
		}
		$tot_str = $tot_str.join(":", @all_gt)."\t";
	}
	$tot_str = $tot_str."\n";
	return $tot_str;
}