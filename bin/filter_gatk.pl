#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

my $vcf = CMD::vcf2->new('file'=>$ARGV[0] );

my @header = split/\n/,$vcf->{header_str};
my $count = 1;
foreach my $header (@header) {    
    if ($header =~ /##FORMAT=<ID=GT/) {
        print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Segment genotype\">\n";
    }
	elsif ($header =~ /##INFO/ && $count == 1) {
		print $header."\n";
		print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
		print "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
		print "##INFO=<ID=gatkCN,Number=1,Type=Integer,Description=\"estimated copy number 0-5 gatk\">\n";
		$count ++;
	}
    else {
        print $header."\n";
    }

}
while ( my $a = $vcf->next_var() ) {

    #next if ($a->{GT}->[0]->{QS} < 1000);
    next if ($a->{GT}->[0]->{GT} == 0);
    #print Dumper($a);
    ## create genotypes matching manta/tiddit
    if ($a->{GT}->[0]->{GT} == 1) {
        if ($a->{GT}->[0]->{CN} == 1) {
            $a->{GT}->[0]->{GT} = "0/1";
        }
        elsif ($a->{GT}->[0]->{CN} == 0) {
            $a->{GT}->[0]->{GT} = "1/1";
        }
        else {
            $a->{GT}->[0]->{GT} = "0/1";
        }
    }
	my $type = $a->{ALT};
	$type =~ s/[>|<]//g;
    $a->{INFO}->{SVTYPE} = $type;
	$a->{INFO}->{SVLEN} = $a->{INFO}->{END} - $a->{POS} +1;
	$a->{INFO}->{gatkCN} = $a->{GT}->[0]->{CN};
	push (@{$a->{INFO_order}},"SVTYPE");
	push (@{$a->{INFO_order}},"SVLEN");
	push (@{$a->{INFO_order}},"gatkCN");
    my $str = vcfstr($a, []);
    print $str;
}


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
