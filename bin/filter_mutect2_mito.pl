#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

# Variants are only interesting if they exist in proband. A 0/0 in proband is non-informative
# Also require atleast 50 basepairs depth for sight.
# Tag variants with a placeholder genotype quality to make downstream loqusdb accept variant


my $line = 1;
my $index;
my $vcf = CMD::vcf2->new('file'=>$ARGV[0] );
my $proband = $ARGV[1];
print $vcf->{header_str};
while ( my $vcf_rec = $vcf->next_var() ) {    
    if ($line == 1) {
        my $count = 0;
        foreach my $ind (@{ $vcf_rec->{GT} }) {
            if ($ind->{_sample_id} eq $proband) {
                $index = $count;
            }
            $count++;
        }
    }
    ## variant exist in proband
    if ($vcf_rec->{GT}->[$index]->{GT} =~ /1/) {
        ## variant has at least 50 depth    
        if ($vcf_rec->{GT}->[$index]->{DP} > 50) {
            my $quality = 99;
            my $af_thres = 0.99;
            my $gt = genotype_quality($vcf_rec->{GT},$quality,$af_thres);
            $vcf_rec->{GT} = $gt;
            push (@{$vcf_rec->{FORMAT}},"GQ");
            my $vcfstr = vcfstr($vcf_rec,[]);
            print $vcfstr;
        }
    }
    $line++;
}

sub genotype_quality {
    # Add placeholder 99 GQ to each individual
    # Also update genotype to 1/1 if AF>99%
    my $gt = shift;
    my $quality = shift;
    my $af_thres = shift;
    my $c = 0;
    foreach my $ind ( @{ $gt }) {
       $gt->[$c]->{GQ} = $quality;
       if ($ind->{AF} >= $af_thres) {
            $gt->[$c]->{GT} = '1/1';
       }
       $c++;
    }
    return $gt;
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