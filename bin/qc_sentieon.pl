#!/usr/bin/env perl
use warnings;
use strict;
use Data::Dumper;
use JSON;# qw( encode_json );
use Getopt::Long;

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# General QC-script for sentieon-data. Takes two arguments: SAMPLE-ID, TYPE(panel or wgs)
# PANEL: Requires the following algos from sentieon qc:
#        --algo MeanQualityByCycle mq_metrics.txt \\
#        --algo QualDistribution qd_metrics.txt \\
#        --algo GCBias --summary gc_summary.txt gc_metrics.txt \\
#        --algo AlignmentStat aln_metrics.txt \\
#        --algo InsertSizeMetricAlgo is_metrics.txt \\
#        --algo HsMetricAlgo --targets_list $target_intervals --baits_list $target_intervals hs_metrics.txt \\
#        --algo CoverageMetrics --cov_thresh 1 --cov_thresh 10 --cov_thresh 30 --cov_thresh 100 --cov_thresh 250 --cov_thresh 500 cov_metrics.txt
# WGS: Requires the following algos from sentieon qc
#       --algo MeanQualityByCycle mq_metrics.txt \\
#        --algo QualDistribution qd_metrics.txt \\
#        --algo GCBias --summary gc_summary.txt gc_metrics.txt \\
#        --algo AlignmentStat aln_metrics.txt \\
#        --algo InsertSizeMetricAlgo is_metrics.txt \\
#        --algo WgsMetricsAlgo wgs_metrics.txt \\
#
# BOTH options also need the dedup_metrics from the locus_collector + dedup steps
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


my $SID;
my $type;

my $align_metrics_file;
my $insert_file;
my $dedup_metrics_file;
my $metrics_file;
my $gcsummary_file;
my $coverage_file_summary;
my $coverage_file;

GetOptions(
    "SID=s"                   => \$SID,
    "type=s"                  => \$type,
    "align_metrics_file=s"    => \$align_metrics_file,
    "insert_file=s"           => \$insert_file,
    "dedup_metrics_file=s"    => \$dedup_metrics_file,
    "metrics_file=s"          => \$metrics_file,
    "gcsummary_file=s"        => \$gcsummary_file,
    "coverage_file_summary=s" => \$coverage_file_summary,
    "coverage_file=s"         => \$coverage_file
);

unless ($SID && $type && $align_metrics_file && $insert_file && $dedup_metrics_file && $metrics_file && $gcsummary_file) {
    die "Usage: $0 --SID <sample_id> --type <panel|wgs> --align_metrics_file <file> --insert_file <file> --dedup_metrics_file <file> --metrics_file <file> --gcsummary_file <file> [--coverage_file file] [--coverage_file_summary file]\n"
}

unless ($type eq "panel" || $type eq "wgs") { die "--type must be 'panel' or 'wgs'" }
unless (-f $align_metrics_file) { die "--align_metrics_file does not point to a file" }
unless (-f $insert_file)        { die "--insert_file does not point to a file" }
unless (-f $dedup_metrics_file) { die "--dedup_metrics_file does not point to a file" }
unless (-f $metrics_file)       { die "--metrics_file does not point to a file" }
unless (-f $gcsummary_file)     { die "--gcsummary_file does not point to a file" }

if ($type eq "panel") {
    unless ($coverage_file && $coverage_file_summary) {
        die "If running a panel, --coverage_file and --coverage_file_summary must be provided"
    }
    unless (-f $coverage_file) { die "--coverage_file_summary does not point to a file" }
    unless (-f $coverage_file_summary) { die "--coverage_file_summary does not point to a file" }
}

my %pct_above_x;
my $median;
my $pct_above_panel;
if ($type eq "panel") {

    ($pct_above_panel, $median) = coverage_calc();
    %pct_above_x = %$pct_above_panel;
}

my %results;

if ($type eq "wgs") {
    my ( $sum, %quartiles, $pct25_obs, $pct50_obs, $pct75_obs );
    open( HS, $metrics_file );
    while( <HS> ) {
  
        if( /^\#SentieonCommandLine/ ) {
	    <HS>;
            my $vals = <HS>;
            my @a = split /\t/, $vals;
            $results{'median_cov'} = $a[3];
            $results{'sd_coverage'} = $a[2];
            $results{'mean_coverage'} = $a[1];
            $pct25_obs = ( $a[0] + 1 ) / 4; # get observations from different quartiles
            $pct50_obs = ( $a[0] + 1 ) / 2; # get observations from different quartiles
            $pct75_obs = ( ( $a[0] + 1 ) * 3 ) / 4; # get observations from different quartiles
            $pct_above_x{'1'} = $a[12];
            $pct_above_x{'5'} = $a[13];
            $pct_above_x{'10'} = $a[14];
            $pct_above_x{'15'} = $a[15];
            $pct_above_x{'20'} = $a[16];
            $pct_above_x{'25'} = $a[17];
            $pct_above_x{'30'} = $a[18];
            $pct_above_x{'40'} = $a[19];
            $pct_above_x{'50'} = $a[20];
            $pct_above_x{'60'} = $a[21];
            $pct_above_x{'70'} = $a[22];
            $pct_above_x{'80'} = $a[23];
            $pct_above_x{'90'} = $a[24];
            $pct_above_x{'100'} = $a[25];

	}elsif( /^(\d+)\s+(\d+)/ ){
            $sum += $2;
            if( $sum >= $pct25_obs and not $quartiles{ 'R_25' } ){
                $quartiles{ 'R_25' } = $1;
            }
            if( $sum >= $pct50_obs and not $quartiles{ 'R_50' } ){
                $quartiles{ 'R_50' } = $1;
            }
            if( $sum >= $pct75_obs and not $quartiles{ 'R_75' } ){
                $quartiles{ 'R_75' } = $1;
            }
        }
    }
    close HS;
    $results{'iqr'} = ( $quartiles{ 'R_75' } - $quartiles{ 'R_25' } );

    open( GC, $gcsummary_file );
    while( <GC> ) {
        if( /^\#SentieonCommandLine/ ) {
	    <GC>;
            my $vals = <GC>;
            my @a = split /\t/, $vals;
            $results{'at_drop'} = $a[5];
            $results{'gc_drop'} = $a[6];
        }
    }
    close GC;

}
elsif ($type eq "panel") {
    open( HS, $metrics_file );
    while( <HS> ) {
        if( /^\#SentieonCommandLine/ ) {
	    <HS>;
            my $vals = <HS>;
            my @a = split /\t/, $vals;
            $results{'pct_on_target'} = $a[18];
            #print "pct_on_target: $a[18]\n";
            $results{'fold_enrichment'} = $a[25];
            #print "fold_enrichment: $a[25]\n";
            $results{'mean_coverage'} = $a[22];
            #print "median_coverage: $a[22]\n";
            $results{'fold_80'} = $a[32];
            #print "fold_80: $a[32]\n";
            $results{'at_drop'} = $a[48];
            $results{'gc_drop'} = $a[49];
	    }
    }
    close HS;
    $results{'median_cov'} = $median;
}
else { print STDERR "panel or wgs\n"; exit;}


## INSSIZE ##
open( INS, $insert_file );
while( <INS> ) {
	if( /^\#SentieonCommandLine/ ) {
	    <INS>;
	    my $vals = <INS>;
	    my @a = split /\t/, $vals;
        my $ins_size = sprintf "%.0f", $a[4];
	    $results{'ins_size'} = $ins_size;
        #print "ins_size: $a[4]\n";
        my $ins_size_dev = sprintf "%.0f", $a[5];
	    $results{'ins_size_dev'} = $ins_size_dev;
        #print "ins_size_dev: $a[5]\n";
	}
}
close INS;


## DEDUP ##

open( DEDUP, $dedup_metrics_file );
while( <DEDUP> ) {
    if( /^\#SentieonCommandLine/ ) {
	    <DEDUP>;
	    my $vals = <DEDUP>;
	    my @a = split /\t/, $vals;
	    $results{'dup_reads'} = $a[6];
        #print "dup_reads: $a[6]\n";
	    $results{'num_reads'} = $a[2];
        #print "num_reads: $a[2]\n";
        $results{'dup_pct'} = $a[8];
        #print "dup_pct: $a[8]\n";
        my $mapped = $a[2]-$a[4];
        $results{'mapped_reads'} = $mapped;
        #print "mapped_reads: $mapped\n";
	}
}
close DEDUP;

## ALIGMENTMETRICS ##

open( ALIGN, $align_metrics_file );
while( <ALIGN> ) {
    if( /^\#SentieonCommandLine/ ) {
	    <ALIGN>;
	    my $vals = <ALIGN>;
	    my @a = split /\t/, $vals;
	    $results{'pf_mismatch_rate'} = $a[12];
	    $results{'pf_error_rate'} = $a[13];
	}
}
close ALIGN;

sub coverage_calc {

    my @cov;
    open( COV, $coverage_file );
    while( <COV> ) {

        unless( /^Locus/ ) {
            my $vals = <COV>;
            my @a = split /\t/, $vals;
            push @cov,$a[1];
        }
        
    }
    my @sorted_cov = sort @cov;

    my $len_cov = scalar(@cov);
    my $median;
    if (($len_cov / 2) =~ /\d+\.\d+/) { 
        my $median_index = ($len_cov / 2) + 0.5;
        $median = $sorted_cov[$median_index-1];
    }
    else {
        my $median_index = ($len_cov / 2);
        $median = ($sorted_cov[$median_index-1] + $sorted_cov[$median_index-1]) / 2;
    }


    ## cov_% ##
    my %pct_above_x;
    open( COV_THRESH, $coverage_file_summary );
    while( <COV_THRESH> ) {
        if( /^sample_id/ ) {
            my $vals = <COV_THRESH>;
            my @a = split /\t/, $vals;
            chomp @a;
            $pct_above_x{'1'} = $a[6];
            $pct_above_x{'10'} = $a[7];
            $pct_above_x{'30'} = $a[8];
            $pct_above_x{'100'} = $a[9];
            $pct_above_x{'250'} = $a[10];
            $pct_above_x{'500'} = $a[11];
        }
    }
    close COV_THRESH;
    return \%pct_above_x, $median;

}

$results{'pct_above_x'} = \%pct_above_x;
$results{'sample_id'} = $SID;
my $json = JSON->new->allow_nonref;
print $json->pretty->encode( \%results );
