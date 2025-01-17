#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

our %opt;
GetOptions( \%opt, 'cnv=s', 'pb=s', 'help' );
my $SCORE = 10;

print_usage("ERROR: --cnv and --pb required!\n") unless $opt{cnv} and $opt{pb};

unless( -e $opt{cnv} ) { print "ERROR: VCF file not found\n"; exit; }

if( defined($opt{cnv}) ){

    if( $opt{cnv} =~ /gz$/ ){
        open( VCF, "zcat $opt{cnv} |" ) or die "gunzip $opt{cnv}: $!";
    }else{
        open( VCF, "$opt{cnv}" ) or die "gunzip $opt{cnv}: $!";
    }
    my %labels;
    while ( my $row = <VCF> ){
        chomp( $row );
        my @cols = split(/\t/, $row);
        next if( $row=~/^##/ );
        if( $row=~/^#CHROM/ ){ 
            for my $index(0..$#cols){ $labels{ $cols[$index] } = $index; }
            if( not defined$labels{ $opt{pb} } ){ die "$opt{pb} is not recognized in file!\n"; }
            print "track db=\"hg38\" name=\"cnv_regions\" description=\"cnvs extracted from vcf using cnv2bed.pl\"\n"; 
        }else{
            my ( $svtype, $end, $score, $rank );
            my $color = '192,192,192';
            next if( $row=~/^Y/ );
            next if( $row=~/^M/ );
            next if( $cols[ $labels{ $opt{pb} } ] =~/^\.\// );
            if( $row =~/SVTYPE=(\w+);/ ) { $svtype = $1; next if( $svtype eq 'BND' ); }
            if( $row =~/END=(\d+);/ ){ $end = $1; }
            if( $row =~/RankScore=$opt{pb}:([\d-.]+);/ ){ $score = $1; next if( $score <= $SCORE ); }
            if( $row =~/RankResult=([\d-]+)\|/ ){ $rank = $1; }
            if( $svtype =~/DEL/ ){ $color = '204,0,0'; }elsif( $svtype =~/DUP/ ){ $color = '0,0,153'; }
            print "chr$cols[0]\t".( $cols[1] - 1 )."\t$end\t\"$svtype,score=$score\"\t0\t\.\t0\t0\t$color\n";
        }
    } close VCF;
}


sub print_usage {
    print "$_[0]\n\n" if $_[0];
    print "USAGE: cnv2bed.pl --cnv <VCF FILE>\n\n";
    print "    --cnv     FILE       Path to CNV VCF file \(required\)\n\n";
    print "    --pb      SAMPLE     Proband sample name \(required\)\n\n";
    print "    --help               Will print this message\n\n";
    exit(0);
} # print_usage
