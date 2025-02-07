#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

my %opt = ();
GetOptions( \%opt, 'merged_panel_sv_vcf=s' );

my $vcf = CMD::vcf2->new('file'=>$opt{merged_panel_sv_vcf} );

my @header = split/\n/,$vcf->{header_str};
foreach my $header (@header) {
        print $header."\n";
}

while ( my $a = $vcf->next_var() ) {

    my @str = split/\t/,$a->{vcf_str};
    ## Print first 7 columns of vcf
    print join("\t",@str[0..6])."\t";

    my @INFO = split/;/,$str[7];

    print join(';',@INFO)."\t";
    print $str[8]."\t";

    ## 0/0 genotypes cannot be loaded as artefacts, if 0/0 set as 0/1
    my $gtmod = $str[9];
    $gtmod =~ s/^0\/0/0\/1/;
    print $gtmod;
    print "\n";

}
