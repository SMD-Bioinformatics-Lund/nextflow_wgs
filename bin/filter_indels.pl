#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;



my $vcf = CMD::vcf2->new('file'=>$ARGV[0] );

print $vcf->{header_str};
while ( my $a = $vcf->next_var() ) {

    
    my $AF = ( $a->{INFO}->{CSQ}->[0]->{'gnomADg_AF'} or 0 );
    my @AF = split '&',$AF;
    if ( scalar(@AF) >= 2 ) {
        $AF = findmax(@AF);
    }
    if ($AF <= 0.05 ) {
        print $a->{vcf_str};
        print "\n";
    }
   
}

sub findmax {
    my @in = @_;
    my $high = 0;
    foreach my $val (@in) {
        if ($val && $val ne '.') {
           if($val > $high) {
               $high = $val;
           }
       }
       else {

       }
    }
    return $high;
}