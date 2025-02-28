#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use lib dirname(__FILE__);
use vcf2;
use strict;
use Data::Dumper;

my $vcf = CMD::vcf2->new( 'file' => $ARGV[0] );

print $vcf->{header_str};
while ( my $a = $vcf->next_var() ) {

    if ( $a->{INFO}->{CLNSIG} ) {

        #print Dumper($a);
        if ( $a->{GT}->[0]->{GT} eq '0/0' ) {
            next;
        }
        if (   $a->{INFO}->{CLNSIG} eq 'Pathogenic'
            || $a->{INFO}->{CLNSIG} eq 'Likely_pathogenic'
            || $a->{INFO}->{CLNSIG} eq 'Pathogenic/Likely_pathogenic' )
        {
            my @str = split /\t/, $a->{vcf_str};
            print join( "\t", @str[ 0 .. 6 ] ) . "\t";
            print "AC=" . $a->{INFO}->{AC};
            print ";AF=" . $a->{INFO}->{AF};
            print ";DP=" . $a->{INFO}->{DP};
            print ";MQ=" . $a->{INFO}->{MQM};

            print "\t";

            print "GT:AD:DP:GQ:PGT:PID:PL\t";

            print $a->{GT}->[0]->{GT} . ":"
              . "100,100" . ":" . "200" . ":" . "0" . ":" . "0" . ":" . "0"
              . ":" . "0";

            print "\n";
        }

    }
    else {
        next;
    }

}

#GT:AD:DP:GQ:PGT:PID:PL
