#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use lib dirname (__FILE__);
use vcf2;
use strict;
use Data::Dumper;

my $vcfheader = $ARGV[0];
my $id = $ARGV[1];


my $alu = "ALU.final_comp.vcf";
my $line1 = "LINE1.final_comp.vcf";
my $sva = "SVA.final_comp.vcf";

my %qcval = (
    'lc' => 'Low Complex Region',
    's25' => 'Greater than 25.0% of samples do not have data',
    'hDP' => '"More than the expected number of discordant pairs at this site are also split',
    'PASS' => 'PASS',
    'rSD' => "Ratio of LP to RP is greater than 2.0 standard deviations"
);
my %assess = (
    0 => 'No overlapping reads at site',
    1 => 'Imprecise breakpoint due to greater than expected distance between evidence',
    2 => 'discordant pair evidence only -- No split read information',
    3 => 'left side TSD evidence only',
    4 => 'right side TSD evidence only',
    5 => 'TSD decided with split reads -> highest possible quality'
);

my @files;
push @files,$alu,$line1,$sva;

my @non_empty;
foreach my $file (@files) {
    my $size = -s $file;
    if ($size < 1) { next; }

    my $vcf = CMD::vcf2->new('file'=>$file);
    my $out = $file."mod";
    open (OUT, '>' ,$out);
    my @header = split/\n/,$vcf->{header_str};
    my $count = 1;
    foreach my $header (@header) {
        $header =~ s/\'\\\|\'/pipe-sign/g;
        
        if ($header =~ /^##INFO/ && $count == 1) {
            print OUT '##INFO=<ID=END,Number=.,Type=Integer,Description="END position set to start position for insertions">'."\n";
            print OUT '##INFO=<ID=SCOUT_CUSTOM,Number=.,Type=String,Description="Custom annotations for scout">'."\n";
            $count++;
        }
        if ($header =~ /^#CHROM/) {
            print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$id\n";
        }
        else {
            print OUT $header."\n";
        }

    }

    while ( my $a = $vcf->next_var() ) {

        my @str = split/\t/,$a->{vcf_str};
        my $start = $a -> {POS};

        ## Filter variants in introns and invalid positions
        #if ($a->{INFO}->{INTERNAL} =~ /INTRONIC|null|PROMOTER/) { next; }
        ## Filter variants due to low quality for sample
        if ($a->{FILTER} =~ /ac0/) { next; }

        ## Print first 6 columns of vcf
        print OUT join("\t",@str[0..5])."\t";
        ## All variants PASS, anything above is filtered out
        ## anyhting below is included
        print OUT "PASS\t";
 

        my @alt = split/:/,$a->{ALT};
        my $alt = $alt[2];
        $alt =~ s/>//g;
        my $event = $a->{INFO}->{INTERNAL};
        $event =~ s/\,/_/g;
        
        my $meinfo = $a->{INFO}->{MEINFO};
        $meinfo =~ s/\,/ /g;

        my @scoutcustom;
        my @newinfo;
        my @qcf = (split ';', $a->{FILTER});
        my $qcf = '';
        foreach my $s (@qcf) {
            if ($s eq 'lc') {
                $qcf = 'lc';
            }
            elsif ($s eq 'PASS') {
                $qcf = 'PASS';
            }
            else {
                $qcf = $qcf[0];
            }
        }
        push @newinfo,"SVTYPE=INS";
        push @newinfo,"END=$start";
        push @scoutcustom,"SCOUT_CUSTOM=Repeat Element|".$alt;
        push @scoutcustom,"Subtype|".$meinfo;
        push @scoutcustom,"Consequence|".$event;
        push @scoutcustom,"Target Site Duplication|".$a->{INFO}->{TSD};
        push @scoutcustom,"Filter|".$qcval{$qcf};
        push @scoutcustom,"Assess|".$assess{$a->{INFO}->{ASSESS}};
        print OUT join(';',@newinfo).";";
        print OUT "MELT_RANK=".$a->{INFO}->{ASSESS}.";";
        print OUT "MELT_QC=".$qcf.";";
        print OUT join(',',@scoutcustom)."\t";

        print OUT "GT"."\t".$a->{GT}->[0]->{GT};
        #print OUT join("\t",@str[8..$#str]);
        print OUT "\n";

    }

    close OUT;
    push @non_empty,$out;

}

## If no variants in either of three VCFs, print empty VCF with only header
if (scalar(@non_empty) == 0) { 
    my $command = `cat $vcfheader > $id\.melt.merged.vcf`;
    my $sed = `sed -i 's/FORMAT/FORMAT\t$id/' $id\.melt.merged.vcf`;
}
## Concat all with values
else {
    my $vcfs = join(' ',@non_empty);
    my $command = `vcf-concat $vcfs | vcf-sort -c > $id\.melt.merged.vcf`;
}
