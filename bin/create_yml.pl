#! /usr/bin/perl -w
#use MongoDB;
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use JSON;
use List::MoreUtils qw(uniq);

my %opt = ();
GetOptions(
    \%opt,
    'g=s',
    'd=s',
    'p=s',
    'out=s',
    'genome=s',
    'antype=s',
    'ped=s',
    'assay=s',
    'files=s',
    'panelsdef=s',
    'extra_panels=s'
) or usage();

my @required = ('g', 'd', 'out', 'files', 'ped', 'panelsdef', 'assay', 'antype');

my @missing;
foreach my $param (@required) {
    push @missing, $param unless exists $opt{$param};
}

if (@missing) {
    print STDERR "Missing required options: " . join(", ", @missing) . "\n";
    exit 1
}

### Define scout institute per assay and/or analysis. ###
my %assays = (
    'oncov1-0' => {
        'rankm' => 'OncoV1.0-SNV-RM-v5.0',
        'svrankm' => 'OncoV1.0-SV-Panel-RM-v1.0',
        'capture_kit' => 'oncov1-0',
        'screening' => {
            'institute' => 'onkogenetik',
            'institute_owner' => 'onkogenetik',
        },
        'predictive' => {
            'institute' => 'onkogenetik',
            'institute_owner' => 'onkogenetik'
        },
        'PARP_inhib_normals' => {
            'institute' => 'onkogenetik',
            'institute_owner' => 'onkogenetik'
        }
    },
    'oncov2-0' => {
        'rankm' => 'OncoV2.0-SNV-RM-v5.0',
        'svrankm' => 'OncoV2.0-SV-Panel-RM-v1.0.2',
        'capture_kit' => 'oncov2-0',
        'screening' => {
            'institute' => 'onkogenetik',
            'institute_owner' => 'onkogenetik',
        },
        'predictive' => {
            'institute' => 'onkogenetik',
            'institute_owner' => 'onkogenetik'
        },
        'PARP_inhib_normals' => {
            'institute' => 'onkogenetik',
            'institute_owner' => 'onkogenetik'
        }
    },
    'wgs-hg38' => {
        'rankm' => '5.2.1_SNV',
        'svrankm' => '5.2_SV',
        'ph' => {
            'institute' => 'klingen_38',
            'institute_owner' => 'klingen_38'
        },
        'ahus' => {
            'institute' => 'ahus',
            'institute_owner' => 'ahus',
            'rstreshold' => -10000
        },
        'kit' => {
            'institute' => 'KIT',
            'institute_owner' => 'KIT'
        },
        'wgsvalid' => {
            'institute' => 'klingen_38',
            'institute_owner' => 'klingen-genlista'
        }
    },
    'wgs_dev' => {
        'rankm' => '5.2.1_SNV',
        'svrankm' => '5.2_SV',
        'ph' => {
            'institute' => 'klingen-genlista',
            'institute_owner' => 'klingen_38'
        },
        'ahus' => {
            'institute' => 'ahus',
            'institute_owner' => 'ahus'
        },
        'wgsvalid' => {
            'institute' => 'klingen_38',
            'institute_owner' => 'klingen-genlista'
        }
    },
    'dev' => {
        'rankm' => '5.2.1_SNV(dev placeholder)',
        'svrankm' => '5.2_SV(dev placeholder)',
        'wgs' => {
            'institute' => 'klingen_38',
            'institute_owner' => 'klingen-genlista'
        },
        'onco' => {
            'institute' => 'onkogenetik',
            'institute_owner' => 'oncogen'
        },
        'mody' => {
            'institute' => 'MODY',
            'institute_owner' => 'klingen-genlista',
        },
        'cf' => {
            'institute' => 'cf-test',
            'institute_owner' => 'klingen-genlista',
        },
        'ph' => {
            'institute' => 'klingen_38',
            'institute_owner' => 'klingen-genlista'
        },
        'hemato' => {
            'institute' => 'myeloid',
            'institute_owner' => 'klingen-genlista',
        },
        'ahus' => {
            'institute' => 'ahus',
            'institute_owner' => 'klingen-genlista'
        }
    },
    'myeloid-const' => {
        'rankm' => 'SNV-RM-v5.0',
        'svrankm' => 'SV-Panel-RM-v1.0',
        'hemato' => {
            'institute' => 'myeloid',
            'institute_owner' => 'myeloid',
        },
    },
    'mody-cftr-aatv1-0' => {
        'rankm' => 'SNV-RM-v5.0',
        'svrankm' => 'SV-Panel-RM-v1.0',
        'capture_kit' => 'mody-ctfr-aatv1-0',
        'mody' => {
            'institute' => 'MODY',
            'institute_owner' => 'MODY',
        },
        'cf' => {
            'institute' => 'cf-test',
            'institute_owner' => 'cf',
        }
    }
    
);


### Analysis type ###
my $antype = "wgs";
if ($opt{antype}) { $antype = $opt{antype}; }

### Genome build ###
my $genome = "38";
if ($opt{genome}) { $genome = $opt{genome}; }

### Assay ###
my $assay = "wgs";
my $analysis = "";
if ($opt{assay}) { 
    my @a_a = split/,/,$opt{assay};
    $assay = $a_a[0];
    if ($a_a[1]) {
        $analysis = $a_a[1];
    }
    elsif ($opt{d} eq 'ahus') { ## beginning of stinking mess, deadline for fix 2021-03-01
        $analysis = 'ahus';
    }
    else { $analysis = 'ph';}
    
}
print STDERR "assay: $assay analysis: $analysis\n";
### Group ###
### Proband ### Could differ from group, needed to select correct eklipse image
### Clarity-ID ###
my @g_c;
if (defined $opt{g}) { 
    @g_c = split/,/,$opt{g};
    unless (scalar(@g_c) == 2) {
        print STDERR "need group-id,clarity-id\n";
        exit; 
    }
}
else {
    print STDERR "need group-id,clarity-id\n";
    exit;       
}
my $group = $g_c[0];
my $clarity_id = $g_c[1];

### Read ped, save individuals ####################
my $files = $opt{files};
open (INFO, $files) or die "Cannot open $files\n";
my %INFO;
my @bams;
# For mother affected, father affected also print files associated to these
my @inher_patterns;
while ( <INFO> ) {

    my @tmp = split/\s+/,$_;
    my $category = $tmp[0];
    my $subcat = $tmp[1];
    my $filepath = $tmp[2];
    if ($category eq "BAM") {
        $INFO{BAM}->{$subcat} = $filepath;
    }
    elsif ($category eq "TISSUE") {
        $INFO{TISSUE}->{$subcat} = $filepath;
    }
    elsif ($category eq "mtBAM") {
        $INFO{mtBAM}->{$subcat} = $filepath;
    }
    elsif ($category eq "D4") {
        $INFO{D4}->{$subcat} = $filepath;
    }
    elsif ($category eq "IMG") {
        $INFO{IMG}->{$subcat} = $filepath;
    }
    elsif ($category eq "STR_IMG") {
        $INFO{STR_IMG}->{$subcat} = $filepath;
    }
    elsif ($category eq "SV" or $category eq "SVc" or $category eq "SNV" or $category eq "MADDE") {
        if ($category eq "SNV") {
            push @inher_patterns,$subcat;
        }
        $INFO{$category}->{$subcat} = $filepath;
    }
    else {
        $INFO{$category} = $subcat;
    }

}
close INFO;
#my $info_json = to_json(\%INFO, { pretty => 1, indent => 4 });
#print STDERR ($info_json);
####################################################

my $kit = "Intersected WGS"; ## placeholder, does not change for panels
my $diagnosis = $opt{d};

### Read ped, save individuals ####################
my $PED = $opt{ped};
open (PED, $PED) or die "Cannot open $PED\n";
my @ped;
while ( <PED> ) {
    push @ped, $_;
}
close PED;
####################################################

### get genlist ###
my $institute = "klingen";
my $institute_owner = "klingen";
if ($opt{assay}) { 

    unless (exists $assays{$assay}) {
        die "Error: Assay '$assay' not found.\n";
    }

    unless (exists $assays{$assay}{$analysis}) {
        die "Error: Analysis '$analysis' for assay '$assay' not found.\n";
    }

    unless (exists $assays{$assay}{$analysis}{institute}) {
        die "Error: 'institute' key missing or undefined for assay '$assay' and analysis '$analysis'.\n";
    }

    $institute = $assays{$assay}{$analysis}{institute};
    $institute_owner = $assays{$assay}{$analysis}{institute_owner};
    if ($assays{$assay}{capture_kit}) {
        $kit = $assays{$assay}{capture_kit};
    }
}
my $gene_panels = get_genelist($institute);

####################################################


### PRINT YAML PER INHER_PATTERN ####
foreach my $ind (@inher_patterns) {
    my $family = $group;
    my $out = "$group.yaml";
    unless ( $ind eq "proband") {
        $out = "$group.yaml".".".$ind;
        $family = $group."_".$ind;
    }
    
    ### Open out, default $group.yaml, fix for ma and fa! ###
    if ($opt{out}) { 
        $out = $opt{out};
        unless ( $ind eq "proband") {
            $out = $out.".".$ind;
        }
    }
    print("Writing Scout yaml to '$out'\n");
    open (OUT, '>', $out);

    print OUT "---\n";
    ### ASSAY DECIDE OWNER? ####
    print OUT "owner: $institute_owner\n";
    print OUT "family: '$family'\n";
    print OUT "lims_id: '$clarity_id'\n";
    print OUT "samples: \n";

    ### MATCH ped inidividuals with bams ###
    foreach my $sample (@ped) {
        my @pedline = split/\t/,$sample;
        print OUT "  - analysis_type: $antype\n";
        print OUT "    sample_id: '$pedline[1]'\n";
        print OUT "    sample_name: '$pedline[1]'\n";
        print OUT "    mother: '$pedline[3]'\n";
        print OUT "    father: '$pedline[2]'\n";
        print OUT "    capture_kit: $kit\n";

        if ($pedline[5] == 1) {
            if ($ind eq "ma" and $pedline[3] eq "0" and $pedline[2] eq "0" and $pedline[4] == 2) {
                print OUT "    phenotype: affected\n";
            }
            elsif ($ind eq "fa" and $pedline[3] eq "0" and $pedline[2] eq "0" and $pedline[4] == 1) {
                print OUT "    phenotype: affected\n";
            }
            else {
                print OUT "    phenotype: unaffected\n";
            }
            
        }
        elsif ($pedline[5] == 2) {
            print OUT "    phenotype: affected\n";
        }
        else { print STDERR "not a valid phenotype!\n" }
        if ($pedline[4] == 1) {
            print OUT "    sex: male\n";
        }
        elsif ($pedline[4] == 2) {
            print OUT "    sex: female\n";
        }
        else { print STDERR "not a valid sex!\n" }
        if ($INFO{BAM}{$pedline[1]}) {
            print OUT "    bam_path: $INFO{BAM}{$pedline[1]}\n";
        }
        if ($INFO{mtBAM}{$pedline[1]}) {
            print OUT "    mt_bam: $INFO{mtBAM}{$pedline[1]}\n";
        }
        if ($INFO{D4}{$pedline[1]}) {
            print OUT "    d4_file: $INFO{D4}{$pedline[1]}\n";
        }
    }
    ##########################################

    ### Print optional variables ###

    ## If trio and both SV and SNV calling has been done, SVc should be in INFO-file
    ## This contains both SV and SNV info
    if ($INFO{SVc}) {
        my @tmp = split/,/,$INFO{SVc}{$ind};
        print OUT "vcf_snv: $tmp[1]\n"; 
        print OUT "vcf_sv: $tmp[0]\n";
    }
    ## If SNV single, check for SNV, if exist, check for SV
    elsif ($INFO{SNV}) {
        print OUT "vcf_snv: $INFO{SNV}{$ind}\n";
        if ($INFO{SV}) {
            print OUT "vcf_sv: $INFO{SV}{$ind}\n";
        } 
    }
    ## If only SV calling
    elsif ($INFO{SV}) {
        print OUT "vcf_sv: $INFO{SV}{$ind}\n";
    }
    else {
        print STDERR "need at least one VCF, SV/SNV"; exit;
    }
    if ($INFO{STR}) {
        if ($analysis ne 'kit') {
            print OUT "vcf_str: $INFO{STR}\n";
        }
    }
    if ($INFO{SMN}) {
        print OUT "smn_tsv: $INFO{SMN}\n";
    }

    if ($INFO{MADDE}) {
        print OUT "madeline: $INFO{MADDE}{$ind}\n";
    }

    if ($INFO{PEDDY}) {
        my @tmp = split/,/,$INFO{PEDDY};
        print OUT "peddy_ped: $tmp[1]\n";
        print OUT "peddy_check: $tmp[0]\n";
        print OUT "peddy_sex: $tmp[2]\n";
    }

    ### Print gene panels ###
    print OUT "gene_panels: [";
    print OUT join ",", @$gene_panels;
    print OUT "]\n";
    if ($diagnosis eq "pediatrics") {
        print OUT "default_gene_panels: []\n";
    }
    else {
        my @panels = split /\+/, $diagnosis;
        my $panels_str = '"'. join('","', @panels). '"';
        print OUT "default_gene_panels: [$panels_str]\n";
    }
    print OUT "rank_model_version: ".$assays{$assay}{rankm}."\n";
    print OUT "sv_rank_model_version: ".$assays{$assay}{svrankm}."\n";
    if ($assays{$assay}{$analysis}{rstreshold}) {
        print OUT "rank_score_threshold: ". $assays{$assay}{$analysis}{rstreshold}."\n";
    }
    else {
        print OUT "rank_score_threshold: -1\n";
    }
    print OUT "human_genome_build: $genome\n";


    ## If IMGage is available
    my %img = ( 
        'overviewplot' => {
            'desc' => "Genome overview plot, UPD and ROH", 
            'width' => '2000',
            'height' => '1000'
            },
        'eklipse' => {
            'desc' => "Circular mitochondrial plot, Eklipse", 
            'width' => '750',
            'height' => '750'
            },
        'haplogrep' => {             
            'desc' => "Mitochondrial haplotypes, Haplogrep", 
            'width' => '750',
            'height' => '1000'
            },
        'STR' => {             
            'desc' => "Reviewer plot for STR loci", 
            'width' => '500',
            'height' => '100'
            }
    );
    # only print header if there are any images for case, pydantic will crash otherwise
    if ($INFO{IMG} || $INFO{STR_IMG}) {
        print OUT "custom_images:\n";
    }
    if ($INFO{IMG}) {
        print OUT "  case:\n";
        foreach my $img_type (sort keys %{ $INFO{IMG} }) {
            print OUT "    $img_type:\n";
            print OUT "      - title: $INFO{IMG}{$img_type}\n";
            print OUT "        description: $img{$img_type}{desc}\n";
            print OUT "        width: $img{$img_type}{width}\n";
            print OUT "        height: $img{$img_type}{height}\n";
            print OUT "        path: $INFO{IMG}{$img_type}\n";
        }

    }
    if ($INFO{STR_IMG}) {
        print OUT "  str:\n";
        foreach my $img_type (sort keys %{ $INFO{STR_IMG} }) {
            print OUT "    - title: $INFO{STR_IMG}{$img_type}\n";
            print OUT "      description: $img{STR}{desc}\n";
            print OUT "      width: $img{STR}{width}\n";
            print OUT "      height: $img{STR}{height}\n";
            print OUT "      path: $INFO{STR_IMG}{$img_type}\n";
        }
    }

    close OUT;

}



sub get_genelist {
    my $institute = shift;
    
    die "Error: Institute parameter is required\n" unless defined $institute;
    die "Error: 'panelsdef' option is not defined\n" unless exists $opt{panelsdef};


    # Allow outside-institute panels
    my @extra_panels;
    if ($opt{extra_panels}) {
        my $extra_panels_file = $opt{extra_panels};

        my $json_text = '';
        open(my $fh, '<', $extra_panels_file) or die "Cannot open '$extra_panels_file'";
        while (my $line = <$fh>) {
            $json_text .= $line;
        }
        close($fh);

        my $extra_panels_per_institute = decode_json($json_text);

        if (exists $extra_panels_per_institute -> {$institute}) {
            @extra_panels = @{ $extra_panels_per_institute -> {$institute} };
        } else {
            print("No extra panels matching institute '${institute}'\n");
        }
    }

    my $panels_dump = $opt{panelsdef};
    my $json_data;
    my @ok_panels;
    open (JSON, $panels_dump) or die "Cannot open '$panels_dump'";
    while (<JSON>) {
        $json_data = decode_json($_);
    }
    close(JSON);


    foreach my $panel (@{$json_data}) {

        my $panel_name = $panel->{panel_name};

        # Can be a single value or an array reference
        my $panel_institutes_ref = $panel->{institute};

        ## This is because some gene panels are assigned to two different institutes
        ## NOT SUPPORTED by Scout, it's a fluke that it even works
        my @panel_institutes;
        if (ref $panel_institutes_ref eq 'ARRAY') {
            @panel_institutes = @{ $panel_institutes_ref };
        } else {
            @panel_institutes = ( $panel_institutes_ref );
        }

        foreach my $panel_inst (@panel_institutes) {
            my $is_extra_panel = grep { $_ eq $panel_name } @extra_panels;
            if ($panel_inst eq $institute || $is_extra_panel) {

                push @ok_panels, $panel_name;
                # One match is enough, stopping the loop
                last;
            }
        }

    }

    @ok_panels = uniq(@ok_panels);

    # print(Dumper(@ok_panels) . "\n");

    return \@ok_panels;
}
