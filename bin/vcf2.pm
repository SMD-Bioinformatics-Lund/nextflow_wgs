package CMD::vcf2;
use strict;
use Data::Dumper;
use Exporter qw(import);
our @EXPORT_OK = qw( new next_line parse_vcf parse_metainfo parse_variant parse_genotype parse_info ); 


sub new {
    my( $class, @args ) = @_;
    my $self = {@args};
    bless $self, $class;
    return $self->_open();
}

sub _open {
    my( $self ) = @_;
    my $fn = $$self{file};

    if( is_gzipped($fn) ) {
        open( $$self{fh}, "zcat $fn |" ) or die "Could not open file $fn";
    }
    else {
        open( $$self{fh}, $fn ) or die "Could not open file $fn";;
    }

    my @head ;
    my $full_header_str;

    while( readline($$self{fh}) ) {
    chomp;

    # Skip empty lines
    next if /^\s*$/;

    $full_header_str .= $_."\n" if /^#/;
    
    # Header row with meta data.
    if( /^##/ ) {
        my( $type, $meta, $val ) = parse_metainfo( $_ );
        if( $type ne "NONE" ) {
        $$self{meta}->{$type}->{$meta->{ID}} = $meta;
        }
        else {
        $$self{meta}->{$meta} = $val;
        }
    }

    # Header with column description.
    elsif( /^#/ ) {
        $_ =~ s/^#//;
        @head = split /\t/;
        last;
    }
    }

    @{$$self{head}} = @head;
    $$self{samples} = splice @head, 9;

    $$self{header_str} = $full_header_str;

    return $self;
}

sub next_var {
    my ($self) = @_;
    
    my $line = readline($$self{fh});
    if ($line) {
    chomp $line;

    my $variant = parse_variant( $line, $$self{head}, $$self{meta} );

    if( $variant->{CHROM} ) {
        return $variant if defined $variant;
    }
    }
    return 0;
}

# Parse VCF file and return
sub parse_vcf {
    
    my $fn = shift;
    my %vcf_meta;
    my @vcf_data;
    my @head ;
    my $full_header_str;

    my $vcf_fh;
    if( is_gzipped($fn) ) {
#    if( $fn =~ /\.gz/ ) {
    open( $vcf_fh, "zcat $fn |" );
    }
    else {
    open( $vcf_fh, $fn );
    }
    
    while( <$vcf_fh> ) {
    chomp;

    # Skip empty lines
    next if /^\s*$/;

    $full_header_str .= $_."\n" if /^#/;
    
    # Header row with meta data.
    if( /^##/ ) {
        my( $type, $meta ) = parse_metainfo( $_ );
        $vcf_meta{$type}->{$meta->{ID}} = $meta if defined $type;
    }

    # Header with column description.
    elsif( /^#/ ) {
        #print STDERR $_;
        $_ =~ s/^#//;
        @head = split /\t/;
    }

    # Actual variant data
    else {
        die "Malformed VCF: No column description header." unless @head;
        my $variant = parse_variant( $_, \@head, \%vcf_meta );
        if( $variant->{CHROM} ) {
        push @vcf_data, $variant if defined $variant;
        }
    }
    }

    my @samples = splice @head, 9;
    $vcf_meta{header_str} = $full_header_str;
    
    return \%vcf_meta, \@vcf_data, \@samples;
}




# Parse VCF meta info line (only FORMAT and INFO)
sub parse_metainfo {
    my $comment = shift;
    
    $comment =~ s/^##//;
    my( $type, $data ) = ( $comment =~ /^(.*?)=(.*)$/ );


    if( $type eq "FORMAT" or $type eq "INFO" or $type eq "SAMPLE" or $type eq "FILTER" ) {
    $data = remove_surrounding( $data, '<', '>' );
    my( $pairs, $order ) = keyval( $data, '=', ',' );
    return $type, $pairs, 0;
    } 
    elsif( $type and $data ) {
    return "NONE", $type, $data;
    }
    

    return undef, undef, undef;
}


# Parse VCF variant line
sub parse_variant {
    my( $var_str, $head, $meta ) = @_;
    my @var_data = split /\t/, $var_str;
    my %var;

    $var{ vcf_str } = $var_str;
    
    # First seven fields
    for ( 0..6 ) {
        $var{ $head->[$_] } = $var_data[$_];
    }

    # Eigth field, INFO
    ( $var{ INFO }, $var{INFO_order} ) = parse_info( $var_data[7] );

    # Parse VEP annotation field, if any
    if( $var{ INFO }->{ CSQ } ) {
        $var{ _CSQstr } = $var{INFO}->{CSQ};
        $var{ INFO }->{ CSQ } = parse_VEP_CSQ( $var{INFO}->{CSQ}, $meta->{INFO}->{CSQ} );
    }

    my @FORMAT = split /:/, $var_data[8];
    $var{ FORMAT} = \@FORMAT;
    
    # Genotypes for each sample
    for ( 9 .. (@var_data-1) ) {
        push @{$var{ GT }}, parse_genotype( $var_data[8], $var_data[$_], $head->[$_] );
    }

    return \%var;
}


# Parse genotype field of VCF
sub parse_genotype {
    my( $format, $data, $sample_id ) = @_;

    my @format = split ':', $format;
    my @data   = split ':', $data;

    my %gt;
    @gt{@format} = @data;
    $gt{_sample_id} = $sample_id;

    return \%gt;
}


# Parse info column of VCF file
sub parse_info {
    my $str = shift;
    my( $info, $order ) = keyval( $str, "=", ";" );

    return $info, $order;
}


sub parse_VEP_CSQ {
    my( $CSQ_var, $CSQ_meta ) = @_;

    $CSQ_meta->{Description} =~ /Consequence annotations from Ensembl VEP\. Format: (.*?)$/;

    my @field_names = split /\|/, $1;
    
    my @transcripts = split /,/, $CSQ_var;

    my @data_transcripts;
    foreach my $transcript_CSQ ( @transcripts ) {
    my @values = split /\|/, $transcript_CSQ;

    my %data;
    for( 0 .. $#field_names ) {
        if( $field_names[$_] eq "Consequence" ) {
        my @conseq_array = split '&', $values[$_];
        $data{ $field_names[$_] } = \@conseq_array;
        }
        else {
        $data{ $field_names[$_] } = ( $values[$_] or "" );
        }
    }

    push( @data_transcripts, \%data )
    }
    return \@data_transcripts;
}

# Removes character(s) defined in arg2 if first in string, and arg3 if last in string.
sub remove_surrounding {
    my( $str, $before, $after ) = @_;
    $str =~ s/^$before//;
    $str =~ s/$after$//;
    return $str;
}


# Parse string with key value pairs. Return hash. 
#  * Keys and values separated by 2nd argument. 
#  * Pairs separated by 3rd argument
#  * Handles commas in values if surrounded by double quotes
sub keyval {
    my( $str, $keyval_sep, $pair_sep ) = @_;

    my @pair_str = split /$pair_sep/, $str;
    my %pairs;
    my @order;
    foreach( @pair_str ) {

    # If key-value separator exists, save the value for the key
    if( /$keyval_sep/ ) {
        my( $key, $val ) = split /$keyval_sep/;
        $pairs{$key} = $val;
        push @order, $key;
    }

    # Otherwise treat the whole string as a flag and set it to one (true).
    else {
        $pairs{$_} = 1;
    }
    }
    return \%pairs, \@order;
}

sub is_gzipped {
    my $fn = shift;

    my $file_str = `file -L $fn`;
    return 1 if $file_str =~ /gzip compressed/;
    return 0;
}


1;
