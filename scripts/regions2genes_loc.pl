#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename; # used for retrieving script name (see "fileparse")
use List::Util qw(sum);


#use lib '/data/common_modules/ensembl69/ensembl/modules';

#use lib '/ccagc/home/stef/code/ENSEMBL_API/ensembl74/ensembl/modules/';
use lib '/net/nfs/PAT/home/stef/code/ENSEMBL_API/ensembl74/ensembl/modules/';
#use lib '/ccagc/home/stef/code/ENSEMBL_API/bioperl-1.2.3/';
#use lib '/net/nfs/PAT/home/stef/code/ENSEMBL_API/bioperl-1.2.3/';

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

## when run with option -coding these transcript biotypes are skipped
my $SKIP_REGEX = 'pseudogene|retained_intron|nonsense_mediated_decay|processed_transcript|ambiguous_orf';
#my @CHRS = qw( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M MT Mt mt x y );
my $SCRIPT_NAME = (fileparse($0))[0];

my %opt = (
  'help'          => undef,
  'bed'           => undef,
  'species'       => 'Homo_sapiens',
  'loc_serv'      => 'zeus',
  #'web_serv'      => 'ensembldb.ensembl.org', # alt: useastdb.ensembl.org
  'verbose'       => undef,
  'no_local'      => undef,
  'prot_coding'   => undef,
  'debug'         => undef,
);

my $usage = <<END;

  For any number of bed file(s):
    - adds all genes in those regions to extra column

  Example usage:
    $SCRIPT_NAME -o my_analysis -bed regions.bed

   Required params:
     -out_base|o     output file name
     -bed            bed file (chr<TAB>start<TAB>end)

   Debug options:
     -verbose     load ensembl registry with verbose = 1

   Example input file:
   ------------------------------------
     # commented lines are ignored
     chr1 54321 64321 column4

   Example output file:
   ------------------------------------
     # commented lines are printed back
     chr1 54321 64321 column4 GENE1,GENE2

END


# ======================================================
# Options / init
# ======================================================
die $usage if @ARGV == 0;
GetOptions (
  'h|help'   => \$opt{help},
  'bed|b=s@' => \$opt{bed},
  'out|o=s'  => \$opt{out},
  'log|l=s'  => \$opt{log},
  'debug'    => \$opt{debug},
  'prot|p'   => \$opt{ prot_coding },
  'serv=s'   => \$opt{ web_serv },
)
or die $usage;
die $usage if $opt{help};
die "[ERROR] Missing input: pls provide output base (-o)\n" unless $opt{ out };
die "[ERROR] Missing input: pls provide bed file (-bed)\n" unless $opt{ bed };
die "[ERROR] Missing input: pls provide log file (-log)\n" unless $opt{ log };

my $registry = 'Bio::EnsEMBL::Registry';
## create slice-adaptor and gene-adaptor
my ( $sa, $ga, $ga_of, $dbentry_a ); # of = otherfeatures (is needed for refseq ids NM_ input, so use no_local)

## first try local db
loadEnsemblAdaptors( $opt{loc_serv}, 'anonymous', $registry, \$sa, \$ga, \$ga_of, \$dbentry_a ) unless $opt{ no_local };
## if local failed, try online
#loadEnsemblAdaptors( $opt{web_serv}, 'anonymous', $registry, \$sa, \$ga, \$ga_of ) unless $sa and $ga and $ga_of;
#loadEnsemblAdaptors( $opt{web_serv}, 'anonymous', $registry, \$sa, \$ga, \$ga_of, \$dbentry_a ) unless $sa and $ga;
## only continue if adaptor loaded
#die "[ERROR] (some) ensembl adaptors not loaded...\n" unless $sa and $ga and $ga_of;
die "[ERROR] (some) ensembl adaptors not loaded...\n" unless $sa and $ga;

# ======================================================
## start analysis
# ======================================================
my $out_file = $opt{out}.'.geneAnn.bed';
my $log_file = $opt{log}; # $opt{out}.'geneAnn.log'
open my $fh_out, ">$out_file" or die "$!: @_\n";
open my $fh_log, ">$log_file" or die "$!: @_\n";

printMsg( 'i', "...OK adaptors loaded (ensembl_version: ".software_version().")", *STDOUT, $fh_log );

my $genes_in   = 'INPUT_BED_FILE='.join(',', @{$opt{ bed }} );
my $script_n   = 'SCRIPT_NAME='.(fileparse($0))[0];
my $ensembl_v  = 'ENSEMBL_VERSION='.software_version();
my $host_n     = 'HOST='.$sa->dbc->host;
my $db_n       = 'DATABASE='.$sa->dbc->dbname;
my $species_s  = 'SPECIES='.$opt{species};
my $analysis_n = 'ANALYSIS_NAME='.$opt{ out };
my $settings   = '##ANNOTATOR_SETTINGS: '.join(';', $script_n, $host_n, $ensembl_v, $species_s, $genes_in )."\n";
#my @comment_lines = ($design_n, $script_n, $ensembl_v, $host_n, $species_s, $genes_in );
my @comment_lines = ($analysis_n, $script_n, $ensembl_v, $host_n, $db_n, $genes_in );

## print info block to each filehandle output
foreach ( $fh_out, $fh_log, *STDOUT ){
    print $_ join( "\n", map( '## '.$_, @comment_lines) )."\n";
}

## load genes to include exons from
my $total_genes_printed = 0;
foreach my $f ( @{$opt{bed}} ){
    printMsg( 'i', "Reading input file [$f]...", *STDOUT );
    #parseGeneFile( \%ids_to_include, $f, \@id_array );
    open IN, '<', $f or die "Unable to open infile [$f]\n";
    while (<IN>){
    	chomp;
    	next if ( $_ =~ /^#|chrom|Chrom/ ); # comments/header
    	next if ( $opt{ debug } and $. > 5 );

    	my ($chr, $start, $end) = (split("\t", $_))[0,1,2];
    	$chr =~ s/chr//;
    	my $region_string = $chr.':'.$start.'-'.$end;
    	print "[INFO] Reading for region: ".$region_string."\n";
    	#my $slice = $sa->fetch_by_location( $region_string,'chromosome');
    	my $slice = $sa->fetch_by_region( 'chromosome', $chr, $start, $end );

		my @final_genes = ();
		foreach my $gene ( @{ $slice->get_all_Genes } ) {
		    my $gene_name = $gene->external_name;
		    my $gene_biotype = $gene->biotype();
		    #my $gene_id = $gene->stable_id if $gene;

	    	#my $chr  = $gene->seq_region_name if $gene;
	    	#my $synonym_string = getGeneSynonyms( $gene );
	    	next if ( $opt{ prot_coding } and ($gene_biotype ne 'protein_coding') );
	    	push( @final_genes, $gene_name );
		}
		my $genes_string = join( ",", @final_genes );
		$total_genes_printed += scalar(@final_genes);
		print $fh_out join( "\t", $_, $genes_string)."\n";
    }
    close IN;
}

printMsg( 'i', "Total of $total_genes_printed genes printed...", *STDOUT );

close $fh_out;
close $fh_log;
printMsg( 'i', "DONE", *STDOUT );



# ======================================================
## SUBROUTINES
# ======================================================
sub printMsg{
    my ($type, $msg, @fhs) = @_;
    my $pre = '';
    if   ( $type eq 'w' ){ $pre = "[---WARNING---] "; }
    elsif( $type eq 'e' ){ $pre = "[----ERROR----] "; }
    elsif( $type eq 'i' ){ $pre = "[INFO] "; }
    else{ die "[printMsg] wrong type input\n"; }
    foreach ( @fhs ){
	print $_ $pre.$msg."\n";
    }
}

sub printHashInfo{
    my ($hash, $title) = @_;
    warn "--- $title ---\n";
    while ( my ($key,$val) = each( %{$hash} ) ){
	warn join("\t", $key, $val)."\n";
    }
}

sub getGeneSynonyms{
    my ($gene) = @_;
    my $dbes = $gene->get_all_DBEntries('HGNC');
    foreach my $dbe ( @{$dbes} ) {
        if ( $dbe->dbname() eq 'HGNC' ){
	    return( join(',', @{$dbe->get_all_synonyms}) ) if scalar( @{$dbe->get_all_synonyms} );
        }
    }
    return( '.' );
}

sub arrayMean{
    my $array = shift;
    my $mean = sum(@$array)/@$array;
    return( sprintf("%.2f", $mean) );
}

sub loadEnsemblAdaptors{
    my ($serv, $user, $reg, $sa, $ga, $ga_of, $da) = @_;
    printMsg( 'i', "Will now try to connect to ensembl db at [$serv]", *STDOUT );
    my $species = 'Homo_sapiens';
    $reg->load_registry_from_db( -host => $serv, -user => $user, -verbose => $opt{verbose}, -species => $species );
    $$ga  = $reg->get_adaptor( $species, 'core', 'Gene');
    $$sa  = $reg->get_adaptor( $species, 'core', 'slice' );
    $$ga_of  = $reg->get_adaptor( $species, 'otherfeatures', 'Gene');
    $$da  = $reg->get_adaptor( $species, 'core', 'dbentry');
}
