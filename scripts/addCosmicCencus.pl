#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename; # used for retrieving script name (see "fileparse")

my $SCRIPT_NAME = (fileparse($0))[0];

my %opt = (
  'help'    => undef,
  'file'    => undef,
  'col'     => 4,
  'sep'     => ',',
  'census'  => '/net/nfs/PAT/home/matias/data/ref/cosmic/hg19_okt2015/CosmicMutantExport.tsv',
  'debug'   => undef,
);

my $usage = <<END;

  For every line in input file:
    - reads a certain column (col = $opt{col})
    - splits this string into genes by seperator (sep = $opt{sep})
    - checks presence of each in COCMIC Census list ($opt{census})
    - prints back all plus extra column with census genes

  Example minimal usage:
    -bed input.bed

  Options:
     -col|c     column containing the genes
     -sep|s     character seperating the genes
     -out|o     optional output filename

  Example input file:
  ------------------------------------
     # commented lines are ignored
     chr1 54321 64321 PDE4DIP
     chr2 54321 64321 SOME_GENE
     chr3 54321 64321 SOME_GENE,PDE4DIP,FLI1

  Example output file:
  ------------------------------------
     # commented lines are printed back
     chr1 54321 64321 PDE4DIP PDE4DIP
     chr2 54321 64321 SOME_GENE
     chr3 54321 64321 SOME_GENE,PDE4DIP,FLI1 PDE4DIP,FLI1

END


# ======================================================
# Options / init
# ======================================================
die $usage if @ARGV == 0;
GetOptions (
  'h|help'   => \$opt{help},
  'bed|b=s'  => \$opt{bed},
  'col|c=s'  => \$opt{col},
  'sep|s=s'  => \$opt{sep},
  'out|o=s'  => \$opt{out},
  'census=s' => \$opt{census},
  'debug'    => \$opt{debug},
)
or die $usage;
die $usage if $opt{help};
die "[ERROR] Missing file: census list does not exist? ($opt{census})\n" unless -f $opt{ census };
die "[ERROR] Missing input: pls provide bed file (-bed)\n" unless $opt{ bed };

# ======================================================
## start analysis
# ======================================================

## header line cosmic census file
#Gene Symbol,Name,Entrez GeneId,Chr,Chr Band,Somatic,Germline,Tumour Types(Somatic),Tumour Types(Germline),Cancer Syndrome,Tissue Type,Molecular Genetics,Mutation Types,Translocation Partner,Other Germline Mut,Other Syndrome,Synonyms
my %censusGenes = ();
open CENSUS, '<', $opt{census} or die "Unable to open infile [$opt{census}]\n";
while (<CENSUS>){
  chomp;
  next if ( $_ =~ /^Gene Symbol/ ); # header
  next if ( $opt{ debug } and $. > 5 );
  my ($gene_name) = (split(/\t/, $_))[0];
  $censusGenes{ $gene_name } = 1;
}
close CENSUS;
#print Dumper( \%censusGenes );

my $out_file = $opt{bed};
$out_file =~ s/\.bed/_CosmicCensus.bed/;
$out_file = $opt{out} if defined($opt{out}) and $opt{out} ne '';
open my $fh_out, ">$out_file" or die "$!: @_\n";

## open bed file and read column with genes
## add new column with genes that are present in CENSUS list
open BED, '<', $opt{bed} or die "Unable to open infile [$opt{bed}]\n";
while (<BED>){
  chomp;
	next if ( $_ =~ /^#|chrom|Chrom/ ); # comments/header
	next if ( $opt{ debug } and $. > 5 );
	my ($genes_string) = (split("\t", $_))[$opt{col}-1] || ('');
  my @genes = split( $opt{sep}, $genes_string);
  my @genes_present = ();
  foreach my $gene_name (@genes){
    push( @genes_present, $gene_name ) if defined $censusGenes{ $gene_name };
    #print "YESSSSS\n" if defined $censusGenes{ $gene_name };
  }
  print $fh_out join( "\t", $_, join( $opt{sep}, @genes_present))."\n";
}
close BED;
close $fh_out;
printMsg( 'i', "DONE. Output in $out_file", *STDOUT );


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
