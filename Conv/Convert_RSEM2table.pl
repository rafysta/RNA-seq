#!/usr/bin/perl
# Jul 3, 2017 ver.1

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV != 8 or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [RSEM file] -m [sample mapping] -o [output file] -a [annotation file]\n";
}
my %opt;
getopts("i:m:o:a:", \%opt);

my $FILE_RSEM = $opt{i};
my $FILE_mapping = $opt{m};
my $FILE_out = $opt{o};
my $FILE_annotation = $opt{a};

#---------------------------------------
# mapping
#---------------------------------------
my %num2name;
my $fh_map = IO::File->new($FILE_mapping) or die "cannot open $FILE_mapping: $!";
while($_ = $fh_map->getline()){
	s/\r?\n//;
	my ($num, $name) = split /\t/;
	$num2name{$num} = $name;
}
$fh_map->close();


#---------------------------------------
# gene annotation
# Transcript	chr	strand	Geene	Symbol
#---------------------------------------
my %title_anno = (
	'Transcript' => 0,
	'chr' => 1,
	'strand' => 2,
	'Gene' => 3,
	'Symbol' => 4
);


my %geneAno;
my $fh_annotation = IO::File->new($FILE_annotation) or die "cannot open $FILE_annotation: $!";
$fh_annotation->getline();
while($_ = $fh_annotation->getline()){
	s/\r?\n//;
	my @data = split /\t/;
	$geneAno{$data[$title_anno{'Transcript'}]} = \@data;
}
$fh_annotation->close();


sub SplitTitle{
	my ($char) = @_;
	$char =~ s/\r?\n//;
	my @data = split /\t/, $char;
	return @data;
}

#---------------------------------------
# read RSEM output
#---------------------------------------
my $fh_RSEM = IO::File->new($FILE_RSEM) or die "cannot open $FILE_RSEM: $!";
my @titles = &SplitTitle($fh_RSEM->getline());
shift @titles;
my @titles_new = ('Gene');

foreach my $t(@titles){
	my $barcode = '';
	my $sample = 0;
	if($t =~ m/(\w+)_Sample_(.+)/){
		$barcode = $1;
		$sample = $2;
	}
	if(exists $num2name{$sample}){
		my $name = $num2name{$sample};
		my $name_w_num = $name . ':' . $barcode;
		push @titles_new, $name_w_num;
	}else{
		print "$t is unknown sample\n";
	}
}

my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";
$fh_out->print(join("\t", @titles_new) . "\n");

while($_ = $fh_RSEM->getline()){
	s/\r?\n//;
	my ($transcript, @data) = split /\t/;
	
	my @transcripts = split /,/, $transcript;
	my %GeneList;
	foreach my $t(@transcripts){
#		$t =~ s/\.\d+//;
		if(exists $geneAno{$t}){
			$GeneList{$geneAno{$t}->[$title_anno{'Gene'}]}++;
		}
	}
	
	my $gene = "NA";
	foreach my $g(sort {$GeneList{$b} <=> $GeneList{$a}} keys %GeneList){
		$gene = $g;
		last;
	}
	if($gene eq "NA"){
		print "no gene name found for $transcript\n";
	}else{
		$fh_out->print(join("\t", $gene,  @data) . "\n");
	}
}
$fh_RSEM->close();
$fh_out->close();












