#!/usr/bin/perl
# 2017/05/19 DESeq2の結果のGENEファイルにGeneの情報を追加する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use File::Basename;
use Carp qw(croak);
$| = 0;


my $FILE_geneInformation = '/wistar/noma/Data/Human_seq/hg38/transcript_ID_and_gene_ID.txt';



if(@ARGV != 4  or $ARGV[0] eq '--help'){
	die "Usage : $0 -i [DESeq2 output] -o [output file]\n";
}
my %opt;
getopts("o:i:", \%opt);
my $FILE_in = $opt{i};
my $FILE_out = $opt{o};


#---------------------------------------
# get gene information
#---------------------------------------
my %Data;
my $fh_database = IO::File->new($FILE_geneInformation) or die "cannot open $FILE_geneInformation: $!";
$fh_database->getline();
while($_ = $fh_database->getline()){
	s/\r?\n//;
	my ($geneID, $transcriptID, $name, $description) = split /\t/;
	$Data{$geneID} = [$name, $description];
}
$fh_database->close();


#---------------------------------------
# read input file
#---------------------------------------
my $fh_out = IO::File->new($FILE_out, 'w') or die "cannot write $FILE_out: $!";
my $fh_in = IO::File->new($FILE_in) or die "cannot open $FILE_in: $!";
$fh_out->print($fh_in->getline());
while($_ = $fh_in->getline()){
	s/\r?\n//;
	my ($ID, @others) = split /\t/;
	my ($geneID, @others2) = split /_/;
	if(exists $Data{$geneID}){
		my ($symbol, $description) = @{$Data{$geneID}};
		$fh_out->print(join("\t", $geneID, $description, $symbol, @others) . "\n");
	}
}
$fh_in->close();
$fh_out->close();






