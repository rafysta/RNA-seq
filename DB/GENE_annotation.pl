#!/usr/bin/perl
# 2017/03/22 Gene fileを作成する

use strict;
use warnings;
use IO::File;
use Getopt::Std;
use Carp qw(croak);
$| = 0;

if(@ARGV != 6 or $ARGV[0] eq '--help'){
	die "Usage : $0 -g [Gene for IPA] -m [Mapped ID by IPA] -x [organism]\n";
}
my %opt;
getopts("g:m:x:", \%opt);
my $FILE_gene = $opt{g};
my $FILE_ipa= $opt{m};
my $ORGANISM = $opt{x};

my $FILE_mapping;
if($ORGANISM eq 'human'){
	$FILE_mapping = '/wistar/noma/Data/Human_seq/hg19/hg19_transcriptID_geneID_mapping.txt';
}elsif($ORGANISM eq 'mouse'){
	$FILE_mapping = '/wistar/noma/Data/Mouse_seq/mm10/mm10_transcriptID_geneID_mapping.txt';
}


# Transcript and Gene ID connection
my %Trasnscripts;
{
	my $fh_map = IO::File->new($FILE_mapping) or die "cannot open $FILE_mapping: $!";
	$fh_map -> getline();
	while($_ = $fh_map->getline()){
		s/\r?\n//;
		my ($gene, $transcript) = split /\t/;
		my @split_gene = split /_/, $gene;
		$gene = $split_gene[0];
		push @{$Trasnscripts{$gene}}, $transcript;
	}
	$fh_map->close();
}


# Symbol and Gene ID connection
my %Symbol;
{
	my $fh_ipa = IO::File->new($FILE_ipa) or die "cannot open $FILE_ipa: $!";
	while($_ = $fh_ipa->getline()){
		s/\r?\n//;
		unless(m/^E/){
			next;
		}
		my ($gene, $symbol, $description) = split /\t/;
		$Symbol{$gene} = [$symbol, $description];
	}
	$fh_ipa->close();
}



my $fh_gene = IO::File->new($FILE_gene) or die "cannot open $FILE_gene: $!";
my $title = $fh_gene->getline();
$title =~ s/\r?\n//;
my @titles = split /\t/, $title;
shift @titles;
print "GeneID\tTranscriptID\tDescription\tSymbol\t" . join("\t", @titles) . "\n";
while($_ = $fh_gene->getline()){
	s/\r?\n//;
	my ($gene, @other) = split /\t/;
	my $transcript = join(",", @{$Trasnscripts{$gene}});
	unless(exists $Symbol{$gene}){
		next;
	}
	my ($symbol, $description) = @{$Symbol{$gene}};

	my @output;
	push @output, $gene, $transcript, $description, $symbol, @other;
	print join("\t", @output) . "\n";
}
$fh_gene->close();



