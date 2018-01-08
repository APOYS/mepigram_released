use strict;
use warnings;

my $WORDLENGTH = shift;
my $GENOMEDIR = shift;

if(! defined $WORDLENGTH or ! defined $GENOMEDIR)
	{
	print "You must specifiy a word length (e.g. 8) and a genome dir (a dir containing only genome fasta sequnce files - idealy one chr per file)\n";
	exit;
	}

my (%nt2counts, $total);
my @genomeFiles = glob "$GENOMEDIR/*";
FILE: for my $file (@genomeFiles)
	{
	print "Reading $file\n";
	my $seq;
	open IN, "<$file" or die "can't open $file - $!\n";
	while(<IN>)
		{
		unless(/^>/)
			{
			chomp;
			$seq .= $_;
			}
		}
	close IN;
	$seq = uc $seq;
	my $length = length $seq;
	for(my $i = 1; $i <= ($length - $WORDLENGTH); $i++)
		{
		my $word = substr $seq, $i, $WORDLENGTH;
		if($word !~ /[NRM]+/)
			{
			$nt2counts{$word}++;
			$total++;
			}
		}
	$seq = revComp($seq);
	for(my $i = 1; $i <= ($length - $WORDLENGTH); $i++)
		{
		my $word = substr $seq, $i, $WORDLENGTH;
		if($word !~ /[NRM]+/)
			{
			$nt2counts{$word}++;
			$total++;
			}
		}
	}

open OUT, ">background_typeEF-$WORDLENGTH.tsv" or die "can't write to background file - $!\n";
print OUT "TOTAL\t$total\n";
for my $eightMer (keys %nt2counts)
	{
	print OUT "$eightMer\t$nt2counts{$eightMer}\n";
	}
close OUT;

sub revComp
{
my $revComp = shift;
$revComp =~ tr/ACGTEF/TGCAFE/;
$revComp = reverse $revComp;
return $revComp;
}

