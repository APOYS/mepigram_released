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
	$seq = uc $seq; #upper case
	my $length = length $seq;
	for(my $i = 1; $i <= ($length - $WORDLENGTH); $i++)
		{
		my $word = substr $seq, $i, $WORDLENGTH;
		if($word !~ /[NRM]+/) #doesnt contain N or R or M
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

open OUT, ">background_typeE-$WORDLENGTH.tsv" or die "can't write to background file - $!\n";
print OUT "TOTAL\t$total\n";
for my $eightMer (keys %nt2counts)
	{
	print OUT "$eightMer\t$nt2counts{$eightMer}\n";
	}
close OUT;

sub revComp
{
my $revComp = shift;
$revComp =revComp2($revComp);
return $revComp;
}
sub revComp2
{
	my $word = shift;
	my $revComp = '';
	my %transdict = ('A' => 'T','T' =>'A','C'=>'G','G'=>'C','N'=>'N');
	my $index= 0;
	while ($index <length($word)){
		my $char = substr($word, $index, 1);
		if ($char eq 'E'){
			if (substr($word, $index +1, 1) eq 'G'){
				$revComp.='GE';
				$index+=2;
			
			}
			else {
				$revComp.='G';
				$index++;
			}
		}
		else{
			
		$revComp.=$transdict{$char};
		$index++;
		}
	}

	$revComp = reverse $revComp;
	return $revComp;
}


