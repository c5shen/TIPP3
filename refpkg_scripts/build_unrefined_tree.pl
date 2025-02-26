use File::Spec::Functions qw(rel2abs);
use File::Basename;
use strict;
use File::Path;
use File::Temp qw/ tempfile tempdir /;


my $species_mapping = $ARGV[0];
my $unrefined_tree = $ARGV[1];
my $output = $ARGV[2];

my %map = %{read_mapping($species_mapping, ",")};
rename_tree($unrefined_tree,\%map,$output);


sub trim {
  my $string = $_[0];
  if (not defined($string)) {
      return undef;
  }
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

sub read_mapping {
    my $input_file = $_[0];
    my $delimiter = $_[1];
    
    my %hash = ();
    open(INPUT, "<$input_file");
    my $line = <INPUT>;    
    while ($line = trim(<INPUT>."")) {  
      my @results = split($delimiter, trim($line));
      if (not defined $hash{$results[1]}) {
        $hash{$results[1]} = $results[0];
      } else {
        $hash{$results[1]} .= ",".$results[0];
      }
    }
    close(INPUT);
    return \%hash;
}


sub rename_tree {
  my $tree_file = $_[0];
  my $map = %{$_[1]};
  my $output_file = $_[2];

  open(INPUT, "<$tree_file");
  open(OUTPUT, ">$output_file");
  while (my $result = <INPUT>) {
    $result = trim($result);
    $result =~ s/([,\(])([^(,:\)]+)/$1$map{$2}/g;
    print OUTPUT "$result\n";
  }
  close(INPUT);
  close(OUTPUT);
}
