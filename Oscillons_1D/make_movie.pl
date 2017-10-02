#!/usr/bin/perl -w
#######################################################################
# make_movie
#######################################################################
use Getopt::Std;
$USAGE=<<'EoF';
 Usage: make_movie [-ah] variable_name [ time ]
    Writes gnuplot script to graph all datafiles for variable "variable_name".
    The optional argument "time" specifies the lenght of the pause for each
    graph.  The option -a includes plots for an analytical solution; the option
    -h prints this message.
EoF
#######################################################################
# get option
#######################################################################
my %options=();
getopts("ah", \%options) || die "$USAGE";
die "$USAGE" if defined $options{h}; 
#######################################################################
# count arguments
#######################################################################
$num_args = $#ARGV + 1;
if ($num_args < 1) { die "$USAGE" };
#######################################################################
# check if time is provided as argument 
#######################################################################
if ($num_args == 2) { $time = $ARGV[1]; } 
else { $time = 0.5; } 
#######################################################################
# read variable name
#######################################################################
$variable = $ARGV[0];
#######################################################################
# make_movie
#######################################################################
opendir(DIR,".");
if (open(OUT,">movie_$variable")) { print "Created file movie_$variable...\n"; }
$file_name = $variable . "_dis_";
@files = readdir(DIR);
close(DIR);
$file_number = 0;
while (@files) {
      $file = shift(@files);
      if ($file =~ /$file_name*/ ) {
	  if (defined $options{a}) {
	      print OUT "plot '$file' u 1:2 w l, '$file' u 1:3 w l\npause '$time'\n";
	  } else {
	      print OUT "plot '$file' u 1:3 w l\npause '$time'\n";
	  }
	  $file_number += 1;
      }
}
print "Found $file_number files for variable $variable.\n";
close(OUT);
