#!/usr/bin/perl -w

# Author: C. Freysoldt, freyso@fhi-berlin.mpg.de
#use warnings;
use strict;
#use diagnostics;

my $nocheck = (@ARGV > 1 && $ARGV[0] eq '--nocheck');
      
shift if $nocheck;

if (@ARGV != 1)  {
   print STDERR <<USAGE;
   This tool prints the tcsh complete specification for the add-on
   given. Currently, this is an experimental version and doesn't
   work for all cases. Notably, arguments are NOT recognized.

   The add-on must understand the --help option, and give standardized
   output (technically speaking: the SxCLI format). Use the --nocheck
   option to suppress the CLI check.

   Usage:

   $0 [--nocheck] <S/PHI/nX add-on>

USAGE
   exit 0;
}

my $addon = shift;

if (! $nocheck)  {
   system "$addon --opts | grep 'SVN Tag' 2>&1 > /dev/null";

   if ($? != 0)  {
      print STDERR "$addon is no SxCLI tool!\n";
      exit 2;
   }
}

open HELP, "$addon --help |";

do { $_ = <HELP>; } until (eof(HELP) || /Usage:/); 

my %options;
my @tokens;
while (!eof(HELP) )  {
   $_ = <HELP>;
   chomp; s/^\s*//;
   push @tokens, (split);
}
   
my $first;

shift @tokens; 

LOOP: {
do {
   $_ = shift @tokens;
   $_ = shift @tokens if ($_ eq "|"); # '|' separating alternative groups
   s/\[//;
   my $forceFlag = /\]/;
   s/\]//;
   s/[{}]//g;
# print STDERR "$_\n";
   if (/</)  {
      # Ignore arguments (no marks)
      while (! />/ ) { die if (@tokens == 0); shift @tokens; };
      next LOOP;
   }
   my @marks = split /\|/;

   if (defined $first)  {
      last LOOP if ($marks[$[] eq $first);
   } else  {
      $first = $marks[$[];
   }

   my $what;
   if ($forceFlag)  {
      $what = "flag";
   } else {
      $_ = shift @tokens;
      # get <...> part 
      if (/</) { while (! />/ ) { $_ .= ' '.shift @tokens; } }
      if (/<(.*)>/)  {
         $what = $1;
      } else {
         unshift @tokens, $_;
         $what = "flag";
      }
   }


   for my $mark (@marks)  {
      $options{$mark} = $what;
   }
      
   
} until (@tokens == 0);
}

close HELP;

my $complete = "complete $addon ";

my (@single, @double);
for (keys %options)  {
   push @single, $_ if /^-[^-]/;
   push @double, $_ if /^--/;
}
$complete .= "'C/--/(@double)/' 'C/-/(@single)/' ";

for my $key (keys %options)  {
   my $what = $options{$key};

   if ($what =~ /file/)  {
      $complete .= "'n/$key/f/' ";
   } elsif ($what eq "flag")  {
      $complete .= "'n/$key/x:$key is a flag./' ";
   } else {
      $complete .= "'n/$key/x:<$what> expected/' ";
   }
}
$complete =~ s/ *$//;

print "$complete\n";


