#!/usr/bin/perl -w
##############################################################################
#                                       
#   The ab-initio simulation package  S F H I n g X      
#                                     http://www.sfhingx.de
# 
#            by Sixten Boeck et al.   boeck@sfhingx.de
# 
#############################################################################

use strict;

# --- prototypes
sub parseAndOverlayX($$$@);
sub error ($$$);

sub checkLeadingTabs ($$);
sub checkLineWidth   ($$);
sub checkAssignments ($$);
sub checkDerefer     ($$);

# --- configuration
my $exitOnError = 0;


# --- Check commine line
if ($#ARGV != 0 ) {
   print "Usage: sxstylecheck file.cpp\n";
   exit;
}
my $file = $ARGV[0];

# --- Check whether file exists
if ( ! -e $file )  {
   die "File $file does not exist.\n";
}
my $exitCode = 0;
my $nErrors = 0;

# --- Parse file line by line
my (@lines, @origLines);
open FP, "<$file" || die "Can't open file $file.\n";
while (<FP>)  { push @lines, $_; }
close FP;
@origLines = @lines;  # keep a copy for better error message printout.

# --- Replace C/C++ comments as well as strings with placeholders
@lines = parseAndOverlayX ('\/\*', '\*\/', '!', @lines);   # /* ... */
@lines = parseAndOverlayX ('\/\/', '$',    '@', @lines);   # // comment
@lines = parseAndOverlayX ('"', '"', ,     '#', @lines);   # "strings"

# --- Perform style checks
my $i;
my $lineNr = 1;
foreach $i (@lines)  { 
   checkLeadingTabs ($i, $lineNr);
   checkLineWidth   ($i, $lineNr);
   checkAssignments ($i, $lineNr);
   checkDerefer     ($i, $lineNr);
   $lineNr++;
}


# --- print number of all errors in file
print "File $file contains $nErrors style guide violation(s).\n";

exit 0;
#exit $exitCode;




# ---------------------------------------------------------------------------
# Parse file line by line and replace comments/strings with placeholders.
# Examples for start end end tokens:
#    '\/\*'   '\*\/'      -> C-like comment:     /* ... */    
#    '\/\.'   '$'         -> C++-like comment:   // comment
#    '"'      '"'         -> string:             "hello world"
sub parseAndOverlayX($$$@)
{
   my ($startToken, $endToken, $placeholder, @lines) = @_;
   my @res;

   my ($i, $l, $line, $left, $rightPos);
   my $overlayX = 0;
   foreach $l (@lines)  {
      $_ = $l;
      if ($overlayX)  {
         if (/$endToken/) { # /* ... */ comment block, end in new line
            $line = "";
            for ($i=0; $i < length($`)+length($&); $i++) { 
               $line .= $placeholder; 
            }
            $line .= $';
            $overlayX = 0;
         }  else  {    # eat up entire line
            $line = "";
            for ($i=0; $i < length($l)-1; $i++) { $line .= $placeholder; }
            $line .= "\n";
         }
      }  else  {
         $line = $l;
         while ( /$startToken/ )  {        # /* ... */ comment block, begin
            $overlayX = 1;
            $line = $`;
            for ($i=0; $i < length($&); $i++) { $line .= $placeholder; }
            if ( $' =~ /$endToken/ )  { # /* ... */ block end in same line
               $overlayX = 0;
               $rightPos = length($`) + length($&);
               for ($i=0; $i < $rightPos; $i++) { $line .= $placeholder; }
               $line .= $';
            } else {                # /* ... begin of multiline comment block
               for ($i=0; $i < length($')-1; $i++) { $line .= $placeholder; }
               $line .= "\n";
            }
            $_ = $line;
         } 
      }
      push @res, $line;
   }
   return @res;
}




# ---------------------------------------------------------------------------

# Prompt an error message
sub error ($$$)
{
   my $lineNr = shift;
   my $column = shift;
   my $mesg   = shift;
   my $output = "$file:$lineNr: " . $mesg . "\n";
   my $c;

   # --- print error message
   $exitCode = 1;
   if ($exitOnError)  {
      die $output;
   }  else  {
      print $output;
   }

   # --- print error position indicator ('^')
   print $origLines[$lineNr-1];
   for ($c=0; $c < $column; $c++)  {  print " "; }
   print "^\n";

   print "---\n";

   $nErrors++;
 }



# ---------------------------------------------------------------------------

# Rule: Leading tabulators are not allowed.
sub checkLeadingTabs ($$)
{
   my ($line, $number) = @_;
   if ( $line =~ /^(\s*\t)/ )  {
      error ($number, length $1,
             "Leading tabulators are not allowed.");
   }
   
}



# Rule: The maximum width is restricted to 80 columns
sub checkLineWidth ($$)
{
   my ($line, $number) = @_;
   # --- replace tabs width 3 spaces
   while ($line =~ s/\t/\ \ \ /g) {  }

   my $size = length ($line) - 2;
   if ($size > 80)  {
      error ($number, 80,
             "The maximum width is restricted to 80 columns (width=$size).");
   }
}

# Rule: Assignments are aligned. Put spaces between operators
sub checkAssignments ($$)
{
   my ($line, $number) = @_;

   if ( $line =~ /^(\s*\w+=)/  || $line =~ /(^\s*\w+\s*=)\w+/)  {
      error ($number, length($1)-1,
             "Assignments are aligned. Put spaces between operators.");
   }
}

# Rule: Addressing members of objects works without whitespaces
sub checkDerefer ($$) 
{
   my ($line, $number) = @_;
   $_ = $line;
   if (  /(^.*\w+\s+[\.])\s+\w+/ 
      || /(^.*\w+[\.])\s+\w+/  
      || /(^.*\w+\s+[\.])\w+/ )  
   {
      error ($number, length($1)-1,
             "Dereferring static members only without whitespaces.");
   }

   if (  /(^.*\w+\s+->)\s+\w+/ 
      || /(^.*\w+->)\s+\w+/  
      || /(^.*\w+\s+->)\w+/ )  
   {
      error ($number, length($1)-2,
             "Pointer to members only without whitespaces.");
   }

}



