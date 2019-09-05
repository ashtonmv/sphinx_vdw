#!/usr/bin/perl -w

#use diagnostics;
my $vimfile = "$ENV{HOME}/.vim/syntax/sphinx.vim";

if (@ARGV == 0)  {
   printf STDERR <<USAGE;
   This parses S/PHI/nX type definition files, tries to extract keywords
   that are missing in $vimfile and prints
   the vim syntax commands to add the missing keywords.
   The keywords recognized are lines like
   
   NAME { type="..."; ...

   which often occur in std files. The output has to be checked
   and eventually wrapped by hand.

   USAGE:

   $0 [--vimfile vimfile] <stdfiles>

USAGE
   exit;
}

if ($ARGV[0] eq "--vimfile")  {
   shift @ARGV;
   $vimfile = shift @ARGV;
   print STDERR "vimfile = $vimfile\n";
   if ( ! -f "$vimfile" ) {
      print STDERR "$vimfile does not exist!\n";
      $vimfile = undef;
   }
}

my (@names, @group, @flag, @var);

for my $file (@ARGV) {
   if (! -e $file)  {
      printf STDERR "Can't find file '$file'.\n";
      next ;
   }
   open INPUT, 'sed -ne\'s/\s*\([A-Za-z]\+\) *{ *type *= *"\([^"]*\).*/\2 \1/p\' '.$file 
              ."| sort -k1 | uniq |";
   while (<INPUT>) {
      my ($type, $name) = split;
      if ($name ne "" && grep ({/$name/} @names) == 0)  {
         system "grep -q -e '$name ' -e '$name\$' $vimfile" if defined $vimfile;
         if (! defined $vimfile || $? != 0) {
            push @names, $name;
            if ($type eq "group")  {
               push @group, $name;
            } elsif ($type eq "flag")  {
               push @flag, $name;
            } elsif ($type =~ m/real|int|string|file|enum|vector|matrix/)  {
               push @var, $name;
            } else {
               print STDERR "$file: '$name' has unknown type '$type'.\n";
            }
         }
      }
   }
   close INPUT;
   printf "\n\" from $file\n" if (@group || @flag || @var);
   $"=" ";
   printf "syn keyword SxGroup     @group\n" if (@group);
   printf "syn keyword SxAttrib    @flag nextgroup=sxReqSemicolon\n" if (@flag);
   printf "syn keyword SxVariable  @var\n" if (@var);
   @group = @var = @flag = ();
}
print "\n";





