#include <SxTree.h>
#include <SxCLI.h>

int main (int argc, char **argv)
{
   SxCLI cli (argc, argv);
   cli.finalize ();

   SxTree<int> tree;

   // insert values
   tree << 1 << 2 << 3;
   // insert sub trees
   SxTree<int> subTree10, subTree100;
   subTree10 << 10 << 20 << 30;
   subTree100 << 100 << 200 << 300;
   
   tree << subTree10 << subTree100;

   tree.print ();
   return 0;
}
